% -------------------------------------------------------------------------
% Compare PCG normal-equation and LSQR augmented-system solutions of the
% Sharma et al. abdominal QSM objective using identical inputs and weights.
% The original script_abdominal_qsm_ref10.m remains unchanged.
% -------------------------------------------------------------------------

clear
clc

% -------------------------------------------------------------------------
% Paths and acquisition parameters
% -------------------------------------------------------------------------

restoredefaultpath

script_directory = fileparts(mfilename('fullpath'));

% running a section of an unsaved file executes a temporary copy (Editor_xxxxx
% in tempdir), so fall back to the current folder, which should be the repository
if ~isfile(fullfile(script_directory, 'data', 'liver_mask.mat'))
    script_directory = pwd;
end

assert(isfile(fullfile(script_directory, 'data', 'liver_mask.mat')), 'data/ not found in %s: run the script from the repository folder', script_directory)

addpath(genpath(fullfile(script_directory, 'functions')))
addpath(genpath(fullfile(script_directory, 'external', 'NIfTI_20140122')))

addpath(script_directory)

% input_mat = '/autofs/space/daisy_002/users/yuting/20261004_bch/data4QSM/imoco_softgating_TV_PVbin.mat';
% input_mat = '/autofs/space/daisy_002/users/yuting/20261004_bch/PTbin/imoco_softgating_TV_PTbin_nrep1.mat';
input_mat = fullfile(script_directory, 'data', 'exp100926_mgre_PT_imoco_softgating_fw.mat');

orientation_nifti_directory = '/autofs/space/daisy_002/users/yuting/20260605_liver/DICOM/';

output_directory = fullfile(script_directory, 'qsm_output');

TE_seconds = [1.0200 2.4400 3.8600 5.2800 6.7000 8.1200] * 1e-3;
FOV_mm = [350, 250, 350];
B0_strength_t = 2.89;

% Optional MAT file containing a logical variable named fat_reference_roi.
fat_reference_roi_mat = '';


% -------------------------------------------------------------------------
% Load coil sensitivity for masking background
% -------------------------------------------------------------------------

tic
    load(fullfile(script_directory, 'data', 'exp100926_mgre_fb_coilsensmask.mat'))
toc


% -------------------------------------------------------------------------
%% Load and orient the input arrays
% -------------------------------------------------------------------------


fprintf('Loading %s\n', input_mat)
input = load(input_mat);

required_variables = {'img_cplx', 'water', 'fat', 'r2s', 'fieldHz'};


% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------

imagesc3d2(input.fat, size(input.fat)/2, 1, [90,90,90], [-0,1e-1]),  
imagesc3d2(input.water, size(input.fat)/2, 2, [90,90,90], [-0,1e-1]),  
imagesc3d2(input.r2s, size(input.fat)/2, 3, [90,90,90], [-0,1e2]),  
imagesc3d2(input.fieldHz, size(input.fat)/2, 4, [90,90,90], [-1e2,1e2]),  

for t = 1:6
    % imagesc3d2(input.img_cplx(:,:,:,t), size(input.fat)/2, 10+t, [90,90,90], [-0,5e-2]),  
    imagesc3d2(angle(input.img_cplx(:,:,:,t)), size(input.fat)/2, 10+t, [90,90,90], [-pi,pi]),  
end


% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------


for variable_index = 1:numel(required_variables)
    if ~isfield(input, required_variables{variable_index})
        error('Missing required variable %s in %s.', ...
            required_variables{variable_index}, input_mat)
    end
end

if size(input.img_cplx, 4) ~= numel(TE_seconds)
    error('Expected %d echoes, but img_imoco contains %d.', numel(TE_seconds), size(input.img_cplx, 4))
end

complex_images = permute(input.img_cplx, [3, 2, 1, 4]);

water = permute(input.water, [3, 2, 1]);
fat = permute(input.fat, [3, 2, 1]);

r2star_hz = permute(input.r2s, [3, 2, 1]);
field_hz = permute(input.fieldHz, [3, 2, 1]);

matrix_size = size(field_hz);
voxel_size_mm = FOV_mm ./ matrix_size;

B0_direction = [0, 0, 1];   % qsm_ref10_resolve_B0_direction(orientation_nifti_directory) gives exactly [0 0 1]


fat_reference_roi = [];

if ~isempty(fat_reference_roi_mat)
    roi_input = load(fat_reference_roi_mat);
    if ~isfield(roi_input, 'fat_reference_roi')
        error('%s must contain fat_reference_roi.', fat_reference_roi_mat)
    end
    fat_reference_roi = logical(roi_input.fat_reference_roi);
end



% -------------------------------------------------------------------------
%% reorient espirit mask, perform total segmentator for bounding box: liver
% -------------------------------------------------------------------------

use_espirit = 1;
use_imfill = 1;   % to fill in holes inside coil sensitivity mask

if use_espirit
    mask_esp_crop = receive(:,:,:,1)~=0;
    
    % it seems espirit mask needs to be cropped due to readout oversampling
    mask_esp_crop = mask_esp_crop(1+end/2-70:end/2+70,:,:);
    mask2use = permute(mask_esp_crop, [3,2,1]);
else
    % coil sensitivity mask not available, use MIP mask
    mask2use = mask_A;
end

imagesc3d2(water .* mask2use, matrix_size/2, 100, [90,90,90], [-0,1e-2]), setGcf(.5)
imagesc3d2(imfill(mask2use,'holes'), matrix_size/2, 101, [90,90,90], [-0,1e-2]), setGcf(.5)
imagesc3d2(mask2use, matrix_size/2, 102, [90,90,90], [-0,1e-2]), setGcf(.5)

if use_imfill
    % this seems to be needed to remove zeroed out voxels at the top of the
    % liver which otherwise propogate to qsm with lbv/smv erosion
    mask2use = imfill(mask2use,'holes');
end


% -------------------------------------------------------------------------
%% mask out slices without liver: union of TotalSegmentator (total_mr) livers
% -------------------------------------------------------------------------

% on the echo-2 magnitude, 6-echo RSS and water image, S-I (dim3) extent + 10 mm.
% liver_mask and liver_bbox were computed once with qsm_ref10_liver_bbox
% (TotalSegmentator total_mr) and are loaded here; liver_mask also sets the
% W_op normalization below
load(fullfile(script_directory, 'data', 'liver_mask.mat'), 'liver_mask', 'liver_bbox')


% -------------------------------------------------------------------------
%% bounding box mask
% -------------------------------------------------------------------------

use_bounding_box = 0;    % set to 1 to use bounding box in addition to coil sensitivity mask
select_only_slices = 0;  % set to 1 to use bounding box only along slice axis 
use_manual_mask = 0;     % set to 1 to remove dead voxels (not needed with bounding box in all axes)

mask_bbox = ones(size(liver_mask));

if use_bounding_box
    % select slices based on bounding box
    mask_bbox(:,:,1:liver_bbox(3,1)-1,:,:) = 0;
    mask_bbox(:,:,liver_bbox(3,2)+1:end,:,:) = 0;


    if ~select_only_slices
        % use bounding box in all 3 axes
        mask_bbox(1:liver_bbox(1,1)-1,:,:) = 0;
        mask_bbox(liver_bbox(1,2)+1:end,:,:) = 0;
        
        mask_bbox(:,1:liver_bbox(2,1)-1,:,:) = 0;
        mask_bbox(:,liver_bbox(2,2)+1:end,:,:) = 0;
    end
end

if use_manual_mask
    % zero out dead voxels manually
    mask_bbox(:,10,:) = 0;
    mask_bbox(:,6,:) = 0;
    mask_bbox(:,92,:) = 0;
end


imagesc3d2(water .* liver_mask, matrix_size/2, 100, [90,90,90], [-0,3e-2]), setGcf(.5)
imagesc3d2(water .* mask_bbox, matrix_size/2, 101, [90,90,90], [-0,3e-2]), setGcf(.5)

imagesc3d2(water .* mask2use, matrix_size/2, 1, [90,90,90], [-0,3e-2]), setGcf(.5)
imagesc3d2(fat .* mask2use, matrix_size/2, 2, [90,90,90], [-0,3e-2]), setGcf(.5)

imagesc3d2(mask2use, matrix_size/2, 3, [90,90,90], [-0,1]), setGcf(.5)
imagesc3d2(mask_bbox, matrix_size/2, 4, [90,90,90], [-0,1]), setGcf(.5)

imagesc3d2(r2star_hz .* mask2use, matrix_size/2, 4, [90,90,90], [-0,100]), setGcf(.5), colormap hot
imagesc3d2(field_hz .* mask2use, matrix_size/2, 5, [90,90,90], [-100,100]),


% -------------------------------------------------------------------------
%% combine bounding box mask (if any) with coil sens mask
% -------------------------------------------------------------------------

mask2use = mask2use .* mask_bbox;

water = water .* mask2use;
fat = fat .* mask2use;

field_hz = field_hz .* mask2use;
r2star_hz = r2star_hz .* mask2use;

imagesc3d2(field_hz .* mask2use, matrix_size/2, 6, [90,90,90], [-100,100]),


% -------------------------------------------------------------------------
%% LBV background removal (Zhou et al., NMR Biomed 2014) on the ESPIRiT mask
% -------------------------------------------------------------------------

% Used only by the hand-built and fat-referenced LSQR below; the PCG/LSQR
% calls to qsm_ref10_reconstruct still use field_hz.
% depth=0 (no coarse grids) iterates the full-resolution grid to tol, as in
% script_buda_qsm_revsharp_v1.m. On a synthetic test with this mask it gave
% 0.026 Hz local-field error (vs 0.11 Hz for MEDI defaults) in ~50 s.
% peel=0 keeps the subcutaneous fat used for fat referencing; padding is not
% needed although the mask touches the array faces (same result padded).

lbv_tol = 0.01;
lbv_depth = 0;
lbv_peel = 0;
lbv_N1 = 100;
lbv_N2 = 100;
lbv_N3 = 100;

% LBV is linear, so Hz in -> Hz out
tic
    field_lbv_hz = LBV(field_hz, mask2use, matrix_size, voxel_size_mm, lbv_tol, lbv_depth, lbv_peel, lbv_N1, lbv_N2, lbv_N3);
toc

% LBV zeroes the outer mask layer even with peel = 0
mask_lbv = field_lbv_hz ~= 0;


imagesc3d2(field_hz, matrix_size/2, 5, [90,90,90], [-20,20]*10),
imagesc3d2(field_lbv_hz, matrix_size/2, 6, [90,90,90], [-20,20]*10),
imagesc3d2(field_hz - field_lbv_hz, matrix_size/2, 7, [90,90,90], [-20,20]*10),


% -------------------------------------------------------------------------
%% build operator: W 
% -------------------------------------------------------------------------

% combined mask 2 use
mask2use_final = mask2use .* mask_lbv;

% weighting operator -> rsos normalized at 99th percentile
complex_images = complex_images .* mask2use_final;

W_op = sum(abs(complex_images).^2,4).^.5;
W_op = W_op / prctile(W_op(liver_mask),99);    % use liver mask so that W_op is agnostic to the selection of the bounding box 

imagesc3d2(W_op, matrix_size/2, 1, [90,90,90], [-0,1]),

[mean(W_op(:)), median(W_op(:))]
 

% -------------------------------------------------------------------------
% L operator
% -------------------------------------------------------------------------

% SMV kernel: 6 mm sphere on the factor-2 padded grid (Ref10), built the
% same way as make_kernels in qsm_ref10_reconstruct.m

smv_radius_mm = 6;
pad_size = round(2 * matrix_size);

% voxel-centre coordinates in mm, origin at the grid centre
x_mm = single((-floor(pad_size(1)/2):ceil(pad_size(1)/2)-1) * voxel_size_mm(1));
y_mm = single((-floor(pad_size(2)/2):ceil(pad_size(2)/2)-1) * voxel_size_mm(2));
z_mm = single((-floor(pad_size(3)/2):ceil(pad_size(3)/2)-1) * voxel_size_mm(3));

r2_mm = reshape(x_mm.^2, [], 1, 1) + reshape(y_mm.^2, 1, [], 1) + reshape(z_mm.^2, 1, 1, []);

sphere = single(r2_mm <= smv_radius_mm^2);
sphere = sphere / sum(sphere(:));

% k-space kernels: SMV, and L = I - SMV (zero at DC)
SMV_k = real(fftn(ifftshift(sphere)));
L_k = 1 - SMV_k;
L_k(1,1,1) = 0;

% apply L to the field: centre-pad, filter in k-space, crop back
idx = arrayfun(@(d) floor((pad_size(d) - matrix_size(d))/2) + (1:matrix_size(d)), 1:3, 'UniformOutput', false);

field_ppm = field_hz .* mask2use_final / (42.58 * B0_strength_t);

field_pad = zeros(pad_size, 'single');
field_pad(idx{:}) = field_ppm;

L_field_ppm = real(ifftn(L_k .* fftn(field_pad)));
L_field_ppm = L_field_ppm(idx{:});

imagesc3d2(L_field_ppm, matrix_size/2, 15, [90,90,90], [-1,1]/2),


% -------------------------------------------------------------------------
%% dipole kernel
% -------------------------------------------------------------------------

kernel = dipole_kernel(pad_size, voxel_size_mm, B0_direction);

test = real(ifftn(kernel .* fftn(field_pad)));
test = test(idx{:});

imagesc3d2(test, matrix_size/2, 16, [90,90,90], [-1,1]/2),


% -------------------------------------------------------------------------
%% create edge weight
% -------------------------------------------------------------------------

fat_fraction = fat ./ (eps + water + fat);
fat_fraction = fat_fraction .* mask2use_final;

fat_thre = 0.8;
fat_mask_thresholded = fat_fraction > 0.8;

imagesc3d2(fat_fraction, matrix_size/2, 16, [90,90,90], [-0,1]),
imagesc3d2(fat_mask_thresholded, matrix_size/2, 17, [90,90,90], [-0,1]),


% RSoS image edges -> top 20% (inside the foreground mask)
edge_R2star = qsm_ref10_gradient(W_op, 'notransp');

mask_R2star_edge = zeros([matrix_size,3],'logical');

for t = 1:3
    tmp = abs(edge_R2star(:,:,:,t));
    % edge_top20_thre = prctile(tmp(:),80);   % over the whole FOV: zeros outside the mask made this ~60% of tissue, 100% with the 3D box
    edge_top20_thre = prctile(tmp(mask2use_final > 0),80);   % top 20% of the masked tissue

    mask_R2star_edge(:,:,:,t) = abs(edge_R2star(:,:,:,t)) > edge_top20_thre;

    imagesc3d2(mask_R2star_edge(:,:,:,t), matrix_size/2, 20+t, [90,90,90], [-0,1]),
end


% chem sat edges -> top 20%
TE_index_chemShift = 1;     % 1 ms, similar to 1.2 ms in ref10 (1.16 ms at 2.89T ideally)

% paper seems to use magnitude of source images
edge_chemShift = qsm_ref10_gradient(abs(complex_images(:,:,:,TE_index_chemShift)) .* mask2use_final, 'notransp');


mask_chemShift_edge = zeros([matrix_size,3],'logical');

for t = 1:3
    tmp = abs(edge_chemShift(:,:,:,t));
    % edge_top20_thre = prctile(tmp(:),80);   % over the whole FOV: zeros outside the mask made this ~60% of tissue, 100% with the 3D box
    edge_top20_thre = prctile(tmp(mask2use_final > 0),80);   % top 20% of the masked tissue

    mask_chemShift_edge(:,:,:,t) = abs(edge_chemShift(:,:,:,t)) > edge_top20_thre;

    imagesc3d2(mask_chemShift_edge(:,:,:,t), matrix_size/2, 30+t, [90,90,90], [-0,1]),
end

% combined mask
mask_R2star_chemShift_combo = (mask_R2star_edge + mask_chemShift_edge) > 0;



% binary C: clip fat non-edge weight from 2 to 1, i.e. C = ~edge | fat80
% (fat edges keep weight 1, non-fat edges stay 0)
use_binary_C = false;

C_op = ~mask_R2star_chemShift_combo + fat_mask_thresholded;

if use_binary_C
    C_op = min(C_op, 1);
end

WC_op = W_op .* C_op;

for t = 1:3
    imagesc3d2(C_op(:,:,:,t), matrix_size/2, 10+t, [90,90,90], [-0,2]),
end


% -------------------------------------------------------------------------
%% create fat mask as in velikina et al., MRM 2023, Eq. 8
% -------------------------------------------------------------------------

fat_mask_v1 = (fat_fraction .* mask2use_final) > 0.9;

fat_mask_v2 = r2star_hz < 300;       % has not effect, R2star < 100 Hz
fat_mask_v3 = fat_mask_v1 & fat_mask_v2;

imagesc3d2(fat_mask_v1, matrix_size/2, 17, [90,90,90], [-0,1]),
imagesc3d2(fat_mask_v2, matrix_size/2, 18, [90,90,90], [-0,1]),
imagesc3d2(fat_mask_v3, matrix_size/2, 19, [90,90,90], [-0,1]),


% -------------------------------------------------------------------------
%% tfi preconditioner
% -------------------------------------------------------------------------

Pb_values = [1, 5, 30, 80];  % reported values in the paper
Pb_index2use = 4;

P_op = zeros(matrix_size);

mask_tfi_v1 = mask2use_final;

% mask_tfi_v2 = r2star_hz < 30;    % from liu et al tfi paper   -> excludes liver from the foreground
mask_tfi_v2 = r2star_hz < 100;       

mask_tfi_v3 = mask_tfi_v1 & mask_tfi_v2;       % inside this mask P=1 (foreground)

P_op(mask_tfi_v3) = 1;
P_op(~mask_tfi_v3) = Pb_values(Pb_index2use);  

imagesc3d2(mask_tfi_v1, matrix_size/2, 7, [90,90,90], [-0,1]),
imagesc3d2(mask_tfi_v2, matrix_size/2, 8, [90,90,90], [-0,1]),
imagesc3d2(mask_tfi_v3, matrix_size/2, 9, [90,90,90], [-0,1]),
imagesc3d2(W_op, matrix_size/2, 10, [90,90,90], [-0,1]),
imagesc3d2(P_op, matrix_size/2, 11, [90,90,90], [-0,max(P_op(:))]),


% -------------------------------------------------------------------------
%% perform LSQR with fat referecing, LBV pre-processing and TFI precond
% -------------------------------------------------------------------------

% LSQR recon
lsqr_iter = 1000;
lsqr_tol = 1e-3;

use_lbv = 1;     % set to 1 to use LBV pre-processing
use_laplace = 0; % set to 0 to remove laplace from forward model (only dipole kernel used)
use_precond = 1; % set to 1 to use TFI preconditioning

% apply L to the field: centre-pad, filter in k-space, crop back
idx = arrayfun(@(d) floor((pad_size(d) - matrix_size(d))/2) + (1:matrix_size(d)), 1:3, 'UniformOutput', false);

if use_lbv
    % field_ppm = field_lbv_hz / (42.58 * B0_strength_t);   % LBV tissue field
    field_ppm = field_lbv_hz .* mask2use_final / (42.58 * B0_strength_t);   % LBV tissue field
else
    field_ppm = field_hz / (42.58 * B0_strength_t);
end

field_pad = zeros(pad_size, 'single');
field_pad(idx{:}) = field_ppm;

if use_laplace
    L_k2use = L_k;
else
    L_k2use = 1;
end

L_field_ppm = real(ifftn(L_k2use .* fftn(field_pad)));
L_field_ppm = L_field_ppm(idx{:});

if use_lbv
    % rhs = W_op .* mask_lbv .* L_field_ppm;   % data weight masked: LBV field is 0 outside mask_lbv
    rhs = W_op .* mask2use_final .* L_field_ppm;   % data weight masked: LBV field is 0 outside mask_lbv
else
    rhs = W_op .* L_field_ppm;   % data weight masked: LBV field is 0 outside mask_lbv
end

rhs2use = cat(1, rhs(:), zeros([4*prod(matrix_size),1], 'single'));   % was 3*

param = [];
if use_lbv  
    % param.W_op = W_op .* mask_lbv;
    param.W_op = W_op .* mask2use_final;
else
    param.W_op = W_op;
end

param.WC_op = WC_op;

if use_precond
    param.P_op = P_op;
else
    param.P_op = 1;    
end

param.pad_size = pad_size;
param.matrix_size = matrix_size;
param.Lk_dipole_kernel = L_k2use .* kernel;
param.lambda = 1e-2;

param.fat_mask = fat_mask_v3;      % logical, matrix_size
param.lambda_fat = 1e-2;

param.print_every = 50;

% to provide initial guess:
chi_init = [];

if ~isempty(chi_init)
    x0 = chi_init ./ param.P_op;   % lsqr solves for y, chi = P_op .* y
    x0 = x0(:);
else
    x0 = [];
end

tic
    clear apply_ref10_fatref_tfi_counted   % reset lsqr iteration counter
    [res, flag, relres, iter, resvec, lsvec] = lsqr(@apply_ref10_fatref_tfi_counted, rhs2use, lsqr_tol, lsqr_iter, [], [], x0, param);
toc

% flag = 0 -> converged
disp(['flag: ', num2str(flag), ' relres: ', num2str(relres) , ' iter: ', num2str(iter)])


chi_fres_tfi = param.P_op .* reshape(res, matrix_size);
chi_fres_tfi = chi_fres_tfi - mean(chi_fres_tfi(param.fat_mask));   % fat reference


imagesc3d2(chi_fres_tfi, matrix_size/2, 3, [90,90,90], [-1,1]/1, [], 'w/ LBV pre-proc, w/ TFI precond, w/o Laplace'),
imagesc3d2(chi_fres_tfi .* mask2use_final, matrix_size/2, 23, [90,90,90], [-1,1]/1, [], 'w/ LBV pre-proc, w/ TFI precond, w/o Laplace'),


genNiigz(chi_fres_tfi .* mask2use_final, voxel_size_mm,    fullfile(output_directory, 'qsm_imoco_PT_lbv_FatReg_noLap_tfi_v1.nii'))


% -------------------------------------------------------------------------
%% MEDI with SMV and fat regularization (Kang et al., MRM 2023, Eq. 3)
% -------------------------------------------------------------------------

% MEDI+0 from the chi-separation MEDI toolbox (as in script_abdominal_xd_flash.m),
% with the fat mask in place of the CSF mask, which gives Kang's fat term.
% - LBV field and mask2use_final, as in the LSQR cells above
% - SMV radius 5 mm (MEDI default, Kang gives none)
% - N_std = 1/RSoS, so the data weight follows W_op (MEDI normalizes it to mean 1)
% - lambda = 1000, lambda_CSF = 200: Kang's lambda1 = 1e-3 = 1/lambda and
%   lambda2 = 0.1 = lambda_CSF/(2 lambda)
% - zero padding to the 2x grid used above; merit off, percentage 0.9 (MEDI defaults)
% MEDI erodes Mask by the SMV but not Mask_CSF, and zeroes chi outside the
% eroded mask before referencing to the Mask_CSF mean, so the fat mask is
% restricted to the eroded mask here.

medi_toolbox_directory = fullfile(script_directory, 'external', 'MEDI_toolbox', 'functions');
medi_directory = fullfile(pwd, 'medi_v2');      % MEDI reads its inputs from a file and writes results/ to the current folder

medi_radius_mm = 5;
% medi_lambda = 1000;
medi_lambda = 500;
% medi_lambda_fat = 200;
medi_lambda_fat = 100;

delta_TE = TE_seconds(2) - TE_seconds(1);       % MEDI works in radians at the echo spacing
CF = 42.58e6 * B0_strength_t;                    % Hz, same gyromagnetic ratio as the ppm conversion above

% it seems MEDI does not use iFreq, but RDF instead. to match Kang et al,
% LBV'ed field map is provided as input
RDF = field_lbv_hz .* mask2use_final * 2*pi * delta_TE;

assert(all(W_op(mask2use_final > 0) > 0), 'RSoS must be positive inside the mask')

N_std = ones(matrix_size);                       % finite outside the mask: MEDI multiplies N_std by Mask
N_std(mask2use_final > 0) = 1 ./ W_op(mask2use_final > 0);

disp(['|RDF| > pi (phase wraps in the MEDI data term) in ', num2str(nnz(abs(RDF) > pi)), ' of ', num2str(nnz(mask2use_final)), ' mask voxels'])

if ~isfolder(medi_directory)
    mkdir(medi_directory)
end

medi_input_file = fullfile(medi_directory, 'RDF_medi_v2.mat');

% the chi-separation toolbox shares function names (SMV, dipole_kernel,
% MEDI_L1, ...) with LIBRARY/MEDI_toolbox, so it is on the path only here
saved_path = path;
addpath(genpath(medi_toolbox_directory))

current_directory = pwd;

try
    % SMV-eroded mask exactly as inside MEDI_L1 (padded grid, > 0.999)
    sphere_k_medi = single(sphere_kernel(2*matrix_size, voxel_size_mm, medi_radius_mm));
    mask_smv_medi = SMV(padarray(double(mask2use_final), matrix_size, 'post'), sphere_k_medi) > 0.999;
    mask_smv_medi = mask_smv_medi(1:matrix_size(1), 1:matrix_size(2), 1:matrix_size(3));

    fat_mask_medi = fat_mask_v3 & mask_smv_medi;

    medi_input = struct('iFreq', field_hz .* mask2use_final * 2*pi * delta_TE, 'RDF', RDF, 'N_std', N_std, 'iMag', W_op, ...
        'Mask', double(mask2use_final), 'matrix_size', matrix_size, 'voxel_size', voxel_size_mm, 'delta_TE', delta_TE, 'CF', CF, 'B0_dir', B0_direction, 'Mask_CSF', fat_mask_medi);

    save(medi_input_file, '-struct', 'medi_input')

    cd(medi_directory)

    tic
        chi_medi = MEDI_L1('filename', medi_input_file, 'lambda', medi_lambda, 'smv', medi_radius_mm, 'lambda_CSF', medi_lambda_fat, 'zeropad', matrix_size);
    toc

catch medi_error
    cd(current_directory)
    path(saved_path)
    rethrow(medi_error)
end

cd(current_directory)
path(saved_path)

% ppm, zero outside mask_smv_medi, referenced to the mean over fat_mask_medi
chi_medi = double(chi_medi);

disp(['fat voxels: ', num2str(nnz(fat_mask_medi)), ' of ', num2str(nnz(fat_mask_v3)), ' in fat_mask_v3 are inside the SMV-eroded mask'])
disp(['eroded mask: ', num2str(nnz(mask_smv_medi)), ' of ', num2str(nnz(mask2use_final)), ' voxels; NaN/Inf in chi_medi: ', num2str(nnz(~isfinite(chi_medi)))])


imagesc3d2(chi_medi, matrix_size/2, 61, [90,90,90], [-1.2,.4]/1, [], 'MEDI + SMV + fat term (Kang)'),


genNiigz(chi_medi .* mask2use_final, voxel_size_mm, fullfile(output_directory, ['qsm_imoco_PT_medi_smv', num2str(medi_radius_mm), 'mm_lam', num2str(medi_lambda), '_fat', num2str(medi_lambda_fat), '_v1.nii']))


% -------------------------------------------------------------------------
%% compare MEDI with the LSQR solves in the liver
% -------------------------------------------------------------------------

% every map re-referenced to the same fat mask (fat_mask_medi) and compared
% in liver_mask & mask_smv_medi, where all of them have data

mask_compare = liver_mask & mask_smv_medi;

chi_compare = {chi_fres_tfi, chi_medi};
chi_names = {'TFI w/ fat ref, w/o L', 'MEDI + SMV + fat (Kang)'};

for k = 1:numel(chi_compare)
    chi_k = chi_compare{k} - mean(chi_compare{k}(fat_mask_medi));
    q = prctile(chi_k(mask_compare), [25, 50, 75]);
    fprintf('%-24s liver median %6.3f ppm (IQR %6.3f to %6.3f)\n', chi_names{k}, q(2), q(1), q(3))
end

chi_diff = (chi_medi - (chi_fres_tfi - mean(chi_fres_tfi(fat_mask_medi)))) .* mask_smv_medi;

imagesc3d2(chi_diff, matrix_size/2, 61, [90,90,90], [-1,1]/2, [], 'MEDI minus TFI w/ fat ref'),
