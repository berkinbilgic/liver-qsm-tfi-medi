function [ res, tflag ] = apply_ref10_fatref_tfi( in, params, tflag )
%APPLY_REF10_FATREF apply_ref10 plus a fat-uniformity term, so lsqr solves
%
%   min_chi ||W (L D chi - L psi)||^2 + lambda ||W C G chi||^2
%           + lambda_fat ||M_fat (chi - mean_fat(chi))||^2
%
%   The last term is Velikina et al., MRM 2023, Eq. 8, with chi_0 taken as
%   the mean of chi itself inside M_fat (as in MEDI+0). It penalizes
%   variation within fat, not the fat level, so constants stay in the null
%   space: reference chi to mean(chi(fat_mask)) after the solve.
%   Its operator P = M - m*m'/n is symmetric, so forward = transpose.
%
%   Extra params: fat_mask (logical, matrix_size), lambda_fat (scalar).
%   rhs needs prod(matrix_size) more zeros: [W L psi; zeros(4*N,1)].

if strcmp(tflag,'transp')    
    % Transposed operator:

    chi = in(1:prod(params.matrix_size));    
    Wchi = params.W_op .* reshape(chi, params.matrix_size);

    reg = in(1+prod(params.matrix_size):4*prod(params.matrix_size));    
    WCreg = params.WC_op .* reshape(reg, [params.matrix_size,3]);

    fat = reshape(in(1+4*prod(params.matrix_size):end), params.matrix_size);
   

    idx = arrayfun(@(d) floor((params.pad_size(d) - params.matrix_size(d))/2) + (1:params.matrix_size(d)), 1:3, 'UniformOutput', false);
    
    Wchi_pad = zeros(params.pad_size, 'single');
    Wchi_pad(idx{:}) = Wchi;
    
    tmp = real(ifftn(params.Lk_dipole_kernel .* fftn(Wchi_pad)));
    res_dc = tmp(idx{:});
    res_dc = params.P_op .* res_dc;

    res_reg = qsm_ref10_gradient(WCreg, 'transp');
    res_reg = params.P_op .* res_reg;

    res_fat = params.fat_mask .* (fat - mean(fat(params.fat_mask)));
    res_fat = params.P_op .* res_fat;

    res = res_dc(:) + sqrt(params.lambda) * res_reg(:) + sqrt(params.lambda_fat) * res_fat(:);    
else
    
    % Forward operator:    
    chi = params.P_op .* reshape(in, params.matrix_size);

    % apply L*D to the field: centre-pad, filter in k-space, crop back
    idx = arrayfun(@(d) floor((params.pad_size(d) - params.matrix_size(d))/2) + (1:params.matrix_size(d)), 1:3, 'UniformOutput', false);
    
    chi_pad = zeros(params.pad_size, 'single');
    chi_pad(idx{:}) = chi;
    
    tmp = real(ifftn(params.Lk_dipole_kernel .* fftn(chi_pad)));
    res_dc = tmp(idx{:});
    res_dc = params.W_op .* res_dc;

    res_reg = params.WC_op .* qsm_ref10_gradient(chi, 'notransp');

    res_fat = params.fat_mask .* (chi - mean(chi(params.fat_mask)));
    
    res = cat(1, res_dc(:), sqrt(params.lambda)*res_reg(:), sqrt(params.lambda_fat)*res_fat(:));    
end

end
