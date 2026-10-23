function output = qsm_ref10_gradient(input, transpose_flag)
%QSM_REF10_GRADIENT Ref10 3D finite-difference operator G or its transpose.
%
% GRADIENT = QSM_REF10_GRADIENT(IMAGE) returns the [size(IMAGE), 3] array
% of forward differences along each dimension, with a zero high-side
% boundary value. IMAGE = QSM_REF10_GRADIENT(GRADIENT, 'transp') applies
% the exact transpose (a negative discrete divergence).
%
% Differences are per voxel, not per mm. Both local functions are copied
% unchanged from qsm_ref10_reconstruct.m, so lambda values used with this
% operator remain comparable with that solver.

if nargin < 2 || isempty(transpose_flag)
    transpose_flag = 'notransp';
end

if strcmp(transpose_flag, 'notransp')
    output = forward_gradient(input);
elseif strcmp(transpose_flag, 'transp')
    output = gradient_adjoint(input);
else
    error('qsm_ref10:GradientTransposeFlag', ...
        'transpose_flag must be ''notransp'' or ''transp''.');
end
end

function gradient = forward_gradient(image)
% Three forward finite differences with a zero high-side boundary value.
gradient = zeros([size(image), 3], 'like', image);
gradient(1:end-1,:,:,1) = image(2:end,:,:) - image(1:end-1,:,:);
gradient(:,1:end-1,:,2) = image(:,2:end,:) - image(:,1:end-1,:);
gradient(:,:,1:end-1,3) = image(:,:,2:end) - image(:,:,1:end-1);
end

function image = gradient_adjoint(gradient)
% Exact transpose of forward_gradient (a negative discrete divergence).
matrix_size = size(gradient);
matrix_size = matrix_size(1:3);
image = zeros(matrix_size, 'like', gradient);

gx = gradient(:,:,:,1);
image(1,:,:) = image(1,:,:) - gx(1,:,:);
if matrix_size(1) > 2
    image(2:end-1,:,:) = image(2:end-1,:,:) + gx(1:end-2,:,:) - gx(2:end-1,:,:);
end
image(end,:,:) = image(end,:,:) + gx(end-1,:,:);

gy = gradient(:,:,:,2);
image(:,1,:) = image(:,1,:) - gy(:,1,:);
if matrix_size(2) > 2
    image(:,2:end-1,:) = image(:,2:end-1,:) + gy(:,1:end-2,:) - gy(:,2:end-1,:);
end
image(:,end,:) = image(:,end,:) + gy(:,end-1,:);

gz = gradient(:,:,:,3);
image(:,:,1) = image(:,:,1) - gz(:,:,1);
if matrix_size(3) > 2
    image(:,:,2:end-1) = image(:,:,2:end-1) + gz(:,:,1:end-2) - gz(:,:,2:end-1);
end
image(:,:,end) = image(:,:,end) + gz(:,:,end-1);
end
