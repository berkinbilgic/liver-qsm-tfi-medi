function [ output_args ] = genNiigz( img, voxel_size, filename )
%GENNII Summary of this function goes here
%   Detailed explanation goes here

if norm2(imag(img))~=0
    disp('imaginary component detected, taking abs')
    img = abs(img);
end

save_nii(make_nii(single(img), voxel_size, [], 16), filename)

system(['gzip ', filename])



end

