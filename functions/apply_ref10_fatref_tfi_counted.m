function [ res, tflag ] = apply_ref10_fatref_tfi_counted( in, params, tflag )
%APPLY_REF10_FATREF_COUNTED apply_ref10_fatref plus an LSQR iteration counter.
%   With x0 = [], lsqr makes one 'notransp' call per iteration, so counting
%   forward calls gives the iteration number. Prints every
%   params.print_every iterations (default 10).
%   Reset the counter before each solve with:  clear apply_ref10_fatref_counted

persistent n_iter t_start

if isempty(n_iter)
    n_iter = 0;
    t_start = tic;
end

if ~strcmp(tflag,'transp')
    n_iter = n_iter + 1;

    print_every = 10;
    if isfield(params, 'print_every')
        print_every = params.print_every;
    end

    if mod(n_iter, print_every) == 0
        fprintf('lsqr iter %d  (%.1f s)\n', n_iter, toc(t_start));
    end
end

[res, tflag] = apply_ref10_fatref_tfi(in, params, tflag);

end
