function [E, k, E_all] = MeanNewton(M, ecc, opts)
% Requires the symbolic toolbox

    % sets maximum number of iterations (defaults to 200)
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'k_max')
        k_max = 200;
    else
        k_max = opts.k_max;
    end
    
    % determines if all intermediate estimates should be returned
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'return_all')
        return_all = false;
    else
        return_all = opts.return_all;
    end
    
    % sets tolerance (defaults to 10⁻¹⁰)
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'TOL')
        TOL = 1e-10;
    else
        TOL = opts.TOL;
    end


    syms f(E) fdE(E)
    f(E) = M-(E-ecc*sin(E));
    fdE(E) = ecc*cos(E)-1;
end
