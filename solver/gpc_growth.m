function source = gpc_growth(rho, K, cfg)
%GPC_GROWTH Project h*rho with h=0.5*(1-0.1*z)*c.
% The radial nutrient field uses the support of the mean density.
if nargin < 3, cfg = experiment_config(); end
persistent saved_count saved_quad matrix
if isempty(matrix) || saved_count ~= K || saved_quad ~= cfg.n_growth_quadrature
    [z,w] = lgwt(cfg.n_growth_quadrature,-1,1);
    psi = legendre_basis(z,K);
    matrix = psi' * (psi .* (0.5*w(:) .* (0.5*(1-0.1*z(:)))));
    saved_count = K; saved_quad = cfg.n_growth_quadrature;
end
C = nutrient_field(rho(:,:,1),cfg);
projected = reshape(reshape(rho,[],K)*matrix',size(rho));
source = projected .* C;
end
