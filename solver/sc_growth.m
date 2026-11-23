function growth = sc_growth(rho, k, z_nodes, cfg)
%SC_GROWTH Deterministic growth coefficient at one collocation parameter.
if nargin < 4, cfg = experiment_config(); end
growth = 0.5*(1-0.1*z_nodes(k))*nutrient_field(rho,cfg);
end
