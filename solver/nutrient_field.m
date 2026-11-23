function C = nutrient_field(rho, cfg)
%NUTRIENT_FIELD Radial in-vivo nutrient solution for the current support.
grid = cfg.domain(1):cfg.dx:cfg.domain(2);
[X,Y] = meshgrid(grid); radius = sqrt(X.^2+Y.^2);
occupied = rho > 0;
if ~any(occupied(:)), C = cfg.cB*ones(size(rho)); return; end
R = max(radius(occupied)); if R <= eps, R = cfg.dx/2; end
denominator = besselk(0,R)*besseli(1,R) + besselk(1,R)*besseli(0,R);
C = zeros(size(rho)); inside = radius <= R;
C(inside) = cfg.cB*besselk(1,R)/denominator*besseli(0,radius(inside));
C(~inside) = cfg.cB-cfg.cB*besseli(1,R)/denominator*besselk(0,radius(~inside));
end
