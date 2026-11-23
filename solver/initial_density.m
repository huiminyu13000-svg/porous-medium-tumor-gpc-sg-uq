function rho = initial_density(X, Y, z)
%INITIAL_DENSITY Concave radial profile with continuous linear transitions.
% The four thresholds are geometric radii, so r = sqrt(X.^2 + Y.^2).
radius = sqrt(X.^2 + Y.^2);
edges = [0.4 0.5 0.6 0.7] * (1 + 0.3*z);
core = 0.15 - 0.1*z; rim = 0.35 + 0.1*z;
rho = zeros(size(X)); rho(radius <= edges(1)) = core;
mask = radius > edges(1) & radius <= edges(2);
rho(mask) = core + (rim-core) * (radius(mask)-edges(1)) / (edges(2)-edges(1));
rho(radius > edges(2) & radius <= edges(3)) = rim;
mask = radius > edges(3) & radius <= edges(4);
rho(mask) = rim * (edges(4)-radius(mask)) / (edges(4)-edges(3));
end
