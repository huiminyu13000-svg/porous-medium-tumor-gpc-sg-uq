% 基础函数：保留数值更新公式；由 run_experiment 统一调用。
% Original author solver; executable formulas retained.
% Xu'an Dou attribution is retained where present in the original source.
function [rhoGGG,rhoAAA,rhoBBB,rhoDDD] = sc_growth(rho,k,z_nodes)
    % Extract the first 45x45 slice of R
cB = 1;
dx = 0.1;
% dy = 0.1;
a = -2.2; b = 2.2;
[X, Y] = meshgrid(a:dx:b, a:dx:b);
Nx = size(X, 1);
Ny = size(X, 2);
z=z_nodes(:); % Use the same nodes as the saved initial samples.
%  wz = 1/2*wz;
    R_slice = rho;
    % Define the center of the slice in matrix indices (23,23 for MATLAB indexing)
    center_idx = 1/2*(Nx+1);
    % Create coordinate grids for x and y
    [x_grid, y_grid] = meshgrid(((1:Nx) - center_idx) * dx, ((1:Nx) - center_idx) * dx);
    % Find the indices of all points where R_slice is 1
    [row_indices, col_indices] = find(R_slice >0);
    % Calculate the Euclidean distance from the center (0,0) to each point with R_slice == 1
    distances = sqrt((x_grid(sub2ind(size(R_slice), row_indices, col_indices))).^2 + ...
                     (y_grid(sub2ind(size(R_slice), row_indices, col_indices))).^2);
    % The radius of the circle in the slice is the maximum distance found
    R = max(distances);
modelType = 2;
C = zeros(Nx, Ny);
[x, y] = meshgrid(a:dx:b, a:dx:b);
r = sqrt(x.^2 + y.^2);
I0_R = besseli(0, R);  % I0(R)
I1_R = besseli(1, R);  % I1(R)
K0_R = besselk(0, R);  % K0(R)
K1_R = besselk(1, R);  % K1(R)
if modelType == 1
    C(r <= R) = cB * I0_R * besseli(0, r(r <= R)); % r <= R
    C(r > R) = cB; % r > R
elseif modelType == 2
    term1 = cB * K1_R / (K0_R * I1_R + K1_R * I0_R);
%     
%     disp(size(K1_R));
%     disp(size(K0_R));
%        disp(size(I1_R));
%           disp(size(I0_R));
    C(r <= R) = term1 * besseli(0, r(r <= R)); % r <= R
    term2 = cB * I1_R / (K0_R * I1_R + K1_R * I0_R);
    C(r > R) = cB - term2 * besselk(0, r(r > R)); % r > R
end
% GFun = @(C, X, Y) C .* ((0.5-  0.001 * z(k)) .* ((X.^2 + Y.^2) <= R^2) + ...
%                            (1 - 0. * z(k)) .* ((X.^2 + Y.^2) > R^2));
GFun = @(C, X, Y) C .* ((0.5 -  0.05 * z(k)) .* ((X.^2 + Y.^2) <= R^2) + ...
                           (1 - 0. * z(k)) .* ((X.^2 + Y.^2) > R^2));
GAFun = @(X, Y)  (0.3-0.3*z(k));
GBFun = @(X, Y)  (0.2-0.2*z(k));
GDFun = @(X, Y)   ((0.5+0.5*z(k)) .* ((X.^2 + Y.^2) <= R^2) + ...
                           (1 - 0. * z(k)) .* ((X.^2 + Y.^2) > R^2)); 
rhoGGG = zeros(size(X, 1), size(X, 2));
rhoAAA = zeros(size(X, 1), size(X, 2));
rhoBBB = zeros(size(X, 1), size(X, 2));
rhoDDD = zeros(size(X, 1), size(X, 2));
    rhoGGG(:,:) = GFun(C, X, Y);
    rhoAAA(:,:) = GAFun(X, Y);
    rhoBBB(:,:) = GBFun(X, Y);
    rhoDDD(:,:) = GDFun(X, Y);
end
