% 基础函数：保留数值更新公式；由 run_experiment 统一调用。
% Original author solver; executable formulas retained.
% Xu'an Dou attribution is retained where present in the original source.
function rhoGGG = gpc_growth(rho,K)
    % Extract the first 45x45 slice of R
cB = 1;
dx = 0.1;
% dy = 0.1;
a = -2.2; b = 2.2;
[X, Y] = meshgrid(a:dx:b, a:dx:b);
Nx = size(X, 1);
Ny = size(X, 2);
    R_slice = rho(:, :, 1);
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
Nz = 32;  
modelType = 2;
[z,wz] = lgwt(Nz,-1,1);  % 1 x Nz
 wz = 1/2*wz;
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
    C(r <= R) = term1 * besseli(0, r(r <= R)); % r <= R
    term2 = cB * I1_R / (K0_R * I1_R + K1_R * I0_R);
    C(r > R) = cB - term2 * besselk(0, r(r > R)); % r > R
end
% GFun = @(C, X, Y, z) C .* ((1.5 - 0.3 * z) .* ((X.^2 + Y.^2) <= R^2) + ...
%                            (1 - 0.2 * z) .* ((X.^2 + Y.^2) > R^2));
% GFun = @(C, X, Y, z) C .* ((0.5 - 0.001 * z) .* ((X.^2 + Y.^2) <= R^2) + ...
%                            (1 - 0. * z) .* ((X.^2 + Y.^2) > R^2));
GFun = @(C, X, Y, z) C .* ((0.5 - 0.05 * z) .* ((X.^2 + Y.^2) <= R^2) + ...
                           (1 - 0. * z) .* ((X.^2 + Y.^2) > R^2));                       
psi = zeros(Nz, K);
for k = 1:7
    if k == 1
        psi(:, k) = sqrt(1) * 1;  % Psi_0(z) = sqrt(1)
    elseif k == 2
        psi(:, k) = sqrt(3) .* z;  % Psi_1(z) = sqrt(3) * z
    elseif k == 3
        psi(:, k) = sqrt(5) .* (3/2 * z.^2 - 1/2);  % Psi_2(z) = sqrt(5) * (3/2 * z^2 - 1/2)
    elseif k == 4
        psi(:, k) = sqrt(7) .* (5/2 * z.^3 - 3/2 * z);  % Psi_3(z) = sqrt(7) * (5/2 * z^3 - 3/2 * z)
    elseif k == 5
        psi(:, k) = sqrt(9) .* (35/8 * z.^4 - 15/4 * z.^2 + 3/8);  % Psi_4(z)
    elseif k == 6
        psi(:, k) = sqrt(11) .* (63/8 * z.^5 - 35/4 * z.^3 + 15/8 * z);  % Psi_5(z)
    elseif k == 7
        psi(:, k) = sqrt(13) .* (231/16 * z.^6 - 315/16 * z.^4 + 105/16 * z.^2 - 5/16);  % Psi_6(z)
    end
end
G = zeros(size(X, 1), size(X, 2), Nz);
for k = 1:Nz
    G(:,:,k) = GFun(C,X, Y, z(k));
end
rhoGG = zeros(size(X, 1), size(X, 2), K, K);
for i = 1:K
    for j = 1:K
        for k = 1:Nz
            rhoGG(:,:,i,j) = rhoGG(:,:,i,j) + G(:,:,k) .* psi(k,i) .* psi(k,j) .* wz(k);
        end
    end
end
% figure;
% subplot(1, 2, 1);
% imagesc(rhoGG(:,:,1,1)); % psi_1 * psi_1
% colorbar;
% title('Result of GFun .* psi_1 .* psi_1');
% 
% subplot(1, 2, 2);
% imagesc(rhoGG(:,:,1,2)); % psi_2 * psi_2
% colorbar;
% title('Result of GFun .* psi_2 .* psi_2');
rhoGGG = zeros(size(X, 1), size(X, 2), K);
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
         rhoVec = reshape(rho(i, j, :), [1, K]);
%         disp(size(rhoVec));
         rhoMat = reshape(rhoGG(i, j, :, :), [K, K]);
%                 rhoMat = rhoGG(i, j, :, :);
        rhoGGG(i, j, :) = rhoVec * rhoMat;
    end
end
% figure;
% for k = 1:4
%     title(['Layer ' num2str(k)]);
%     xlabel('X');
%     ylabel('Y');
% end
end
