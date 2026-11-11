% 基础函数：保留数值更新公式；由 run_experiment 统一调用。
% Original author solver; executable formulas retained.
% Xu'an Dou attribution is retained where present in the original source.
% The main solver. Please see the document attached for usage and examples.
% Xu'an Dou dxa@pku.edu.cn 2020.9.6
function [rho,p]=solve_gpc_sg(T,m,rho_initial)
%% Parsing the Input
% T = 0.10;
% m = 80;
 K = 4;
dx = 0.1;
dy = 0.1;
dt = 0.001; 
a = -2.2; b = 2.2;
Ntime = floor(T/dt);
Nx = length(a:dx:b) - 1;
Ny = length(a:dx:b) - 1;
% [X, Y] = meshgrid(a:dx:b, a:dx:b);
data.rho_gPC = rho_initial;
 rho_gPC = data.rho_gPC;  % Nx x Ny x Ns
% data = load('h_gPC.mat');
%  h_gPC = data.h_gPC;      % Nx x Ny x Ns
 Ns = size(rho_gPC, 3);
rho = rho_gPC;
Velocity_x = @(X, Y) (X * 0);
Velocity_y = @(X, Y) (Y * 0);
[X, Y] = meshgrid(a:dx:b, a:dx:b);
u_e = Velocity_x(X, Y);
u_e = u_e'; u_e = u_e(2:end-1, 2:end-1);
v_e = Velocity_y(X, Y);
v_e = v_e'; v_e = v_e(2:end-1, 2:end-1);
u_e = repmat(u_e, [1, 1, Ns]);
v_e = repmat(v_e, [1, 1, Ns]);
p = rho.^(m-1) * m / (m-1);
u = -(p(3:end,2:end-1,:) - p(1:end-2,2:end-1,:)) ./ (2*dx);
v = -(p(2:end-1,3:end,:) - p(2:end-1,1:end-2,:)) ./ (2*dy);
U_m = zeros(Nx, Ny, Ns);    % Nx x Ny x Ns
%     disp(size(U_m));
V_m = U_m;                  % Nx x Ny x Ns
uv_pp = U_m; uv_pm = U_m; uv_mm = U_m; uv_mp = U_m;   % Nx x Ny x Ns
vu_pp = U_m; vu_pm = U_m; vu_mm = U_m; vu_mp = U_m;   % Nx x Ny x Ns
V_up = U_m; V_low = U_m;    % Nx x Ny x Ns
U_up = U_m; U_low = U_m;    % Nx x Ny x Ns
b_sol = zeros(Nx * Ny * 2 * Ns, 1);
for i = 1:Ntime
    Build_Matrix;
    [sol,~] = gmres(@(x) fA(x), b_sol(:), 20);
    u_star = reshape(sol(1:(Nx-1)*(Ny-1)*Ns), Nx-1, Ny-1, Ns);
    v_star = reshape(sol((Nx-1)*(Ny-1)*Ns+1:end), Nx-1, Ny-1, Ns);
    Convec_Cor;
end
% rho = rho';
% p = p';
%  
function Build_Matrix
    r_ip = ((rho(3:end, 2:end-1, :) + rho(2:end-1, 2:end-1, :)) ./ 2).^(m-2);  % Nx-2 x Ny-1 x Ns
    r_im = ((rho(1:end-2, 2:end-1, :) + rho(2:end-1, 2:end-1, :)) ./ 2).^(m-2);  % Nx-2 x Ny-1 x Ns
    r_jp = ((rho(2:end-1, 3:end, :) + rho(2:end-1, 2:end-1, :)) ./ 2).^(m-2);   % Nx-1 x Ny-2 x Ns
    r_jm = ((rho(2:end-1, 1:end-2, :) + rho(2:end-1, 2:end-1, :)) ./ 2).^(m-2);  % Nx-1 x Ny-2 x Ns
    U_m = 1 + m * (dt/dx/dx) .* rho(2:end-1, 2:end-1, :) .* (r_ip + r_im);   % Nx-1 x Ny-1 x Ns
    V_m = 1 + m * (dt/dy/dy) .* rho(2:end-1, 2:end-1, :) .* (r_jp + r_jm);   % Nx-1 x Ny-1 x Ns
    U_up = -m * (dt/dx/dx) .* r_ip .* rho(3:end, 2:end-1, :);  % Nx-2 x Ny-1 x Ns
    U_up = U_up(1:end-1, :, :);
    U_low = -m * (dt/dx/dx) .* r_im .* rho(1:end-2, 2:end-1, :);  % Nx-2 x Ny-1 x Ns
    U_low = U_low(2:end, :, :);
    V_up = -m * (dt/dy/dy) .* r_jp .* rho(2:end-1, 3:end, :);  % Nx-1 x Ny-2 x Ns
    V_up = V_up(:, 1:end-1, :);
    V_low = -m * (dt/dy/dy) .* r_jm .* rho(2:end-1, 1:end-2, :);  % Nx-1 x Ny-2 x Ns
    V_low = V_low(:, 2:end, :);
    lam = m * dt / (dx * dy) / 4;
    uv_pp = -lam * (rho(3:end, 2:end-1, :).^(m-2)) .* rho(3:end, 3:end, :);  % Nx-2 x Ny-2 x Ns
    uv_pp = uv_pp(1:end-1, 1:end-1, :);  % i+1, j+1
    uv_pm = lam * (rho(3:end, 2:end-1, :).^(m-2)) .* rho(3:end, 1:end-2, :);  % Nx-2 x Ny-2 x Ns
    uv_pm = uv_pm(1:end-1, 2:end, :);  % i+1, j-1
    uv_mp = lam * (rho(1:end-2, 2:end-1, :).^(m-2)) .* rho(1:end-2, 3:end, :);  % Nx-2 x Ny-2 x Ns
    uv_mp = uv_mp(2:end, 1:end-1, :);  % i-1, j+1
    uv_mm = -lam * (rho(1:end-2, 2:end-1, :).^(m-2)) .* rho(1:end-2, 1:end-2, :);  % Nx-2 x Ny-2 x Ns
    uv_mm = uv_mm(2:end, 2:end, :);  % i-1, j-1
    vu_pp = -lam * (rho(2:end-1, 3:end, :).^(m-2)) .* rho(3:end, 3:end, :);  % Nx-1 x Ny-1 x Ns
    vu_pm = lam * (rho(2:end-1, 3:end, :).^(m-2)) .* rho(3:end, 1:end-2, :);  % Nx-1 x Ny-1 x Ns
    vu_pp = vu_pp(1:end-1, 1:end-1, :);  % i+1, j+1
    vu_pm = vu_pm(1:end-1, 2:end, :);  % i+1, j-1
    vu_mp = lam * (rho(2:end-1, 1:end-2, :).^(m-2)) .* rho(1:end-2, 3:end, :);  % Nx-1 x Ny-1 x Ns
    vu_mm = -lam * (rho(2:end-1, 1:end-2, :).^(m-2)) .* rho(1:end-2, 1:end-2, :);  % Nx-1 x Ny-1 x Ns
    vu_mp = vu_mp(2:end, 1:end-1, :);  % i-1, j+1
    vu_mm = vu_mm(2:end, 2:end, :);  % i-1, j-1
    rhoGGG = gpc_growth(rho,K);
    rhoG = rho.^(m-2) .* rhoGGG;  % Nx x Ny x Ns
    bu = u - m * (dt / dx / 2) .* (rhoG(3:end, 2:end-1, :) - rhoG(1:end-2, 2:end-1, :));  % Nx-2 x Ny-1 x Ns
    bv = v - m * (dt / dy / 2) .* (rhoG(2:end-1, 3:end, :) - rhoG(2:end-1, 1:end-2, :));  % Nx-1 x Ny-2 x Ns
    bu = bu - (U_m - 1) .* u_e;
    bu(1:end-1, :, :) = bu(1:end-1, :, :) - U_up .* u_e(2:end, :, :);
    bu(2:end, :, :) = bu(2:end, :, :) - U_low .* u_e(1:end-1, :, :);
    bv = bv - (V_m - 1) .* v_e;
    bv(:, 1:end-1, :) = bv(:, 1:end-1, :) - V_up .* v_e(:, 2:end, :);
    bv(:, 2:end, :) = bv(:, 2:end, :) - V_low .* v_e(:, 1:end-1, :);
    bu = reshape(bu, (Nx-1)*(Ny-1)*Ns, 1);
    bv = reshape(bv, (Nx-1)*(Ny-1)*Ns, 1);
    b_sol = [bu; bv];
end
function ax = fA(x_tmp)
    u_tmp = reshape(x_tmp(1:(Nx-1)*(Ny-1)*Ns), Nx-1, Ny-1, Ns);
    v_tmp = reshape(x_tmp((Nx-1)*(Ny-1)*Ns+1:end), Nx-1, Ny-1, Ns);
    ax_u = U_m .* u_tmp;  % Nx-1 x Ny-1 x Ns
    ax_u(1:end-1, :, :) = ax_u(1:end-1, :, :) + U_up .* u_tmp(2:end, :, :);
    ax_u(2:end, :, :) = ax_u(2:end, :, :) + U_low .* u_tmp(1:end-1, :, :);
    ax_u(1:end-1, 1:end-1, :) = ax_u(1:end-1, 1:end-1, :) + uv_pp .* v_tmp(2:end, 2:end, :);
    ax_u(1:end-1, 2:end, :) = ax_u(1:end-1, 2:end, :) + uv_pm .* v_tmp(2:end, 1:end-1, :);
    ax_u(2:end, 1:end-1, :) = ax_u(2:end, 1:end-1, :) + uv_mp .* v_tmp(1:end-1, 2:end, :);
    ax_u(2:end, 2:end, :) = ax_u(2:end, 2:end, :) + uv_mm .* v_tmp(1:end-1, 1:end-1, :);
    ax_v = V_m .* v_tmp;  % Nx-1 x Ny-1 x Ns
    ax_v(:, 1:end-1, :) = ax_v(:, 1:end-1, :) + V_up .* v_tmp(:, 2:end, :);
    ax_v(:, 2:end, :) = ax_v(:, 2:end, :) + V_low .* v_tmp(:, 1:end-1, :);
    ax_v(1:end-1, 1:end-1, :) = ax_v(1:end-1, 1:end-1, :) + vu_pp .* u_tmp(2:end, 2:end, :);
    ax_v(1:end-1, 2:end, :) = ax_v(1:end-1, 2:end, :) + vu_pm .* u_tmp(2:end, 1:end-1, :);
    ax_v(2:end, 1:end-1, :) = ax_v(2:end, 1:end-1, :) + vu_mp .* u_tmp(1:end-1, 2:end, :);
    ax_v(2:end, 2:end, :) = ax_v(2:end, 2:end, :) + vu_mm .* u_tmp(1:end-1, 1:end-1, :);
    ax_u = reshape(ax_u, (Nx-1)*(Ny-1)*Ns, 1);
    ax_v = reshape(ax_v, (Nx-1)*(Ny-1)*Ns, 1);
    ax = [ax_u; ax_v];
end
function Convec_Cor
    u_sg = (u_star(2:end,:,:) + u_star(1:end-1,:,:) + u_e(2:end,:,:) + u_e(1:end-1,:,:)) ./ 2;
    v_sg = (v_star(:,2:end,:) + v_star(:,1:end-1,:) + v_e(:,2:end,:) + v_e(:,1:end-1,:)) ./ 2;
    diff_Fx = getdFx(rho, u_sg);
    diff_Fy = getdFy(rho, v_sg);
    F_total = dt/dx * diff_Fx + dt/dy * diff_Fy;
    rhoGGG = gpc_growth(rho,K);
    rho(2:end-1,2:end-1,:) = (rho(2:end-1,2:end-1,:) - F_total) + dt * rhoGGG(2:end-1,2:end-1,:);
    p = rho.^(m-1) * m / (m-1);
    u = -(p(3:end,2:end-1,:) - p(1:end-2,2:end-1,:)) ./ (2*dx);
    v = -(p(2:end-1,3:end,:) - p(2:end-1,1:end-2,:)) ./ (2*dy);
end
end
function rst = getdFx(rho, u_star)
    rho = rho(:, 2:end-1, :);
    diff_l = (rho(2:end-1, :, :) - rho(1:end-2, :, :));
    diff_r = (rho(3:end, :, :) - rho(2:end-1, :, :));
    parial_rho = min(diff_l, diff_r) .* (diff_l > 0) .* (diff_r > 0) + max(diff_l, diff_r) .* (diff_l < 0) .* (diff_r < 0);
    rho_l = rho(2:end-2, :, :) + 1/2 * parial_rho(1:end-1, :, :);
    rho_r = rho(3:end-1, :, :) - 1/2 * parial_rho(2:end, :, :);
F = 1/2 * (rho_l .* u_star + rho_r .* u_star - abs(u_star) .* (rho_r - rho_l));
%Add zero rows
   F(end+1,:,:)=0;F(end+1,:,:)=0; F(2:end,:,:)=F(1:end-1,:,:); F(1,:,:)=0;
    rst = F(2:end,:,:) - F(1:end-1,:,:);
end
function rst = getdFy(rho, v_star)
    rho = rho(2:end-1, :, :);
    diff_l = (rho(:, 2:end-1, :) - rho(:, 1:end-2, :));
    diff_r = (rho(:, 3:end, :) - rho(:, 2:end-1, :));
    parial_rho = min(diff_l, diff_r) .* (diff_l > 0) .* (diff_r > 0) + max(diff_l, diff_r) .* (diff_l < 0) .* (diff_r < 0);
    rho_l = rho(:, 2:end-2, :) + 1/2 * parial_rho(:, 1:end-1, :);
    rho_r = rho(:, 3:end-1, :) - 1/2 * parial_rho(:, 2:end, :);
F = 1/2 * (rho_l .* v_star + rho_r .* v_star - abs(v_star) .* (rho_r - rho_l));
F(:,end+1, :)=0;F(:,end+1, :)=0; F(:,2:end, :)=F(:,1:end-1, :); F(:,1, :)=0;
    rst = F(:, 2:end, :) - F(:, 1:end-1, :);
end
