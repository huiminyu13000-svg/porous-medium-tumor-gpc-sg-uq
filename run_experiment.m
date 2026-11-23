function results = run_experiment(mode)
%RUN_EXPERIMENT Generate initial data, solve SG/SC, and plot the statistics.
% run_experiment computes all output times; 'short' computes T=0.01.
if nargin < 1, mode = 'full'; end
assert(any(strcmp(mode,{'full','short'})), 'Use full or short.');
root = fileparts(mfilename('fullpath')); addpath(root,fullfile(root,'solver'));
cfg = experiment_config(); if strcmp(mode,'short'), cfg.times = cfg.times(1); end
out = fullfile(root,'results'); if ~exist(out,'dir'), mkdir(out); end
grid = cfg.domain(1):cfg.dx:cfg.domain(2); [X,Y] = meshgrid(grid);
[z,w] = lgwt(cfg.n_collocation,-1,1); w = w(:)/2;
psi = legendre_basis(z,cfg.num_modes);
samples = zeros([size(X),numel(z)]);
for j = 1:numel(z), samples(:,:,j) = initial_density(X,Y,z(j)); end
initial = reshape(reshape(samples,[],numel(z))*(psi.*w),[size(X),cfg.num_modes]);
results = struct('settings',cfg,'X',X,'Y',Y,'nodes',z,'weights',w, ...
    'initial_samples',samples,'initial_coefficients',initial);
for j = 1:numel(cfg.times)
    T = cfg.times(j); timer = tic;
    [coefficients,~,info] = solve_gpc_sg(T,cfg.m,initial,cfg);
    results.sg_seconds(j) = toc(timer); results.sg_steps(j) = info.steps;
    results.sg_mean(:,:,j) = coefficients(:,:,1);
    results.sg_std(:,:,j) = sqrt(sum(coefficients(:,:,2:end).^2,3));
    timer = tic; values = zeros(numel(X),numel(z));
    for k = 1:numel(z)
        field = solve_sc(k,T,cfg.m,samples,z,cfg); values(:,k) = field(:);
    end
    results.sc_seconds(j) = toc(timer); mu = values*w;
    results.sc_mean(:,:,j) = reshape(mu,size(X));
    results.sc_std(:,:,j) = reshape(sqrt(max(0,(values-mu).^2*w)),size(X));
    results.mean_l2(j) = cfg.dx*norm(reshape(results.sg_mean(:,:,j)-results.sc_mean(:,:,j),[],1));
    results.std_l2(j) = cfg.dx*norm(reshape(results.sg_std(:,:,j)-results.sc_std(:,:,j),[],1));
    fprintf('T=%g: SG and %d SC nodes completed; mean/std L2 = %.6g / %.6g\n', ...
        T,numel(z),results.mean_l2(j),results.std_l2(j));
end
assert(all(isfinite([results.sg_mean(:);results.sg_std(:);results.sc_mean(:);results.sc_std(:)])), ...
    'Nonfinite computed statistics.');
results.environment = struct('version',version,'computer',computer);
save(fullfile(out,['concave_' mode '.mat']),'-struct','results','-v7');
plot_results(results,out);
end

function plot_results(R,out)
% All figures are generated from the calculation above.
[~,row] = min(abs(R.Y(:,1)-R.settings.slice_coordinate));
f = figure('Color','w','Position',[100 100 1100 430]); colors = lines(numel(R.settings.times));
for panel = 1:2
    subplot(1,2,panel); hold on; labels = {};
    for j = 1:numel(R.settings.times)
        if panel==1, a=R.sg_mean(:,:,j);b=R.sc_mean(:,:,j); else,a=R.sg_std(:,:,j);b=R.sc_std(:,:,j);end
        plot(R.X(1,:),b(row,:),'-','Color',colors(j,:),'LineWidth',1.3);
        plot(R.X(1,:),a(row,:),'o','Color',colors(j,:),'MarkerSize',4);
        labels(end+1:end+2) = {sprintf('SC T=%g',R.settings.times(j)),sprintf('SG T=%g',R.settings.times(j))};
    end
    xlabel('x'); if panel==1,ylabel('Mean');else,ylabel('Standard deviation');end
    title(sprintf('y = %.1f',R.Y(row,1))); legend(labels,'Location','best'); grid on;
end
print(f,fullfile(out,'concave_slices.png'),'-dpng','-r160');
[~,j] = min(abs(R.settings.times-1));
f = figure('Color','w','Position',[100 100 900 700]);
fields = {R.sg_mean(:,:,j),R.sc_mean(:,:,j),R.sg_std(:,:,j),R.sc_std(:,:,j)};
labels = {'SG mean','SC mean','SG standard deviation','SC standard deviation'};
for k = 1:4
    subplot(2,2,k); imagesc(R.X(1,:),R.Y(:,1),fields{k}); axis image; set(gca,'YDir','normal');
    pair = 2*floor((k-1)/2)+1;
    limits = [min([fields{pair}(:);fields{pair+1}(:)]),max([fields{pair}(:);fields{pair+1}(:)])];
    if limits(2)<=limits(1),limits(2)=limits(1)+1;end
    caxis(limits);
    title(sprintf('%s, T=%g',labels{k},R.settings.times(j))); xlabel('x');ylabel('y');colorbar;
end
print(f,fullfile(out,'concave_fields.png'),'-dpng','-r160');
end
