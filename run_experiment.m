% 完整实验入口：读取原设置、调用求解器、重构并保存结果。
function results=run_experiment(mode)
%RUN_EXPERIMENT Run the original gPC-SG solver with the archived experiment settings.
% mode='short' checks the first saved time; 'reference' also recomputes all SC.
% 'reference_short' computes both methods at the first saved time.
% Archived result arrays are used only for comparison, not as computed outputs.
case_name='concave';
if nargin<1,mode='reference';end
assert(any(strcmp(mode,{'cached_reference','short','reference','reference_short'})),'Unknown mode');
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'solver'));out=fullfile(root,'results');if ~exist(out,'dir'),mkdir(out);end
S=load(fullfile(root,'data',[case_name '_compute.mat']));
assert(S.K==4 && S.Nz==numel(S.z) && isequal(size(S.X),[45 45]),'Unexpected archived settings');
T=S.T;if any(strcmp(mode,{'short','reference_short'})),T=T(1);end
% Recompute the initial gPC projection from the original sampled initial fields.
initial=reshape(reshape(S.rho_samples,[],numel(S.z))*(S.psi.*S.wz),size(S.rho_gPC));
assert(max(abs(initial(:)-S.rho_gPC(:)))<1e-12,'Initial projection differs from saved settings');
results=struct('T',T,'X',S.X,'Y',S.Y,'m',80,'K',S.K,'Nz',S.Nz,'dt',.001,'dx',.1,'mode',mode);
results.sc_mean=S.mean_results_collo(:,:,1:numel(T));results.sc_std=S.std_results_collo(:,:,1:numel(T));
for i=1:numel(T)
 timer=tic;
 A=solve_gpc_sg(T(i),80,initial);
 results.sg_seconds(i)=toc(timer);results.sg_mean(:,:,i)=A(:,:,1);
 results.sg_std(:,:,i)=sqrt(sum(A(:,:,2:end).^2,3));
 d=results.sg_mean(:,:,i)-S.mean_sg_results3(:,:,i);e=results.sg_std(:,:,i)-S.std_sg__results3(:,:,i);
 results.saved_mean_max_difference(i)=max(abs(d(:)));results.saved_std_max_difference(i)=max(abs(e(:)));
 if any(strcmp(mode,{'reference','reference_short'}))
  Q=zeros(size(S.rho_samples));
  for j=1:numel(S.z),Q(:,:,j)=solve_sc(j,T(i),80,S.rho_samples,S.z);end
  V=reshape(Q,[],numel(S.z));mu=V*S.wz;
  results.sc_mean(:,:,i)=reshape(mu,size(S.X));results.sc_std(:,:,i)=reshape(sqrt((V-mu).^2*S.wz),size(S.X));
  d=results.sc_mean(:,:,i)-S.mean_results_collo(:,:,i);e=results.sc_std(:,:,i)-S.std_results_collo(:,:,i);
  results.sc_saved_mean_max_difference(i)=max(abs(d(:)));results.sc_saved_std_max_difference(i)=max(abs(e(:)));
 end
 fprintf('%s T=%g: SG computed; saved mean/std max difference %.4g / %.4g\n',case_name,T(i),results.saved_mean_max_difference(i),results.saved_std_max_difference(i));
 save(fullfile(out,[case_name '_computed_' mode '.mat']),'-struct','results','-v7');
end
assert(all(isfinite(results.sg_mean(:)))&&all(isfinite(results.sg_std(:))),'Nonfinite SG result');
f=figure('visible','off','Color','w','position',[100 100 1200 450],'paperpositionmode','auto');
row=20;colors=lines(numel(T));
for q=1:2
 subplot(1,2,q);hold on;
 for i=1:numel(T)
  if q==1,a=results.sg_mean(:,:,i);b=results.sc_mean(:,:,i);else,a=results.sg_std(:,:,i);b=results.sc_std(:,:,i);end
  plot(S.X(1,:),b(row,:),'-','Color',colors(i,:),'DisplayName',sprintf('SC T=%g',T(i)));
  plot(S.X(1,:),a(row,:),'o','Color',colors(i,:),'DisplayName',sprintf('Computed SG T=%g',T(i)));
 end
 xlabel('x');if q==1,ylabel('Mean');else,ylabel('Standard deviation');end;grid on;legend('show');
end
print(f,fullfile(out,[case_name '_computed_' mode '.png']),'-dpng','-r150');close(f);
end
