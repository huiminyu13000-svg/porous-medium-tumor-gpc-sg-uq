% 直接使用原始 concave 数据画图，不运行求解器。
function result=plot_saved
%PLOT_SAVED_GPC Reproduce the saved paper figures without running a solver.
% Curve arrays and slice index come directly from the original workspace.
% The layout and line/marker conventions follow the original test1111 script.
case_name='concave';
root=fileparts(mfilename('fullpath'));out=fullfile(root,'results');
if ~exist(out,'dir'),mkdir(out);end
fields={'X','Y','T','x_vals','x_index','sg_means','sg_stds','coll_means','coll_stds', ...
 'mean_sg_results3','std_sg__results3','mean_results_collo','std_results_collo','mean_errors','std_errors'};
result=load(fullfile(root,'data',['gpc_' case_name '_paper.mat']),fields{:});
for k=1:numel(fields),assert(isfield(result,fields{k}),'Missing saved field: %s',fields{k});end
S=result;colors=[0 .447 .741;.850 .325 .098;.929 .694 .125;.494 .184 .556;.466 .674 .188;.301 .745 .933];
f=figure('visible','off','Color','w','position',[100 100 1200 450],'paperpositionmode','auto');
for q=1:2
 if q==1,ref=S.coll_means;sg=S.sg_means;label='E(\rho)';else,ref=S.coll_stds;sg=S.sg_stds;label='Std(\rho)';end
 subplot(1,2,q);hold on;h=zeros(1,6);labels=cell(1,6);markers={'*','x','^'};
 for k=1:3
  h(k)=plot(S.x_vals,ref(:,k),'-','Color',colors(k,:),'LineWidth',1.0);
  labels{k}=sprintf('Ref %s T=%g',label,S.T(k));
 end
 for k=1:3
  h(k+3)=plot(S.x_vals,sg(:,k),markers{k},'Color',colors(k+3,:),'LineWidth',1.0);
  labels{k+3}=sprintf('%s T=%g',label,S.T(k));
  assert(isequal(get(h(k),'YData')(:),ref(:,k))&&isequal(get(h(k+3),'YData')(:),sg(:,k)));
 end
 xlabel('x');ylabel(label);title([label ' Comparison along x']);grid on;
 legend(h,labels,'Location','northeast');
end
print(f,fullfile(out,[case_name '_slice.png']),'-dpng','-r180');close(f);
% Original statistics and difference fields at the second saved time.
for q=1:2
 if q==1,sg=S.mean_sg_results3(:,:,2);ref=S.mean_results_collo(:,:,2);label='E(\rho)';tag='mean';
 else,sg=S.std_sg__results3(:,:,2);ref=S.std_results_collo(:,:,2);label='Std(\rho)';tag='std';end
 f=figure('visible','off','Color','w','position',[100 100 1000 750],'paperpositionmode','auto');
 subplot(2,2,1);imagesc(sg);axis image;colorbar;title(label);xlabel('x');ylabel('y');
 subplot(2,2,2);surf(sg,'EdgeColor','none');view(3);title([label ' (3D)']);xlabel('x');ylabel('y');
 subplot(2,2,3);imagesc(sg-ref);axis image;colorbar;title(['Error for ' label]);xlabel('x');ylabel('y');
 subplot(2,2,4);surf(sg-ref,'EdgeColor','none');view(3);title(['Error for ' label ' (3D)']);xlabel('x');ylabel('y');
 print(f,fullfile(out,[case_name '_' tag '.png']),'-dpng','-r150');close(f);
end
% Plot the saved error statistics with the original log10 presentation.
f=figure('visible','off','Color','w','paperpositionmode','auto');
plot(S.T,log10(S.mean_errors),'-o',S.T,log10(S.std_errors),'-*','LineWidth',1.5);
xlabel('Time');ylabel('Saved error (log10)');legend({'Mean','Standard deviation'});grid on;
print(f,fullfile(out,[case_name '_saved_error.png']),'-dpng','-r150');close(f);
save(fullfile(out,[case_name '_replayed_fields.mat']),'-struct','result','-v7');
fprintf('Replayed original %s dataset; saved slice index=%d.\n',case_name,S.x_index);
end
