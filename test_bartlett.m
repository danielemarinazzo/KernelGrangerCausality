clear;clc
%%% this script compares three tests for the significance of the kernel GC link 2->1:
%%% 'standard' (corr, T-2 dof), 'bartlett' (corrected for autocorrelation and
%%% conditioning dimension, per direction, corr_pval_bartlett) and 'lambda'
%%% (modified Lambda-test of the whole link, lambda_test_bartlett)
%%% (Cliff et al., Phys. Rev. Research 3, 013145, 2021)
%%% null: two independent AR(1) processes (0.3 and -0.8, as in the paper),
%%% low-pass Butterworth filtered (cutoff pi/2) with increasing order
%%% power: quadratic coupling y->x of increasing strength
%%% figures and results are saved in the folder figures_bartlett
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=200;              %number of realizations (null)
Rp=100;             %number of realizations (power)
T=512;              %number of points of simulated time series
m=2;                %order of the model
forders=[0 4 8 12 16];          %orders of the Butterworth filter 
kern={'p',1;'p',2;'g',2};       %kernels {type, par} ('p' polynomial, 'g' gaussian)
gs=[0 0.05 0.1 0.15 0.2 0.3];   %coupling strengths y->x
fpow=[0 8];         %orders of the Butterworth filter (power)
alpha=0.05;         %significance level
fpp=8;              %filter order shown in the P-P plots
quick=true;         %true: short run (~5 min), false: full run (~1 h)
if quick
    R=30;
    Rp=15;
    forders=[0 8 16];
    gs=[0 0.1 0.2 0.3];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outdir=fullfile(fileparts(mfilename('fullpath')),'figures_bartlett');
if quick, outdir=[outdir '_quick']; end
if ~exist(outdir,'dir'), mkdir(outdir); end
addpath(fileparts(mfilename('fullpath')));    %toolbox on the path
olddir=cd(tempdir);     %filtro.m saves test_fast.mat/test_chol.mat in the current folder
rng(1)
nk=size(kern,1);
nf=length(forders);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% null %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=nan(R,nf,nk,3);   %p of the link 2->1 (1 standard, 2 bartlett, 3 lambda)
eta=nan(R,nf);      %effective sample size / T of the raw pair
for io=1:nf
    [b,a]=filtcoef(forders(io));
    for r=1:R
        z=zeros(T+300,2);
        e=randn(T+300,2);
        for t=2:T+300
            z(t,1)=0.3*z(t-1,1)+e(t,1);
            z(t,2)=-0.8*z(t-1,2)+e(t,2);
        end
        z=filter(b,a,z);
        Y=z(end-T+1:end,:);         %data matrix [n_points n_variables]
        [~,~,nu]=corr_pval_bartlett(Y(:,1),Y(:,2),0);
        eta(r,io)=(nu+2)/T;
        for ik=1:nk
            [~,~,rs,ps]=causality(Y,kern{ik,1},kern{ik,2},m,'Bonferroni','standard');
            [~,~,rb,pb]=causality(Y,kern{ik,1},kern{ik,2},m,'Bonferroni','bartlett');
            P(r,io,ik,1)=padj(rs,ps);
            P(r,io,ik,2)=padj(rb,pb);
            [~,~,~,pl]=causality(Y,kern{ik,1},kern{ik,2},m,'Bonferroni','lambda');
            P(r,io,ik,3)=pl(2,1);
        end
    end
    disp(['null, filter order ' num2str(forders(io)) ' done']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PW=nan(Rp,length(gs),length(fpow),3);
for io=1:length(fpow)
    [b,a]=filtcoef(fpow(io));
    for ig=1:length(gs)
        for r=1:Rp
            e=randn(T+300,2);
            y=filter(1,[1 -0.5],e(:,2));
            x=zeros(T+300,1);
            for t=2:T+300
                x(t)=0.4*x(t-1)+gs(ig)*y(t-1)^2+e(t,1);
            end
            z=filter(b,a,[x y]);
            Y=z(end-T+1:end,:);
            [~,~,rs,ps]=causality(Y,'p',2,m,'Bonferroni','standard');
            [~,~,rb,pb]=causality(Y,'p',2,m,'Bonferroni','bartlett');
            PW(r,ig,io,1)=padj(rs,ps);
            PW(r,ig,io,2)=padj(rb,pb);
            [~,~,~,pl]=causality(Y,'p',2,m,'Bonferroni','lambda');
            PW(r,ig,io,3)=pl(2,1);
        end
    end
    disp(['power, filter order ' num2str(fpow(io)) ' done']);
end
cd(olddir);
save(fullfile(outdir,'results_bartlett.mat'),'P','PW','eta','R','Rp','T','m','forders','kern','gs','fpow','alpha');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col=[0.35 0.35 0.35; 0.00 0.45 0.74; 0.85 0.33 0.10];
lab={'Standard test','Bartlett-corrected test','Modified \Lambda-test'};
ktit=cell(nk,1);
for ik=1:nk
    if kern{ik,1}=='p'
        ktit{ik}=sprintf('Polynomial kernel, p = %g',kern{ik,2});
    else
        ktit{ik}=sprintf('Gaussian kernel, \\sigma = %g',kern{ik,2});
    end
end
ci=1.96*sqrt(alpha*(1-alpha)/R);    %95% CI of the FPR under the null
%%% Fig. 1: false positive rate at alpha vs filter order
f1=newfig([100 100 1100 330]);
tl=tiledlayout(1,nk,'TileSpacing','compact','Padding','compact');
for ik=1:nk
    nexttile;hold on
    fill([forders fliplr(forders)],[repmat(alpha-ci,1,nf) repmat(alpha+ci,1,nf)],[0.9 0.9 0.9],'EdgeColor','none');
    yline(alpha,':k');
    for it=1:3
        plot(forders,mean(P(:,:,ik,it)<alpha,1),'-o','Color',col(it,:),'MarkerFaceColor',col(it,:),'LineWidth',1.5);
    end
    title(ktit{ik},'FontWeight','normal');
    xlabel('Filter order');xlim([forders(1) forders(end)]);ylim([0 max(ylim)]);
    if ik==1, ylabel(sprintf('False positive rate (\\alpha = %.2f)',alpha)); end
    box off;set(gca,'TickDir','out','FontSize',10);
end
lg=legend([{'95% CI',''} lab],'Orientation','horizontal');lg.Layout.Tile='south';legend boxoff
savefig2(f1,outdir,'fig1_fpr_vs_filter');
%%% Fig. 2: P-P plots (actual vs nominal FPR) at filter order fpp
io=find(forders==fpp,1);
al=linspace(0,1,201);
f2=newfig([100 100 1100 360]);
tl=tiledlayout(1,nk,'TileSpacing','compact','Padding','compact');
for ik=1:nk
    nexttile;hold on
    plot([0 1],[0 1],':k');
    for it=1:3
        plot(al,mean(P(:,io,ik,it)<al,1),'Color',col(it,:),'LineWidth',1.5);
    end
    axis square;xlim([0 1]);ylim([0 1]);
    title(ktit{ik},'FontWeight','normal');
    xlabel('Nominal \alpha');
    if ik==1, ylabel('Actual FPR'); end
    box off;set(gca,'TickDir','out','FontSize',10);
end
lg=legend([{''} lab],'Orientation','horizontal');lg.Layout.Tile='south';legend boxoff
title(tl,sprintf('Filter order %d',fpp),'FontSize',11);
savefig2(f2,outdir,'fig2_pp_plots');
%%% Fig. 3: effective sample size (a) and power (b)
f3=newfig([100 100 800 330]);
tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;hold on
mu=mean(eta,1);sd=std(eta,0,1);
fill([forders fliplr(forders)],[mu-sd fliplr(mu+sd)],col(2,:),'FaceAlpha',0.2,'EdgeColor','none');
plot(forders,mu,'-o','Color',col(2,:),'MarkerFaceColor',col(2,:),'LineWidth',1.5);
yline(1,':k');
xlabel('Filter order');ylabel('Effective sample size \eta / T');xlim([forders(1) forders(end)]);
title('(a) Null pair','FontWeight','normal');
box off;set(gca,'TickDir','out','FontSize',10);
nexttile;hold on
ls={'-','--'};      %line style per filter order
h=gobjects(0);hl={};
for io=1:length(fpow)
    for it=1:3
        h(end+1)=plot(gs,mean(PW(:,:,io,it)<alpha,1),ls{io},'Color',col(it,:),'Marker','o','MarkerFaceColor',col(it,:),'LineWidth',1.5);
        hl{end+1}=sprintf('%s, filter %d',lab{it},fpow(io));
    end
end
yline(alpha,':k');
xlabel('Coupling strength y \rightarrow x');ylabel(sprintf('Detection rate (\\alpha = %.2f)',alpha));
title('(b) Power, polynomial kernel p = 2','FontWeight','normal');ylim([0 1]);
legend(h,hl,'Location','southeast','FontSize',8);legend boxoff
box off;set(gca,'TickDir','out','FontSize',10);
savefig2(f3,outdir,'fig3_eta_power');

function p=padj(rr,pp)
%Bonferroni-adjusted p-value of the link 2->1 over the new directions (rr~=0)
r=squeeze(rr(2,1,:));
p=squeeze(pp(2,1,:));
k=r~=0;
if any(k)
    p=min(1,min(p(k))*sum(k));
else
    p=1;
end
end

function [b,a]=filtcoef(fo)
%Butterworth low-pass of order fo, cutoff pi/2 (fo=0: no filter)
if fo>0
    [b,a]=butter(fo,0.5);
else
    b=1;a=1;
end
end

function savefig2(f,outdir,name)
%export as .png (300 dpi), vector .pdf and .fig
exportgraphics(f,fullfile(outdir,[name '.png']),'Resolution',300);
exportgraphics(f,fullfile(outdir,[name '.pdf']),'ContentType','vector');
savefig(f,fullfile(outdir,[name '.fig']));
end

function f=newfig(pos)
%white figure with light theme (dark theme otherwise used in batch mode, R2025a+)
f=figure('Color','w','Position',pos);
try
    theme(f,'light');
end
end
