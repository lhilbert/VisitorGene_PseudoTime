%% --- load analysis results

clear all
load('ActinOverexpression_Results')

dist_prctl_all = zeros(numConds,1);
S5P_mean_all = zeros(numConds,1);
S2P_mean_all = zeros(numConds,1);

dist_prctl_CI = zeros(2,numConds);
S5P_mean_CI = zeros(2,numConds);
S2P_mean_CI = zeros(2,numConds);

prctile_val = 5;
n_boot = 500;

for cc = 1:numConds

    thisCondInd = cc;

    dist_vals = [sortedDistCell{thisCondInd}];
    OP_S5P_vals = [sortedOPIntCell{2}{thisCondInd}];
    OP_S2P_vals = [sortedOPIntCell{1}{thisCondInd}];
    
    dist_prctl_all(cc) =  prctile(dist_vals,prctile_val);
    S5P_mean_all(cc) = mean(OP_S5P_vals);
    S2P_mean_all(cc) = mean(OP_S2P_vals);

    S5P_mean_CI(:,cc) = bootci(n_boot,...
        @(xx)mean(xx),OP_S5P_vals);

    S2P_mean_CI(:,cc) = bootci(n_boot,...
        @(xx)mean(xx),OP_S2P_vals);

    dist_prctl_CI(:,cc) = bootci(n_boot,...
        @(xx)prctile(xx,prctile_val),dist_vals);

    dist_prctl_CI(1,cc) = prctile(dist_vals,2);
    dist_prctl_CI(2,cc) = prctile(dist_vals,8);

end

figure(1)
clf

S2P_min = 1;
S2P_max = 1.3;
S5P_min = 1;
S5P_max = 1.25;
dist_min = 0;
dist_max = 0.8;


%% -- Overview figure



gene_names = {...
    'klf2b','iscub','foxd5',...
    'gadd45ga','zgc::64022'};

treatment_names = {...
    'Uninj','BFP','R62D','WT'};

gene_trajectory_inds = {...
    [1,6,11,16]+0,...
    [1,6,11,16]+1,...
    [1,6,11,16]+2,...
    [1,6,11,16]+3,...
    [1,6,11,16]+4};

for gg = 1:5

subplot(3,5,gg)
errorbar(1:5,S5P_mean_all(gene_trajectory_inds{gg}),...
    +S5P_mean_CI(2,gene_trajectory_inds{gg})...
    -S5P_mean_all(gene_trajectory_inds{gg})',...
    -S5P_mean_CI(1,gene_trajectory_inds{gg})...
    +S5P_mean_all(gene_trajectory_inds{gg})',...
    'ko-','Color',[1.0,0.0,0.0],'LineWidth',1,...
    'MarkerFaceColor',[1,0,0])
title('Pol II Ser5P Int. (a.u.)',...
    'FontWeight','normal')
set(gca,'XTickLabels','','YLim',[0.95,1.3],'XLim',[0.8,5.2],...
    'YTick',0.9:0.1:1.3,'XTick',1:5)

ylabel(gene_names(gg))


subplot(3,5,5+gg)
errorbar(1:5,S2P_mean_all(gene_trajectory_inds{gg}),...
    +S2P_mean_CI(2,gene_trajectory_inds{gg})...
    -S2P_mean_all(gene_trajectory_inds{gg})',...
    -S2P_mean_CI(1,gene_trajectory_inds{gg})...
    +S2P_mean_all(gene_trajectory_inds{gg})',...
    'ko-','Color',[0.5,0.5,0.5],'LineWidth',1,...
    'MarkerFaceColor',[0.5,0.5,0.5])
title('Pol II Ser2P Int. (a.u.)','FontWeight','normal')
set(gca,'XTickLabels','','YLim',[0.95,1.2],'XLim',[0.8,5.2],...
    'YTick',1:0.05:1.2,'XTick',1:5)


subplot(3,5,5+gg)
errorbar(1:5,dist_prctl_all(gene_trajectory_inds{gg}),...
    +dist_prctl_CI(2,gene_trajectory_inds{gg})...
    -dist_prctl_all(gene_trajectory_inds{gg})',...
    -dist_prctl_CI(1,gene_trajectory_inds{gg})...
    +dist_prctl_all(gene_trajectory_inds{gg})',...
    'ko-','Color',[0.0,0.0,0.0],'LineWidth',1,...
    'MarkerFaceColor',[0,0,0])
title('Gene-cluser dist. [\mum]','FontWeight','normal')
set(gca,'XTickLabels',treatment_names,...
    'XTickLabelRotation',-0,'YLim',[0,1.5],'XLim',[0.8,5.2],...
    'YTick',0:0.5:1.5,'XTick',1:5)


%% -- generate output of counts of nuclei and observations

for cc = 1:numConds

    disp(sprintf('%s, %d nuclei, %d genes',...
        sortedCondNames{cc},...
        sortedNumNuclei(cc),numel(sortedDistCell{cc})))

end