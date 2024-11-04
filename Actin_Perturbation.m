%% --- load analysis results

clear all
load('ConditionSortedResults_ActinPerturbation')


%% Overview plots for all analyzed conditions

dist_prctl_all = zeros(numConds,1);
S5P_median_all = zeros(numConds,1);
S2P_median_all = zeros(numConds,1);
contactFraction_all = zeros(numConds,1);

dist_prctl_CI = zeros(2,numConds);
S5P_median_CI = zeros(2,numConds);
S2P_median_CI = zeros(2,numConds);
contactFraction_CI = zeros(2,numConds);

dist_prob = cell(1,numConds);

prctile_val = 5;
n_boot = 500;
dist_support = linspace(0,8,500);
dist_bandwidth = 0.25;
dist_threshold = 0.8;

for cc = 1:numConds

    thisCondInd = cc;

    dist_vals = [sortedDistCell{thisCondInd}];
    OP_S5P_vals = [sortedOPIntCell{2}{thisCondInd}];
    OP_S2P_vals = [sortedOPIntCell{1}{thisCondInd}];

    dist_prctl_all(cc) =  prctile(dist_vals,prctile_val);
    S5P_median_all(cc) = median(OP_S5P_vals);
    S2P_median_all(cc) = prctile(OP_S2P_vals,90);

    S5P_median_CI(:,cc) = bootci(n_boot,...
        @(xx)median(xx),OP_S5P_vals);

    S2P_median_CI(:,cc) = bootci(n_boot,...
        @(xx)prctile(xx,90),OP_S2P_vals);

    dist_prctl_CI(:,cc) = bootci(n_boot,...
        @(xx)prctile(xx,prctile_val),dist_vals);

    this_dist_prob = ksdensity(dist_vals,dist_support,...
        'Support','positive','Bandwidth',dist_bandwidth);
    dist_prob{cc} = this_dist_prob;

    contactFraction_all(cc) = mean(dist_vals<=dist_threshold);
    contactFraction_CI(:,cc) = bootci(n_boot,...
        @(xx)mean(xx<=dist_threshold),dist_vals);


end


%% -- Overview figure

gene_names = {...
    'foxd5','klf2b',...
    'zgc:64022','iscub','gadd45ga'};

treatment_names = {...
    'Uninj','BFP','Actin-WT','Actin-R62D',};

gene_cond_inds = {...
    [1,2,4,3]+8,...
    [1,2,4,3]+0,...
    [1,2,4,3]+16,...
    [1,2,4,3]+4,...
    [1,2,4,3]+12,...
    };


figure(1)
clf

figure(2)
clf

S2P_min = 1;
S2P_max = 1.3;
S5P_min = 1;
S5P_max = 1.25;
dist_min = 0;
dist_max = 0.8;

Dist_Delta = zeros(1,5);
S5P_Delta = zeros(1,5);
S2P_Delta = zeros(1,5);

for gg = 1:5

    figure(1)

    subplot(3,5,gg)
    errorbar(1:4,dist_prctl_all(gene_cond_inds{gg}),...
        +dist_prctl_CI(2,gene_cond_inds{gg})...
        -dist_prctl_all(gene_cond_inds{gg})',...
        -dist_prctl_CI(1,gene_cond_inds{gg})...
        +dist_prctl_all(gene_cond_inds{gg})',...
        'ko-','Color',[0.0,0.0,0.0],'LineWidth',1,...
        'MarkerFaceColor',[0,0,0])
    ylabel('Gene-cluser dist. [\mum]','FontWeight','normal')
    set(gca,'XTickLabels',treatment_names,...
        'XTickLabelRotation',90,'YLim',[0,1.5],'XLim',[0.8,4.2],...
        'YTick',0:0.5:1.5,'XTick',1:4)
    title(sortedCondNames(gene_cond_inds{gg}(1)))
    title(gene_names{gg})


    subplot(3,5,5+gg)
    errorbar(1:4,S5P_median_all(gene_cond_inds{gg}),...
        +S5P_median_CI(2,gene_cond_inds{gg})...
        -S5P_median_all(gene_cond_inds{gg})',...
        -S5P_median_CI(1,gene_cond_inds{gg})...
        +S5P_median_all(gene_cond_inds{gg})',...
        'ko-','Color',[1.0,0.0,0.0],'LineWidth',1,...
        'MarkerFaceColor',[1,0,0])
    ylabel('Pol II Ser5P Int. (a.u.)',...
        'FontWeight','normal')
    set(gca,'XTickLabels',treatment_names,...
        'XTickLabelRotation',90,'YLim',[0.95,1.1],'XLim',[0.8,4.2],...
        'YTick',0.9:0.1:1.3,'XTick',1:5)

    subplot(3,5,10+gg)
    errorbar(1:4,S2P_median_all(gene_cond_inds{gg}),...
        +S2P_median_CI(2,gene_cond_inds{gg})...
        -S2P_median_all(gene_cond_inds{gg})',...
        -S2P_median_CI(1,gene_cond_inds{gg})...
        +S2P_median_all(gene_cond_inds{gg})',...
        'ko-','Color',[0.5,0.5,0.5],'LineWidth',1,...
        'MarkerFaceColor',[0.5,0.5,0.5])
    ylabel('Pol II Ser2P Int. (a.u.)','FontWeight','normal')
    set(gca,'XTickLabels',treatment_names,...
        'XTickLabelRotation',90,'YLim',[0.95,1.5],'XLim',[0.8,4.2],...
        'YTick',1:0.05:1.2,'XTick',1:4)

    Dist_Delta(gg) = ...
        dist_prctl_all(gene_cond_inds{gg}(3))...
        -dist_prctl_all(gene_cond_inds{gg}(4));

    S5P_Delta(gg) = ...
        S5P_median_all(gene_cond_inds{gg}(3))...
        -S5P_median_all(gene_cond_inds{gg}(4));

    S2P_Delta(gg) = ...
        S2P_median_all(gene_cond_inds{gg}(3))...
        -S2P_median_all(gene_cond_inds{gg}(2));

    figure(2)
    subplot(4,5,gg)
    plot([1,1].*dist_threshold,...
        [0,0.5],'k--')
    hold on
    lineStyles = {'y-','k-','b:','r--'};
    for kk = [1,2,4,3]
        plot(dist_support,dist_prob{gene_cond_inds{gg}(kk)},...
            lineStyles{kk},'LineWidth',1)
        hold on
    end
    
    xlabel('Gene-cluster dist. d [\mum]')
    ylabel('p(d)')
    title(contactFraction_all(gene_cond_inds{gg}(1)))
    title(sortedCondNames{gene_cond_inds{gg}(1)})

    subplot(4,5,gg+5)
    plot(sortedIntCell{1}{gene_cond_inds{gg}(2)},...
        sortedIntCell{2}{gene_cond_inds{gg}(2)},...
        'k.')
    set(gca,'YLim',[0,4],'XLim',[0,4])

    subplot(4,5,gg+10)
    plot(sortedIntCell{1}{gene_cond_inds{gg}(4)},...
        sortedIntCell{2}{gene_cond_inds{gg}(4)},...
        'k.')
    set(gca,'YLim',[0,4],'XLim',[0,4])

    subplot(4,5,gg+15)
    plot(sortedIntCell{1}{gene_cond_inds{gg}(3)},...
        sortedIntCell{2}{gene_cond_inds{gg}(3)},...
        'k.')
    set(gca,'YLim',[0,4],'XLim',[0,4])


end

figure(3)
clf

subplot(3,1,1)
xLocVals = [1:4,5.5];
plot(xLocVals,Dist_Delta,'ko',...
    'MarkerFaceColor',[0,0,0])
xlabel('')
ylabel('\Deltad (WT-R62D) [\mum]')
set(gca,'XTick',xLocVals,'XTickLabel',gene_names,...
    'XLim',[0.5,6])
hold on
plot([0.5,6],[0,0],'k-')

subplot(3,1,2)
plot(xLocVals,S5P_Delta,'ko',...
    'MarkerFaceColor',[0,0,0])
xlabel('')
ylabel('\DeltaPol II S5P (WT-R62D) [\mum]')
set(gca,'XTick',xLocVals,'XTickLabel',gene_names,...
    'XLim',[0.5,6])
hold on
plot([0.5,6],[0,0],'k-')

subplot(3,1,3)
plot(xLocVals,S2P_Delta,'ko',...
    'MarkerFaceColor',[0,0,0])
xlabel('')
ylabel('\DeltaPol II S2P (WT-BFP) [\mum]')
set(gca,'XTick',xLocVals,'XTickLabel',gene_names,...
    'XLim',[0.5,6])
hold on
plot([0.5,6],[0,0],'k-')

%% -- generate output of counts of nuclei and observations

for cc = 1:numConds

    disp(sprintf('%s, %d nuclei, %d genes',...
        sortedCondNames{cc},...
        sortedNumNuclei(cc),numel(sortedDistCell{cc})))

end