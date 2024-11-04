%% --- load analysis results

clear all
load('ConditionSortedResults_ActinPerturbation_withQC')


%% Overview plots for all analyzed conditions



dist_prctl_all = zeros(numConds,1);
S5P_mean_all = zeros(numConds,1);
S2P_mean_all = zeros(numConds,1);
contactFraction_all = zeros(numConds,1);
nn_observations = zeros(numConds,1);

dist_prctl_CI = zeros(2,numConds);
S5P_mean_CI = zeros(2,numConds);
S2P_mean_CI = zeros(2,numConds);
contactFraction_CI = zeros(2,numConds);

prctile_val = 5;
n_boot = 1000;

dist_support = linspace(0,8,500);
dist_bandwidth = 0.15;
dist_threshold = 0.7;
dist_prob = cell(1,numConds);

NN_dist_prob = cell(1,numConds);

nucVol_support = linspace(0,1000,500);
nucVol_bandwidth = 50;
nucVol_prob = cell(1,numConds);

for cc = 1:numConds

    thisCondInd = cc;

    nucVol_vals = [sortedNucVol{thisCondInd}];
    nucS5P_vals = [sortedNucInt{2}{thisCondInd}];
    nucMeanNNDist_vals = cellfun(@mean,sorted_S5P_NNDist{thisCondInd});
    ClusterNN_vals =vertcat(sorted_S5P_NNDist{thisCondInd}{:});

    dist_vals = [sorted_paired_distCell{thisCondInd}];
    OP_S5P_vals = [sorted_paired_OPIntCell{2}{thisCondInd}];
    OP_S2P_vals = [sorted_paired_OPIntCell{1}{thisCondInd}];

    nn_observations(cc) = numel(dist_vals);

    S2P_prctile = 50;

    dist_prctl_all(cc) =  prctile(dist_vals,prctile_val);
    S5P_mean_all(cc) = mean(OP_S5P_vals);
    S2P_mean_all(cc) = mean(OP_S2P_vals);

    S5P_mean_CI(:,cc) = bootci(n_boot,...
        @(xx)mean(xx),OP_S5P_vals);

    S2P_mean_CI(:,cc) = bootci(n_boot,...
        @(xx)mean(xx),OP_S2P_vals);

    dist_prctl_CI(:,cc) = bootci(n_boot,...
        @(xx)prctile(xx,S2P_prctile),dist_vals);

    contactFraction_all(cc) = mean(dist_vals<=dist_threshold);
    contactFraction_CI(:,cc) = bootci(n_boot,...
        @(xx)mean(xx<=dist_threshold),dist_vals);

    this_dist_prob = ksdensity(dist_vals,dist_support,...
        'Support','positive','Bandwidth',dist_bandwidth);
    dist_prob{cc} = this_dist_prob;

    this_dist_prob = ksdensity(ClusterNN_vals,dist_support,...
        'Support','positive','Bandwidth',dist_bandwidth);
    NN_dist_prob{cc} = this_dist_prob;

    this_NucVol_prob = ksdensity(nucVol_vals,nucVol_support,...
        'Support','positive','Bandwidth',nucVol_bandwidth);
    nucVol_prob{cc} = this_NucVol_prob;


end


% -- Overview figure

gene_names = {...
    'iscub','gadd45ga','zgc:64022','foxd5','klf2b',...
    };

treatment_names = {...
    'Uninj','BFP','Actin-WT','Actin-R62D',};

gene_cond_inds = {...
    [3,1,4,2]+8,...
    [3,1,4,2]+4,...
    [3,1,4,2]+16,...
    [3,1,4,2]+0,...
    [3,1,4,2]+12,...
    };

target_inds = [...
    1,2,4,3,...
    1,2,4,3,...
    1,2,4,3,...
    1,2,4,3,...
    1,2,4,3];

numGenes = numel(gene_cond_inds);
numTargets = numel(unique(target_inds));

S2P_min = 1.00;
S2P_max = 1.10;
S5P_min = 1.00;
S5P_max = 1.10;
contact_min = 0;
contact_max = 7;

figure(1)
clf

subplot(2,3,1)

treatmentCondInds = [cellfun(@(elmt)elmt([2:4]),gene_cond_inds,...
    'UniformOutput',false)];
treatmentCondInds = [treatmentCondInds{:}];

legend_handles = zeros(1,3);

for gg = 1:numGenes
    plot(S5P_mean_all(gene_cond_inds{gg}([2,3])),...
        S2P_mean_all(gene_cond_inds{gg}([2,3])),...
        'b-','LineWidth',1.5)
    hold on
    plot(S5P_mean_all(gene_cond_inds{gg}([2,4])),...
        S2P_mean_all(gene_cond_inds{gg}([2,4])),...
        'r--','LineWidth',1.5)
    legend_handles(1) = ...
        plot(S5P_mean_all(gene_cond_inds{gg}([2])),...
        S2P_mean_all(gene_cond_inds{gg}([2])),...
        'ko','MarkerFaceColor',[1,1,1]);
    legend_handles(2) = ...
        plot(S5P_mean_all(gene_cond_inds{gg}([3])),...
        S2P_mean_all(gene_cond_inds{gg}([3])),...
        'ks','MarkerFaceColor',[0,0,1]);
    legend_handles(3) = ...
        plot(S5P_mean_all(gene_cond_inds{gg}([4])),...
        S2P_mean_all(gene_cond_inds{gg}([4])),...
        'r^','MarkerFaceColor',[1,0,0]);
    text(S5P_mean_all(gene_cond_inds{gg}([2])),...
        S2P_mean_all(gene_cond_inds{gg}([2])),...
        gene_names{gg},...
        'Color',[0,0,0])
end
xlabel('Pol II Ser5P Int.')
ylabel('Pol II Ser2P Int.')

set(gca,'XLim',[S5P_min,S5P_max],'YLim',[S2P_min,S2P_max])

title('Shift of transcriptional state ')

set(gca,'Box','on')

legend(legend_handles,treatment_names(2:4),...
    'Location','Northwest')

subplot(2,3,2)

treatmentCondInds = [cellfun(@(elmt)elmt(2:4),gene_cond_inds,...
    'UniformOutput',false)];
treatmentCondInds = [treatmentCondInds{:}];

scatter(S5P_mean_all(treatmentCondInds),...
    S2P_mean_all(treatmentCondInds),...
    70,contactFraction_all(treatmentCondInds).*100,...
    'filled','o')

[Rval,Pval] = corrcoef(...
    [S5P_mean_all(treatmentCondInds),...
    S2P_mean_all(treatmentCondInds),...
    contactFraction_all(treatmentCondInds).*100]);

text(1.005,1.085,[...
    sprintf('Correlation coefficients and P values:\n'), ...
    sprintf('(Ser5P, Ser2P): R=%2.2f, P=%1.1e\n',Rval(1,2),Pval(1,2)),...
    sprintf('(Ser5P, Contact %%): R=%2.2f, P=%1.1e\n',Rval(1,3),Pval(1,3)),...
    sprintf('(Ser2P, Contact %%): R=%2.2f, P=%1.1e\n',Rval(2,3),Pval(2,3))],...
    'Color',[0,0,0])

fprintf('Correlation coefficients and P values:\n')
fprintf('(Ser5P, Ser2P): R=%2.2f, P=%1.1e\n',Rval(1,2),Pval(1,2))
fprintf('(Ser5P, Contact %%): R=%2.2f, P=%1.1e\n',Rval(1,3),Pval(1,3))
fprintf('(Ser2P, Contact %%): R=%2.2f, P=%1.1e\n',Rval(2,3),Pval(2,3))

xlabel('Pol II Ser5P Int.')
ylabel('Pol II Ser2P Int.')
colorbar

title('Contact %')

set(gca,'XLim',[S5P_min,S5P_max],'YLim',[S2P_min,S2P_max])
set(gca,'Box','on')
caxis([contact_min,contact_max])




% -- fit 2D polynomial surface

subplot(2,3,3)

inputMatrix = [S5P_mean_all(treatmentCondInds),...
    S2P_mean_all(treatmentCondInds)];
outputVec = contactFraction_all(treatmentCondInds).*100;

modelfun = @(theta,xx) theta(1)...
    +(xx(:,1)-1).*theta(2)...
    +(xx(:,2)-1).*theta(3);

[theta,~,~,~,RR] = ...
    nlinfit(inputMatrix,outputVec,modelfun,[1,1,1]);

Rsquared = 1 - sum(RR.^2)/sum(((outputVec-mean(outputVec)).^2));

grid_N = 150; %150
grid_vals = zeros(grid_N,grid_N);

S2P_vec = linspace(S2P_min,S2P_max,grid_N);
S5P_vec = linspace(S5P_min,S5P_max,grid_N);

for mm = 1:grid_N
    for nn = 1:grid_N

        grid_vals(mm,nn) = modelfun(theta,[S2P_vec(mm),S5P_vec(nn)]);

    end
end

cla
imagesc(S5P_vec,S2P_vec,grid_vals,[contact_min,contact_max])
set(gca,'YDir','normal')
colorbar
hold on

scatter(S5P_mean_all(treatmentCondInds),...
    S2P_mean_all(treatmentCondInds),...
    40,contactFraction_all(treatmentCondInds).*100,...
    'filled','MarkerEdgeColor','k')

xlabel('Mean S5P Int. (a.u.)')
ylabel('Mean S2P Int. (a.u.)')
c = colorbar;
c.Label.String='Contact %';
colormap((parula))
set(gca,'Box','on')

vector_origin = [1.01,1.07];

plot(vector_origin(1),vector_origin(2),'wo',...
    'MarkerSize',10,...
    'MarkerFaceColor',[1,1,1],...
    'MarkerEdgeColor',[0,0,0])

vector_scaleFact = 0.0004;

plot(vector_origin(1)+[0,theta(3).*vector_scaleFact],...
    vector_origin(2)+[0,theta(2).*vector_scaleFact],...
    'w-','LineWidth',1.5)

plot(vector_origin(1)+[0,theta(2).*vector_scaleFact],...
    vector_origin(2)+[0,-theta(3).*vector_scaleFact],...
    'w--','LineWidth',1.5)

text(1.005,1.09,sprintf(...
    'd=%2.2f%%+%2.2f%%\\cdotI_{S2P}+%2.2f%%\\cdotI_{S5P},\nR^2=%5.5f',...
    theta(1),theta(2),theta(3),Rsquared),...
    'Color',[1,1,1])




subplot(2,2,3)

cla

plot([contact_min,contact_max],[0,0],'k-','LineWidth',1)
hold on

legend_handles = zeros(1,numGenes);

for kk = 1:numGenes
    plotContactFraction = ...
        100.*contactFraction_all(gene_cond_inds{kk}([2,4]));
    legend_handles(kk) = plot(plotContactFraction(1).*[1,1],...
        [0,plotContactFraction(2)-plotContactFraction(1)],...
        'o-','LineWidth',1,'MarkerFaceColor','auto');
    hold on
end

xlabel('Contact %, BFP control')
ylabel('\DeltaContact % (R62D actin - BFP control)')
legend(legend_handles,gene_names)
set(gca,'YLim',[-4,+4],'XLim',[contact_min,contact_max])

title('Over-expression of monomeric actin')


subplot(2,2,4)

cla

plot([contact_min,contact_max],[0,0],'k-','LineWidth',1)
hold on

legend_handles = zeros(1,numGenes);

for kk = 1:numGenes
    plotContactFraction = ...
        100.*contactFraction_all(gene_cond_inds{kk}([2,3]));
    legend_handles(kk) = plot(plotContactFraction(1).*[1,1],...
        [0,plotContactFraction(2)-plotContactFraction(1)],...
        'o-','LineWidth',1,'MarkerFaceColor','auto');
    hold on
end

xlabel('Contact %, BFP control')
ylabel('\DeltaContact % (WT actin - BFP control)')
legend(legend_handles,gene_names)
set(gca,'YLim',[-4,+4],'XLim',[contact_min,contact_max])

title('Over-expression of polymerizing actin')


%% -- more figures

figure(2)
clf

figure(3)
clf

figure(4)
clf

figure(5)
clf

Dist_Delta = zeros(1,5);
Dist_Delta_CI = zeros(2,5);
S5P_Delta = zeros(1,5);
S5P_Delta_CI = zeros(2,5);
S2P_Delta = zeros(1,5);
S2P_Delta_CI = zeros(2,5);

for gg = 1:5

    % Distance-Ser2P scatter plots
    figure(2)
    for tt = 1:numTargets
        subplot(numTargets,numGenes,...
            numGenes.*(target_inds(tt)-1)+gg)

        plot(sorted_paired_distCell{gene_cond_inds{gg}(tt)},...
            sorted_paired_OPIntCell{1}{gene_cond_inds{gg}(tt)},...
            'k.')
        hold on
        plot([1,1].*dist_threshold,[0,3],'k--')
        set(gca,'XLim',[0,8],'YLim',[0,3])
        legend(sprintf('Contact: %2.2f%%',...
            100.*contactFraction_all(gene_cond_inds{gg}(tt))))

        title(sortedCondNames{gene_cond_inds{gg}(tt)})
    end

    % Summary statistics treatment effect on gene-cluster visit
    figure(3)

    subplot(3,5,gg)
%     errorbar(1:4,dist_prctl_all(gene_cond_inds{gg}),...
%         +dist_prctl_CI(2,gene_cond_inds{gg})...
%         -dist_prctl_all(gene_cond_inds{gg})',...
%         -dist_prctl_CI(1,gene_cond_inds{gg})...
%         +dist_prctl_all(gene_cond_inds{gg})',...
%         'ko-','Color',[0.0,0.0,0.0],'LineWidth',1,...
%         'MarkerFaceColor',[0,0,0])
%     ylabel('Gene-cluser dist. [\mum]','FontWeight','normal')
%     set(gca,'XTickLabels',treatment_names,...
%         'XTickLabelRotation',90,'YLim',[0,1.5],'XLim',[0.8,4.2],...
%         'YTick',0:0.5:1.5,'XTick',1:4)
%     title(sortedCondNames(gene_cond_inds{gg}(1)))
%     %title(gene_names{gg})
    errorbar(1:4,contactFraction_all(gene_cond_inds{gg}),...
        +contactFraction_CI(2,gene_cond_inds{gg})...
        -contactFraction_all(gene_cond_inds{gg})',...
        -contactFraction_CI(1,gene_cond_inds{gg})...
        +contactFraction_all(gene_cond_inds{gg})',...
        'ko-','Color',[0.0,0.0,0.0],'LineWidth',1,...
        'MarkerFaceColor',[0,0,0])
    ylabel('Contact fraction','FontWeight','normal')
    set(gca,'XTickLabels',treatment_names,...
        'XTickLabelRotation',90,'YLim',[0,0.2],'XLim',[0.8,4.2],...
        'YTick',0:0.05:0.2,'XTick',1:4)
    title(sortedCondNames(gene_cond_inds{gg}(1)))
    %title(gene_names{gg})

    subplot(3,5,5+gg)
    errorbar(1:4,S5P_mean_all(gene_cond_inds{gg}),...
        +S5P_mean_CI(2,gene_cond_inds{gg})...
        -S5P_mean_all(gene_cond_inds{gg})',...
        -S5P_mean_CI(1,gene_cond_inds{gg})...
        +S5P_mean_all(gene_cond_inds{gg})',...
        'ko-','Color',[1.0,0.0,0.0],'LineWidth',1,...
        'MarkerFaceColor',[1,0,0])
    ylabel('Pol II Ser5P Int. (a.u.)',...
        'FontWeight','normal')
    set(gca,'XTickLabels',treatment_names,...
        'XTickLabelRotation',90,'YLim',[0.95,1.1],'XLim',[0.8,4.2],...
        'YTick',0.9:0.1:1.3,'XTick',1:5)

    subplot(3,5,10+gg)
    errorbar(1:4,S2P_mean_all(gene_cond_inds{gg}),...
        +S2P_mean_CI(2,gene_cond_inds{gg})...
        -S2P_mean_all(gene_cond_inds{gg})',...
        -S2P_mean_CI(1,gene_cond_inds{gg})...
        +S2P_mean_all(gene_cond_inds{gg})',...
        'ko-','Color',[0.5,0.5,0.5],'LineWidth',1,...
        'MarkerFaceColor',[0.5,0.5,0.5])
    ylabel('Pol II Ser2P Int. (a.u.)','FontWeight','normal')
    set(gca,'XTickLabels',treatment_names,...
        'XTickLabelRotation',90,'YLim',[0.95,1.15],'XLim',[0.8,4.2],...
        'YTick',0.8:0.05:1.2,'XTick',1:4)

    Dist_Delta(gg) = ...
        +prctile(...
        sorted_paired_distCell{gene_cond_inds{gg}(3)},prctile_val)...
        -prctile(...
        sorted_paired_distCell{gene_cond_inds{gg}(4)},prctile_val);

    WT_Dist_BS = bootstrp(n_boot,...
        @(xx)prctile(xx,prctile_val),...
        sorted_paired_distCell{gene_cond_inds{gg}(3)});

    R62D_Dist_BS = bootstrp(n_boot,...
        @(xx)prctile(xx,prctile_val),...
        sorted_paired_distCell{gene_cond_inds{gg}(4)});

    Dist_Delta_CI(:,gg) = prctile(WT_Dist_BS-R62D_Dist_BS,...
        [0.025,97.25])';

    S5P_Delta(gg) = ...
        S5P_mean_all(gene_cond_inds{gg}(3))...
        -S5P_mean_all(gene_cond_inds{gg}(4));

    WT_S5P_BS = bootstrp(n_boot,...
        @median,...
        sorted_paired_OPIntCell{2}{gene_cond_inds{gg}(3)});

    R62D_S5P_BS = bootstrp(n_boot,...
        @median,...
        sorted_paired_OPIntCell{2}{gene_cond_inds{gg}(4)});

    S5P_Delta_CI(:,gg) = prctile(WT_S5P_BS-R62D_S5P_BS,...
        [0.025,97.25])';


    S2P_Delta(gg) = ...
        S2P_mean_all(gene_cond_inds{gg}(3))...
        -S2P_mean_all(gene_cond_inds{gg}(2));

    WT_S2P_BS = bootstrp(n_boot,...
        @median,...
        sorted_paired_OPIntCell{1}{gene_cond_inds{gg}(3)});

    R62D_S2P_BS = bootstrp(n_boot,...
        @median,...
        sorted_paired_OPIntCell{1}{gene_cond_inds{gg}(4)});

    S2P_Delta_CI(:,gg) = prctile(WT_S2P_BS-R62D_S2P_BS,...
        [0.025,97.25])';






    figure(4)
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
    hold on
    lineStyles = {'y-','k-','b:','r--'};
    for kk = [1,2,4,3]
        plot(nucVol_support,nucVol_prob{gene_cond_inds{gg}(kk)},...
            lineStyles{kk},'LineWidth',1)
        hold on
    end
    
    xlabel('Nucleus volume V [\mum^3]')
    ylabel('p(V)')
    title(contactFraction_all(gene_cond_inds{gg}(1)))
    title(sortedCondNames{gene_cond_inds{gg}(1)})
            set(gca,'XScale','log')


    subplot(4,5,gg+10)
    hold on
    lineStyles = {'y-','k-','b:','r--'};
    for kk = [1,2,4,3]
        plot(dist_support,NN_dist_prob{gene_cond_inds{gg}(kk)},...
            lineStyles{kk},'LineWidth',1)
        hold on
    end
    
    xlabel('Cluster NN dist. d_{NN} [\mum]')
    ylabel('p(d_{NN})')
    title(contactFraction_all(gene_cond_inds{gg}(1)))
    title(sortedCondNames{gene_cond_inds{gg}(1)})


    subplot(4,5,gg+15)
    hold on
    dotStyles = {'yo','ko','bs','r^'};
    for kk = [1,2,4,3]
        nucMeanNNDist_vals = cellfun(@median, ...
            sorted_S5P_NNDist{gene_cond_inds{gg}(kk)});
        nucVol_vals = [sortedNucVol{gene_cond_inds{gg}(kk)}];
        nucS5P_vals = [sortedNucInt{2}{gene_cond_inds{gg}(kk)}];
        plot(nucVol_vals,nucMeanNNDist_vals,'ko',...
            'MarkerSize',2,'MarkerFaceColor',[0,0,0],...
            'MarkerEdgeColor','none')
    end
    ylabel('Cluster NN distance [\mum]')
    xlabel('Nucleus volume V [\mum^3]')
    legend(treatment_names)
    set(gca,'XLim',[0,1000],'YLim',[0,10])


    figure(5)
    subplot(4,5,gg+15)
    plot(sorted_paired_OPIntCell{1}{gene_cond_inds{gg}(2)},...
        sorted_paired_OPIntCell{2}{gene_cond_inds{gg}(2)},...
        'k.')
    set(gca,'YLim',[0,4],'XLim',[0,4])

    subplot(4,5,gg+10)
    plot(sorted_paired_OPIntCell{1}{gene_cond_inds{gg}(4)},...
        sorted_paired_OPIntCell{2}{gene_cond_inds{gg}(4)},...
        'k.')
    set(gca,'YLim',[0,4],'XLim',[0,4])

    subplot(4,5,gg+15)
    plot(sorted_paired_OPIntCell{1}{gene_cond_inds{gg}(3)},...
        sorted_paired_OPIntCell{2}{gene_cond_inds{gg}(3)},...
        'k.')
    set(gca,'YLim',[0,4],'XLim',[0,4])


end

figure(6)
clf

subplot(3,1,1)
xLocVals = [1:4,5.5];

errorbar(xLocVals,Dist_Delta,...
        +Dist_Delta_CI(1,:)-Dist_Delta,...
        -Dist_Delta_CI(2,:)+Dist_Delta,...
        'ko','Color',[0.0,0.0,0.0],'LineWidth',1,...
        'MarkerFaceColor',[0,0,0])

xlabel('')
ylabel('\Deltad (WT-R62D) [\mum]')
set(gca,'XTick',xLocVals,'XTickLabel',gene_names,...
    'XLim',[0.5,6])
hold on
plot([0.5,6],[0,0],'k-')

subplot(3,1,2)

errorbar(xLocVals,S5P_Delta,...
        +S5P_Delta_CI(1,:)-S5P_Delta,...
        -S5P_Delta_CI(2,:)+S5P_Delta,...
        'ko','Color',[0.0,0.0,0.0],'LineWidth',1,...
        'MarkerFaceColor',[0,0,0])
xlabel('')
ylabel('\DeltaPol II S5P (WT-R62D) [\mum]')
set(gca,'XTick',xLocVals,'XTickLabel',gene_names,...
    'XLim',[0.5,6])
hold on
plot([0.5,6],[0,0],'k-')

subplot(3,1,3)
errorbar(xLocVals,S2P_Delta,...
        +S2P_Delta_CI(1,:)-S2P_Delta,...
        -S2P_Delta_CI(2,:)+S2P_Delta,...
        'ko','Color',[0.0,0.0,0.0],'LineWidth',1,...
        'MarkerFaceColor',[0,0,0])
xlabel('')
ylabel('\DeltaPol II S2P (WT-R62D) [\mum]')
set(gca,'XTick',xLocVals,'XTickLabel',gene_names,...
    'XLim',[0.5,6])
hold on
plot([0.5,6],[0,0],'k-')

%% --- per nucleus overview and selection plots

figure(7)
clf

poolInds = {[1,5,9,13,17],[2,6,10,14,18],...
    [4,8,12,16,20],[3,7,11,15,19]};

for pp = 1:4

    subplot(1,4,pp)
    nucMeanNNDist_vals = cellfun(@median, ...
        [sorted_S5P_NNDist{poolInds{pp}}]);
    nucVol_vals = horzcat(sortedNucVol{poolInds{pp}});
    nucS5P_vals = horzcat(sortedNucInt{2}{poolInds{pp}});
    plot(nucVol_vals,nucMeanNNDist_vals,'ko',...
        'MarkerSize',2,'MarkerFaceColor',[0,0,0],...
        'MarkerEdgeColor','none')
    ylabel('Cluster NN distance [\mum]')
    xlabel('Nucleus volume V [\mum^3]')
    legend(treatment_names)
    set(gca,'XLim',[0,1000],'YLim',[0,10])
    title(treatment_names(target_inds(cc)))

end

%% -- generate output of counts of nuclei and observations

for cc = 1:numConds

    disp(sprintf('%s, %d nuclei, %d genes',...
        sortedCondNames{cc},...
        sortedNumNuclei(cc),numel(sortedDistCell{cc})))

end