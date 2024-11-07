%% --- load analysis results

clear all
load('ConditionSortedResults_ActinPerturbation_withQC')


%% Overview plots for all analyzed conditions



dist_prctl_all = zeros(numConds,1);
S5P_mean_all = zeros(numConds,1);
S2P_mean_all = zeros(numConds,1);
S5P_mean_close = zeros(numConds,1);
S5P_mean_far = zeros(numConds,1);
S2P_mean_close = zeros(numConds,1);
S2P_mean_far = zeros(numConds,1);
contactFraction_all = zeros(numConds,1);
nn_observations = zeros(numConds,1);

dist_prctl_CI = zeros(2,numConds);
S5P_mean_CI = zeros(2,numConds);
S5P_mean_close_CI = zeros(2,numConds);
S5P_mean_far_CI = zeros(2,numConds);
S2P_mean_CI = zeros(2,numConds);
S2P_mean_close_CI = zeros(2,numConds);
S2P_mean_far_CI = zeros(2,numConds);
contactFraction_CI = zeros(2,numConds);

prctile_val = 5;
n_boot = 5000;

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
    S5P_mean_close(cc) = mean(OP_S5P_vals(dist_vals<=dist_threshold));
    S5P_mean_far(cc) = mean(OP_S5P_vals(dist_vals>=dist_threshold));
    S2P_mean_close(cc) = mean(OP_S2P_vals(dist_vals<=dist_threshold));
    S2P_mean_far(cc) = mean(OP_S2P_vals(dist_vals>=dist_threshold));

    S5P_mean_CI(:,cc) = bootci(n_boot,...
        @(xx)mean(xx),OP_S5P_vals);
    S5P_mean_close_CI(:,cc) = bootci(n_boot,...
        @(xx)mean(xx),OP_S5P_vals(dist_vals<=dist_threshold));
    S5P_mean_far_CI(:,cc) = bootci(n_boot,...
        @(xx)mean(xx),OP_S5P_vals(dist_vals>dist_threshold));

    S2P_mean_CI(:,cc) = bootci(n_boot,...
        @(xx)mean(xx),OP_S2P_vals);
    S2P_mean_close_CI(:,cc) = bootci(n_boot,...
        @(xx)mean(xx),OP_S5P_vals(dist_vals<=dist_threshold));
    S2P_mean_far_CI(:,cc) = bootci(n_boot,...
        @(xx)mean(xx),OP_S5P_vals(dist_vals>dist_threshold));

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
clim([contact_min,contact_max])




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




% check if contact and transcription changes are temporally segregated

figure(2)
clf

[Rval,Pval] = corrcoef(...
    [contactFraction_all(treatmentCondInds),...
    S5P_mean_close(treatmentCondInds),...
    S5P_mean_far(treatmentCondInds),...
    S2P_mean_close(treatmentCondInds),...
    S2P_mean_far(treatmentCondInds)]);

subplot(2,2,1)

pp = polyfit(contactFraction_all(treatmentCondInds).*100,...
    S5P_mean_close(treatmentCondInds),1);
legend_handle = ...
    plot([contact_min,contact_max],polyval(pp,[contact_min,contact_max]),...
    'r-','LineWidth',1.5);

hold on

errorbar(contactFraction_all(treatmentCondInds).*100,...
    S5P_mean_close(treatmentCondInds),...
    S5P_mean_close(treatmentCondInds)-S5P_mean_close_CI(1,treatmentCondInds)',...
    S5P_mean_close_CI(2,treatmentCondInds)'-S5P_mean_close(treatmentCondInds),...
    'ko')

legend(legend_handle,sprintf('R=%2.2f, P=%4.4f',Rval(1,2),Pval(1,2)),...
    'Location','northwest')
xlabel('Contact %')
ylabel('Pol II Ser5P Int.')
title(sprintf('Gated for d\\leq%2.2f\\mum',dist_threshold),...
    'FontWeight','normal')
set(gca,'XLim',[contact_min,contact_max])

subplot(2,2,2)

pp = polyfit(contactFraction_all(treatmentCondInds).*100,...
    S5P_mean_far(treatmentCondInds),1);
legend_handle = ...
    plot([contact_min,contact_max],polyval(pp,[contact_min,contact_max]),...
    'r-','LineWidth',1.5);

hold on

errorbar(contactFraction_all(treatmentCondInds).*100,...
    S5P_mean_far(treatmentCondInds),...
    S5P_mean_far(treatmentCondInds)-S5P_mean_far_CI(1,treatmentCondInds)',...
    S5P_mean_far_CI(2,treatmentCondInds)'-S5P_mean_far(treatmentCondInds),...
    'ko')

legend(legend_handle,sprintf('R=%2.2f, P=%4.4f',Rval(1,3),Pval(1,3)),...
    'Location','northwest')
xlabel('Contact %')
ylabel('Pol II Ser5P Int.')
title(sprintf('Gated for d>%2.2f\\mum',dist_threshold),...
    'FontWeight','normal')
set(gca,'XLim',[contact_min,contact_max])

subplot(2,2,3)

pp = polyfit(contactFraction_all(treatmentCondInds).*100,...
    S2P_mean_close(treatmentCondInds),1);
legend_handle = ...
    plot([contact_min,contact_max],polyval(pp,[contact_min,contact_max]),...
    'r-','LineWidth',1.5);

hold on

errorbar(contactFraction_all(treatmentCondInds).*100,...
    S2P_mean_close(treatmentCondInds),...
    S2P_mean_close(treatmentCondInds)...
    -S2P_mean_close_CI(1,treatmentCondInds)',...
    S2P_mean_close_CI(2,treatmentCondInds)'...
    -S2P_mean_close(treatmentCondInds),...
    'ko')

legend(legend_handle,sprintf('R=%2.2f, P=%4.4f',Rval(1,4),Pval(1,4)),...
    'Location','northwest')
xlabel('Contact %')
ylabel('Pol II Ser2P Int.')
set(gca,'XLim',[contact_min,contact_max])

subplot(2,2,4)

pp = polyfit(contactFraction_all(treatmentCondInds).*100,...
    S2P_mean_far(treatmentCondInds),1);
legend_handle = ...
    plot([contact_min,contact_max],polyval(pp,[contact_min,contact_max]),...
    'r-','LineWidth',1.5);

hold on

errorbar(contactFraction_all(treatmentCondInds).*100,...
    S2P_mean_far(treatmentCondInds),...
    S2P_mean_far(treatmentCondInds)...
    -S2P_mean_far_CI(1,treatmentCondInds)',...
    S2P_mean_far_CI(2,treatmentCondInds)'...
    -S2P_mean_far(treatmentCondInds),...
    'ko')

legend(legend_handle,sprintf('R=%2.2f, P=%4.4f',Rval(1,5),Pval(1,5)),...
    'Location','northwest')
xlabel('Contact %')
ylabel('Pol II Ser2P Int.')
set(gca,'XLim',[contact_min,contact_max])




%% -- more figures

figure(3)
clf

figure(4)
clf

figure(5)
clf

figure(6)
clf

figure(7)
clf

Dist_Delta = zeros(1,5);
Dist_Delta_CI = zeros(2,5);
S5P_Delta = zeros(1,5);
S5P_Delta_CI = zeros(2,5);
S2P_Delta = zeros(1,5);
S2P_Delta_CI = zeros(2,5);

for gg = 1:5

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

    % Distance-Ser2P scatter plots
    figure(4)
    for tt = 1:numTargets
        subplot(numTargets,numGenes,...
            numGenes.*(target_inds(tt)-1)+gg)

        patch([0,0,dist_threshold,dist_threshold],[0,3,3,0],...
            [0.6,0.6,0.6],'Edgecolor','none')
%         plot([1,1].*dist_threshold,[0,3],...
%             'k-','Color',[0.6,0.6,0.6],'LineWidth',1.5)
        hold on
        points_handle = plot(...
            sorted_paired_distCell{gene_cond_inds{gg}(tt)},...
            sorted_paired_OPIntCell{1}{gene_cond_inds{gg}(tt)},...
            'ko','MarkerSize',1.8,'MarkerFaceColor',[0,0,0]);
        %alpha(points_handle,0.3)
        hold on
        set(gca,'XLim',[0,8],'YLim',[0,3],'Box','on')
        legend(sprintf('Contact: %2.2f%%',...
            100.*contactFraction_all(gene_cond_inds{gg}(tt))))

        xlabel('Gene-cluster distance [\\mum]')
        ylabel('Gene Pol II Ser2P Int.')
        title(sortedCondNames{gene_cond_inds{gg}(tt)})
    end


    


    figure(5)
    subplot(4,5,gg)
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


    subplot(4,5,gg+5)
%     hold on
%     lineStyles = {'y-','k-','b:','r--'};
%     for kk = [1,2,4,3]
%         plot(dist_support,NN_dist_prob{gene_cond_inds{gg}(kk)},...
%             lineStyles{kk},'LineWidth',1)
%         hold on
%     end
%     
%     xlabel('Cluster NN dist. d_{NN} [\mum]')
%     ylabel('p(d_{NN})')
%     title(contactFraction_all(gene_cond_inds{gg}(1)))
%     title(sortedCondNames{gene_cond_inds{gg}(1)})

    hold on
    dotStyles = {'yo','ko','bs','r^'};
    dotColors = {[0,0,0],[0,0,0],[1,0,0],[0,0,1]};
    for kk = [2,4,3]
        OPDist_vals = ...
            [sorted_paired_distCell{gene_cond_inds{gg}(kk)}];
        nucVol_vals = ...
            [sorted_paired_NucVolCell{gene_cond_inds{gg}(kk)}];
        nucS5P_vals = ...
            [sorted_paired_NucIntCell{2}{gene_cond_inds{gg}(kk)}];
        plot(nucS5P_vals,OPDist_vals,'ko',...
            'MarkerSize',2,'MarkerFaceColor',[0,0,0],...
            'MarkerEdgeColor','none')
    end
    
    ylabel('Gene-cluster distance [\mum]')
    xlabel('Nucleus Pol II Ser5P Int.')
    legend(treatment_names)
    set(gca,'XLim',[0,1500],'YLim',[0,10])


    subplot(4,5,gg+10)
    hold on
    dotStyles = {'yo','ko','bs','r^'};
    for kk = [2,4,3]
        OPDist_vals = ...
            [sorted_paired_distCell{gene_cond_inds{gg}(kk)}];
        nucVol_vals = ...
            [sorted_paired_NucVolCell{gene_cond_inds{gg}(kk)}];
        nucS5P_vals = ...
            [sorted_paired_NucIntCell{2}{gene_cond_inds{gg}(kk)}];
        plot(nucVol_vals,OPDist_vals,'ko',...
            'MarkerSize',2,'MarkerFaceColor',[0,0,0],...
            'MarkerEdgeColor','none')
    end
    ylabel('Gene-cluster distance [\mum]')
    xlabel('Nucleus volume V [\mum^3]')
    legend(treatment_names)
    set(gca,'XLim',[0,1000],'YLim',[0,10])


    subplot(4,5,gg+15)
    hold on
    dotStyles = {'yo','ko','bs','r^'};
    for kk = [2,4,3]
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


    
    
    hold on
    dotStyles = {'yo','ko','bs','r^'};
    NucVol_range = [0,1000];
    for kk = 1:3

        figure(6)
        subplot(6,numGenes,2.*numGenes.*(target_inds(kk+1)-2)+gg)
        cla

        plot(NucVol_range,[1,1].*dist_threshold,'k-','Color',[0.6,0.6,0.6])

        hold on

        nucMeanNNDist_vals = cellfun(@median, ...
            sorted_S5P_NNDist{gene_cond_inds{gg}(kk+1)});
        nucVol_vals = [sortedNucVol{gene_cond_inds{gg}(kk+1)}];
        
        plot(nucVol_vals,nucMeanNNDist_vals,'o',...
            'MarkerSize',3,'MarkerFaceColor',[0,0,0],...
            'MarkerEdgeColor','none')

        pp = polyfit(nucVol_vals(isfinite(nucMeanNNDist_vals)),...
            nucMeanNNDist_vals(isfinite(nucMeanNNDist_vals)),1);
        legend_handle = ...
            plot(NucVol_range,...
            polyval(pp,[contact_min,contact_max]),...
            'r-','LineWidth',1.5);

        [Rval,Pval] = corrcoef(...
            [nucVol_vals(isfinite(nucMeanNNDist_vals)),...
            nucMeanNNDist_vals(isfinite(nucMeanNNDist_vals))']);

        if gg == 1 
            ylabel('Cluster NN dist. [\mum]')
        else
            ylabel('')
            set(gca,'YTickLabel',[])
        end
        xlabel('')
        legend(legend_handle,...
            sprintf('R=%2.2f, P=%2.2f',Rval(1,2),Pval(1,2)))
        %legend('Distance threshold','Cluster nearest neighour')
        set(gca,'XLim',NucVol_range,'YLim',[0,10],'Box','on',...
            'XTickLabels',[],'YTick',0:2:10,'XTick',0:200:1000)

        title(sortedCondNames{gene_cond_inds{gg}(kk+1)})


        subplot(6,numGenes,2.*numGenes.*(target_inds(kk+1)-1.5)+gg)
        cla


        plot(NucVol_range,[1,1].*dist_threshold,'k-','Color',[0.6,0.6,0.6])

        hold on

        OPDist_vals = ...
            [sorted_paired_distCell{gene_cond_inds{gg}(kk+1)}];
        OPnucVol_vals = ...
            [sorted_paired_NucVolCell{gene_cond_inds{gg}(kk+1)}];
        
        plot(OPnucVol_vals,OPDist_vals,'o',...
            'MarkerSize',3,'MarkerFaceColor',[1,0,0],...
            'MarkerEdgeColor','none')

        pp = polyfit(OPnucVol_vals,OPDist_vals,1);
        legend_handle = ...
            plot(NucVol_range,...
            polyval(pp,[contact_min,contact_max]),...
            'k-','LineWidth',1.5);

        [Rval,Pval] = corrcoef(...
            [OPnucVol_vals,OPDist_vals]);

        if gg == 1
            ylabel('Gene-cluster dist. [\mum]')
        else
            ylabel('')
            set(gca,'YTickLabel',[])
        end

        xlabel('Nucleus volume V [\mum^3]')
        legend(legend_handle,...
            sprintf('R=%2.2f, P=%2.2f',Rval(1,2),Pval(1,2)))
        %legend('Distance threshold','Gene-cluster')
        set(gca,'XLim',NucVol_range,'YLim',[0,10],'Box','on',...
            'YTick',0:2:10,'XTick',0:200:NucVol_range(2))

    end

end

%% ---

figure(7)
clf

dotColors = {[0,0,0],[0,0,0],[1,0,0],[0,0,1]};

NucVol_range = [0,700];

for kk = 1:3

    allVols = [];

    Vols = [];
    Dists = [];

    OPVols = [];
    OPDists = [];

    for gg = 1:numGenes

        nucMeanNNDist_vals = cellfun(@mean, ...
            sorted_S5P_NNDist{gene_cond_inds{gg}(kk+1)});
        nucVol_vals = [sortedNucVol{gene_cond_inds{gg}(kk+1)}];

        allVols = [allVols;nucVol_vals];

        Vols = [Vols;nucVol_vals(isfinite(nucMeanNNDist_vals))];
        Dists = [Dists,nucMeanNNDist_vals(isfinite(nucMeanNNDist_vals))];

        OPVols = [OPVols;...
            sorted_paired_NucVolCell{gene_cond_inds{gg}(kk+1)}];
        OPDists = [OPDists; ...
            sorted_paired_distCell{gene_cond_inds{gg}(kk+1)}];

    end


    subplot(4,1,1)

    [prob_vec,support_vec] = ksdensity(allVols,...
        'Support','positive');
    plot(support_vec,prob_vec,'k-','LineWidth',1.5,...
        'Color',dotColors{kk+1})

    xlabel('Nuclear volume (\mum^3)')
    ylabel('Probability density')

    set(gca,'XLim',NucVol_range)
    title('All detected nuclei',...
        'FontWeight','normal')


    hold on



    subplot(4,1,2)

    [sortedVols,sortInds] = ...
        sort(OPVols(OPDists<=dist_threshold),'ascend');

    stairs([0;sortedVols],[0:numel(sortedVols)]./numel(sortedVols),...
        'Color',dotColors{kk+1},'LineWidth',1.5)

    hold on

    prctlVols = prctile(sortedVols,[2.5,97.5]);

    for ll = 1:2
        handle_val = plot([1,1].*prctlVols(ll),[0,1],...
            'k:','LineWidth',1.5,'Color',dotColors{kk+1});
    end
    
    if kk == 1
        prctl_handle = handle_val;
    end

    xlabel('Nuclear volume (\mum^3)')
    ylabel('Cumulative probability')
    set(gca,'XLim',NucVol_range)
    title(sprintf('Gene detections with cluster distance d\\leq%2.2f \\mum',...
        dist_threshold),...
        'FontWeight','normal')

    subplot(4,1,3)

    [sortedVols,sortInds] = sort(OPVols,'ascend');

    stairs([0;sortedVols],[0:numel(sortedVols)]./numel(sortedVols),...
        'Color',dotColors{kk+1},'LineWidth',1.5)
    xlabel('Nuclear volume (\mum^3)')
    ylabel('Cumulative probability')
    set(gca,'XLim',NucVol_range)
    title('All gene detections',...
        'FontWeight','normal')


    hold on

    for ll = 1:2
        plot([1,1].*prctlVols(ll),[0,1],...
            'k:','LineWidth',1.5,'Color',dotColors{kk+1})
    end
    

    subplot(4,1,4)

    window_kk = 64;

    [sortedVols,sortInds] = sort(Vols,'ascend');
    sortedDists = Dists(sortInds);
    plot(sortedVols,...
        movmad(sortedDists,window_kk),...
        'k-','Linewidth',1.5,'Color',dotColors{kk+1},...
        'MarkerEdgeColor','none')
    hold on

    for ll = 1:2
        plot([1,1].*prctlVols(ll),[0,1.2],...
            'k:','LineWidth',1.5,'Color',dotColors{kk+1})
    end

    set(gca,'XLim',NucVol_range,'YLim',[0,1.2])
    ylabel('NN Distance MAD (\mum)')
    xlabel('Nuclear volume (\mum^3)')
    title('Pol II Ser5P clusters, mean nearest-neighbor distance',...
        'FontWeight','normal')

end

subplot(4,1,1)

legend(treatment_names([2,4,3]),'Location','northeast')

subplot(4,1,2)

legend(prctl_handle,'2.5 to 97.5-percentile range','Location','northwest')


%% --- whole-nucleus level analyses

figure(8)
clf

minNucVol = 0; % Miminal nuclear volume
maxNucVol = Inf; % Maxinal nuclear volume
n_boot = 1000;

mapLimits = {[0,6,0,3],[0,8,0,6]};

S5P_threshold = 0.0;

plotStyles = {'k-','r--'};

numTargets = numel(unique(target_inds));

% Collect for distribution plots and normalization

Nuc_group_vec = [];
Nuc_Vol_vec = [];
Nuc_S5P_vec = [];
Nuc_S2P_vec = [];
Nuc_ClusterNum_vec_unscaled = [];
Nuc_ClusterNum_vec = [];

Num_S5P_vec = [];
Num_S2P_vec = [];

S5P_Vol_vec = [];
S5P_group_vec = [];

S5P_NN_vec = [];
NN_group_vec = [];

for kk = 1:4
        
    % Collect all data needed to analyze this data set
    
    for gg = 1:numGenes

        Nuc_Vol_vals = [sortedNucVol{gene_cond_inds{gg}(kk)}];
        Nuc_S5P_vals = [sortedNucInt{1}{gene_cond_inds{gg}(kk)}]...
            -[sortedCytoInt{1}{gene_cond_inds{gg}(kk)}];
        Nuc_S2P_vals = [sortedNucInt{2}{gene_cond_inds{gg}(kk)}]...
            -[sortedCytoInt{2}{gene_cond_inds{gg}(kk)}];
        
        tt = target_inds(kk);

        Nuc_group_vec = vertcat(Nuc_group_vec,tt.*ones(size(Nuc_S5P_vals')));
        Nuc_Vol_vec = vertcat(Nuc_Vol_vec,Nuc_Vol_vals);
        Nuc_S5P_vec = vertcat(Nuc_S5P_vec,Nuc_S5P_vals');
        Nuc_S2P_vec = vertcat(Nuc_S2P_vec,Nuc_S2P_vals');

        S5P_Vol_vals = sortedS5PVolCell{gene_cond_inds{gg}(kk)};
        S5P_Vol_vec = vertcat(S5P_Vol_vec,S5P_Vol_vals);
        S5P_group_vec = [S5P_group_vec,...
            tt.*ones(size(S5P_Vol_vals'))];

        S5P_NN_vals = cellfun(@mean,...
            sorted_S5P_NNDist{gene_cond_inds{gg}(kk)});
        S5P_NN_vals = S5P_NN_vals(S5P_NN_vals>0);
        S5P_NN_vec = [S5P_NN_vec,S5P_NN_vals];
        NN_group_vec = [NN_group_vec,tt.*ones(size(S5P_NN_vals))];

    end

end


Nuc_Vol_mean = zeros(1,numTargets);
Nuc_Vol_CI = zeros(2,numTargets);
Nuc_S5P_mean = zeros(1,numTargets);
Nuc_S5P_CI = zeros(2,numTargets);
Nuc_S2P_mean = zeros(1,numTargets);
Nuc_S2P_CI = zeros(2,numTargets);

S5P_Vol_mean = zeros(1,numTargets);
S5P_Vol_CI = zeros(2,numTargets);

S5P_NN_mean = zeros(1,numTargets);
S5P_NN_CI = zeros(2,numTargets);

for tt = 1:numTargets

    targetInclInds = Nuc_group_vec==tt ...
        & Nuc_Vol_vec>=minNucVol & Nuc_Vol_vec<=maxNucVol;

    Nuc_S5P_mean(tt) = mean(Nuc_S5P_vec(targetInclInds));
    Nuc_S5P_CI(:,tt) = ...
        bootci(n_boot,@mean,Nuc_S5P_vec(targetInclInds));

    Nuc_S2P_mean(tt) = mean(Nuc_S2P_vec(targetInclInds));
    Nuc_S2P_CI(:,tt) = ...
        bootci(n_boot,@mean,Nuc_S2P_vec(targetInclInds));

    targetInclInds = Nuc_group_vec==tt;

    Nuc_Vol_mean(tt) = mean(Nuc_Vol_vec(targetInclInds));
    Nuc_Vol_CI(:,tt) = ...
        bootci(n_boot,@mean,Nuc_Vol_vec(targetInclInds));

    targetInclInds = S5P_group_vec==tt;

    S5P_Vol_mean(tt) = mean(S5P_Vol_vec(targetInclInds));
    S5P_Vol_CI(:,tt) = ...
        bootci(n_boot,@mean,S5P_Vol_vec(targetInclInds));

    targetInclInds = NN_group_vec==tt;

    S5P_NN_mean(tt) = mean(S5P_NN_vec(targetInclInds));
    S5P_NN_CI(:,tt) = ...
        bootci(n_boot,@mean,S5P_NN_vec(targetInclInds));


end

Nuc_S5P_mean = Nuc_S5P_mean./mean(Nuc_S5P_vec);
Nuc_S5P_CI = Nuc_S5P_CI./mean(Nuc_S5P_vec);
Nuc_S2P_mean = Nuc_S2P_mean./mean(Nuc_S2P_vec);
Nuc_S2P_CI = Nuc_S2P_CI./mean(Nuc_S2P_vec);

Nuc_S5P_vec = Nuc_S5P_vec./mean(Nuc_S5P_vec);
Nuc_S2P_vec = Nuc_S2P_vec./mean(Nuc_S2P_vec);
 
    
%% --- Distribution plots



subplot(2,5,1)

distributionPlot(Nuc_Vol_vec,'groups',Nuc_group_vec,...
    'xNames',treatment_names(target_inds(1:4)),'color',[0.6,0.6,0.6],...
    'showMM',0,'addSpread',0)
hold on

errorbar(1:numTargets,Nuc_Vol_mean,...
    Nuc_Vol_CI(1,:)-Nuc_Vol_mean,...
    Nuc_Vol_mean-Nuc_Vol_CI(2,:),...
    'k-o','MarkerFaceColor',[1,1,1],...
    'LineWidth',1);

set(gca,'XLim',[0.5,numTargets+0.5],'Box','on','YLim',[-Inf,+Inf],...
    'XTick',1:numTargets,'XTickLabel',treatment_names(target_inds(1:4)),...
    'XLim',[0.3,numTargets+0.8],'XTickLabelRotation',-45)

ylabel('Nucleus volume (\mum^3)')

subplot(2,5,6)

errorbar(1:numTargets,Nuc_Vol_mean,...
    Nuc_Vol_CI(1,:)-Nuc_Vol_mean,...
    Nuc_Vol_mean-Nuc_Vol_CI(2,:),...
    'k-o','MarkerFaceColor',[1,1,1],...
    'LineWidth',1);

set(gca,'XLim',[0.5,numTargets+0.5],'Box','on','YLim',[-Inf,+Inf],...
    'XTick',1:numTargets,'XTickLabel',treatment_names(target_inds(1:4)),...
    'XLim',[0.3,numTargets+0.8],'XTickLabelRotation',-45)

ylabel('Nucleus volume (\mum^3)')





subplot(2,5,2)

distributionPlot(Nuc_S5P_vec,'groups',Nuc_group_vec,...
    'xNames',treatment_names(target_inds(1:4)),'color',[0.6,0.6,0.6],...
    'showMM',0,'addSpread',0)
hold on

errorbar(1:numTargets,Nuc_S5P_mean,...
    Nuc_S5P_CI(1,:)-Nuc_S5P_mean,...
    Nuc_S5P_mean-Nuc_S5P_CI(2,:),...
    'k-o','MarkerFaceColor',[1,1,1],...
    'LineWidth',1);

set(gca,'XLim',[0.5,numTargets+0.5],'Box','on','YLim',[-Inf,+Inf],...
    'XTick',1:numTargets,'XTickLabel',treatment_names(target_inds(1:4)),...
    'XLim',[0.3,numTargets+0.8],'XTickLabelRotation',-45)

ylabel('Nuclear S5P level')


subplot(2,5,7)

errorbar(1:numTargets,Nuc_S5P_mean,...
    Nuc_S5P_CI(1,:)-Nuc_S5P_mean,...
    Nuc_S5P_mean-Nuc_S5P_CI(2,:),...
    'k-o','MarkerFaceColor',[1,1,1],...
    'LineWidth',1);

set(gca,'XLim',[0.5,numTargets+0.5],'Box','on','YLim',[-Inf,+Inf],...
    'XTick',1:numTargets,'XTickLabel',treatment_names(target_inds(1:4)),...
    'XLim',[0.3,numTargets+0.8],'XTickLabelRotation',-45)

ylabel('Nuclear S5P level')






subplot(2,5,3)

distributionPlot(Nuc_S2P_vec,'groups',Nuc_group_vec,...
    'xNames',treatment_names(target_inds(1:4)),'color',[0.6,0.6,0.6],...
    'showMM',0,'addSpread',0)
hold on

errorbar(1:numTargets,Nuc_S2P_mean,...
    Nuc_S2P_CI(1,:)-Nuc_S2P_mean,...
    Nuc_S2P_mean-Nuc_S2P_CI(2,:),...
    'k-o','MarkerFaceColor',[1,1,1],...
    'LineWidth',1);

set(gca,'XLim',[0.5,numTargets+0.5],'Box','on','YLim',[-Inf,+Inf],...
    'XTick',1:numTargets,'XTickLabel',treatment_names(target_inds(1:4)),...
    'XLim',[0.3,numTargets+0.8],'XTickLabelRotation',-45)

ylabel('Nuclear S2P level')

subplot(2,5,8)

errorbar(1:numTargets,Nuc_S2P_mean,...
    Nuc_S2P_CI(1,:)-Nuc_S2P_mean,...
    Nuc_S2P_mean-Nuc_S2P_CI(2,:),...
    'k-o','MarkerFaceColor',[1,1,1],...
    'LineWidth',1);

set(gca,'XLim',[0.5,numTargets+0.5],'Box','on','YLim',[-Inf,+Inf],...
    'XTick',1:numTargets,'XTickLabel',treatment_names(target_inds(1:4)),...
    'XLim',[0.3,numTargets+0.8],'XTickLabelRotation',-45)

ylabel('Nuclear S2P level')



subplot(2,5,4)

distributionPlot(S5P_Vol_vec,'groups',S5P_group_vec,...
    'xNames',treatment_names(target_inds(1:4)),'color',[0.6,0.6,0.6],...
    'showMM',0,'addSpread',0)
hold on

errorbar(1:numTargets,S5P_Vol_mean,...
    S5P_Vol_CI(1,:)-S5P_Vol_mean,...
    S5P_Vol_mean-S5P_Vol_CI(2,:),...
    'k-o','MarkerFaceColor',[1,1,1],...
    'LineWidth',1);

set(gca,'XLim',[0.5,numTargets+0.5],'Box','on','YLim',[-Inf,+Inf],...
    'XTick',1:numTargets,'XTickLabel',treatment_names(target_inds(1:4)),...
    'XLim',[0.3,numTargets+0.8],'XTickLabelRotation',-45)

ylabel('Cluster volume (\mum^3)')

subplot(2,5,9)

errorbar(1:numTargets,S5P_Vol_mean,...
    S5P_Vol_CI(1,:)-S5P_Vol_mean,...
    S5P_Vol_mean-S5P_Vol_CI(2,:),...
    'k-o','MarkerFaceColor',[1,1,1],...
    'LineWidth',1);

set(gca,'XLim',[0.5,numTargets+0.5],'Box','on','YLim',[-Inf,+Inf],...
    'XTick',1:numTargets,'XTickLabel',treatment_names(target_inds(1:4)),...
    'XLim',[0.3,numTargets+0.8],'XTickLabelRotation',-45)

ylabel('Cluster volume (\mum^3)')



subplot(2,5,5)

distributionPlot(S5P_NN_vec','groups',NN_group_vec',...
    'xNames',treatment_names(target_inds(1:4)),'color',[0.6,0.6,0.6],...
    'showMM',0,'addSpread',0)
hold on

errorbar(1:numTargets,S5P_NN_mean,...
    S5P_NN_CI(1,:)-S5P_NN_mean,...
    S5P_NN_mean-S5P_NN_CI(2,:),...
    'k-o','MarkerFaceColor',[1,1,1],...
    'LineWidth',1);

set(gca,'XLim',[0.5,numTargets+0.5],'Box','on','YLim',[-Inf,+Inf],...
    'XTick',1:numTargets,'XTickLabel',treatment_names(target_inds(1:4)),...
    'XLim',[0.3,numTargets+0.8],'XTickLabelRotation',-45)

ylabel('Cluster NN distance (\mum)')

subplot(2,5,10)

errorbar(1:numTargets,S5P_NN_mean,...
    S5P_NN_CI(1,:)-S5P_NN_mean,...
    S5P_NN_mean-S5P_NN_CI(2,:),...
    'k-o','MarkerFaceColor',[1,1,1],...
    'LineWidth',1);

set(gca,'XLim',[0.5,numTargets+0.5],'Box','on','YLim',[-Inf,+Inf],...
    'XTick',1:numTargets,'XTickLabel',treatment_names(target_inds(1:4)),...
    'XLim',[0.3,numTargets+0.8],'XTickLabelRotation',-45)

ylabel('Cluster NN distance (\mum)')



%% -- generate output of counts of nuclei and observations

for cc = 1:numConds

    disp(sprintf('%s, %d nuclei, %d gene-cluster pairs',...
        sortedCondNames{cc},...
        sortedNumNuclei(cc),numel(sorted_paired_distCell{cc})))

end