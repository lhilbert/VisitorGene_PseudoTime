%% --- load analysis results

clear all
load('ConditionSortedResults')

dist_prctl_all = zeros(numConds,1);
S5P_mean_all = zeros(numConds,1);
S2P_mean_all = zeros(numConds,1);

dist_prctl_CI = zeros(2,numConds);
S5P_mean_CI = zeros(2,numConds);
S2P_mean_CI = zeros(2,numConds);

n_boot = 500;

for cc = 1:numConds

    thisCondInd = cc;

    dist_vals = [sortedDistCell{thisCondInd}];
    OP_S5P_vals = [sortedOPIntCell{2}{thisCondInd}];
    OP_S2P_vals = [sortedOPIntCell{1}{thisCondInd}];
    
    dist_prctl_all(cc) =  prctile(dist_vals,10);
    S5P_mean_all(cc) = mean(OP_S5P_vals);
    S2P_mean_all(cc) = mean(OP_S2P_vals);

    S5P_mean_CI(:,cc) = bootci(n_boot,...
        @(xx)mean(xx),OP_S5P_vals);

    S2P_mean_CI(:,cc) = bootci(n_boot,...
        @(xx)mean(xx),OP_S2P_vals);

    dist_prctl_CI(:,cc) = bootci(n_boot,...
        @(xx)prctile(xx,10),dist_vals);

end

figure(1)
clf

S2P_min = 0.95;
S2P_max = 1.3;
S5P_min = 0.9;
S5P_max = 1.5;


%% -- Overview figure



gene_names = {...
    'drll2','foxd5','gadd45ga',...
    'iscub','klf2b','ripply1',...
    'vamp2','zgc::64022'};

stage_names = {...
    'Oblong','Sphere','Dome',...
    '30% Epi.','50% Epi.'};

gene_trajectory_inds = {...
    [1,9,17,25,33]+0,...
    [1,9,17,25,33]+1,...
    [1,9,17,25,33]+2,...
    [1,9,17,25,33]+3,...
    [1,9,17,25,33]+4,...
    [1,9,17,25,33]+5,...
    [1,9,17,25,33]+6,...
    [1,9,17,25,33]+7};

example_gg = 5;

subplot(3,3,1)
errorbar(1:5,S5P_mean_all(gene_trajectory_inds{example_gg}),...
    +S5P_mean_CI(2,gene_trajectory_inds{example_gg})...
    -S5P_mean_all(gene_trajectory_inds{example_gg})',...
    -S5P_mean_CI(1,gene_trajectory_inds{example_gg})...
    +S5P_mean_all(gene_trajectory_inds{example_gg})',...
    'ko-','Color',[1.0,0.0,0.0],'LineWidth',1,...
    'MarkerFaceColor',[1,0,0])
title('Pol II Ser5P Int. (a.u.)',...
    'FontWeight','normal')
set(gca,'XTickLabels','','YLim',[0.9,1.5],'XLim',[0.8,5.2],...
    'YTick',0.9:0.1:1.5)

ylabel(gene_names(example_gg))


subplot(3,3,4)
errorbar(1:5,S2P_mean_all(gene_trajectory_inds{example_gg}),...
    +S2P_mean_CI(2,gene_trajectory_inds{example_gg})...
    -S2P_mean_all(gene_trajectory_inds{example_gg})',...
    -S2P_mean_CI(1,gene_trajectory_inds{example_gg})...
    +S2P_mean_all(gene_trajectory_inds{example_gg})',...
    'ko-','Color',[0.5,0.5,0.5],'LineWidth',1,...
    'MarkerFaceColor',[0.5,0.5,0.5])
title('Pol II Ser2P Int. (a.u.)','FontWeight','normal')
set(gca,'XTickLabels','','YLim',[0.95,1.15],'XLim',[0.8,5.2],...
    'YTick',0.95:0.05:1.15)


subplot(3,3,7)
errorbar(1:5,dist_prctl_all(gene_trajectory_inds{example_gg}),...
    +dist_prctl_CI(2,gene_trajectory_inds{example_gg})...
    -dist_prctl_all(gene_trajectory_inds{example_gg})',...
    -dist_prctl_CI(1,gene_trajectory_inds{example_gg})...
    +dist_prctl_all(gene_trajectory_inds{example_gg})',...
    'ko-','Color',[0.0,0.0,0.0],'LineWidth',1,...
    'MarkerFaceColor',[0,0,0])
title('Gene-cluser dist. [\mum]','FontWeight','normal')
set(gca,'XTickLabels',stage_names,...
    'XTickLabelRotation',-0,'YLim',[0,1.5],'XLim',[0.8,5.2],...
    'YTick',0:0.5:1.5)



subplot(1,3,2)
cla

gene_handles = zeros(1,8);

for gg = 1:8

    inds = gene_trajectory_inds{gg};
    gene_handles(gg) = plot(...
        S5P_mean_all(inds),...
        S2P_mean_all(inds),...
        'LineStyle','-','LineWidth',1.0);
    hold on
    plot(...
        S5P_mean_all(inds(1)),...
        S2P_mean_all(inds(1)),...
        'ks','MarkerSize',10)

end

scatter(S5P_mean_all,S2P_mean_all,40,...
   dist_prctl_all,'filled')

xlabel('Mean S5P Int. (a.u.)')
ylabel('Mean S2P Int. (a.u.)')

set(gca,'XLim',[1,1.4],...
    'YLim',[S2P_min,S2P_max])

legend(gene_handles,gene_names)

colorbar

%% -- Interpolation figure

subplot(1,3,3)

grid_N = 150; %150
grid_vals = zeros(grid_N,grid_N);
averaging_N = numel(S2P_mean_all); %8
averaging_KK = 0.2;

S2P_vec = linspace(S2P_min,S2P_max,grid_N);
S5P_vec = linspace(S5P_min,S5P_max,grid_N);

for mm = 1:grid_N
    for nn = 1:grid_N

        S2P_val = S2P_vec(mm);
        S5P_val = S5P_vec(nn);

        grid_dist_vec = ...
            + (S2P_mean_all-S2P_val).^2./var(S2P_mean_all) ...
            + (S5P_mean_all-S5P_val).^2./var(S5P_mean_all);

        %[~,sort_inds] = sort(grid_dist_vec,'ascend');
        %grid_vals(mm,nn) = mean(Contact_mean_all(sort_inds(1:averaging_N)));

        weights = averaging_KK./(averaging_KK+grid_dist_vec);
        grid_vals(mm,nn) = sum(weights.*dist_prctl_all)./sum(weights);

    end
end

cla
imagesc(S5P_vec,S2P_vec,grid_vals)
set(gca,'YDir','normal')
colorbar
hold on


% --- clustering

S2P_mean_zScore = (S2P_mean_all-mean(S2P_mean_all))...
    ./std(S2P_mean_all);
S5P_mean_zScore = (S5P_mean_all-mean(S5P_mean_all))...
    ./std(S5P_mean_all);

observMatrix = [S2P_mean_zScore,S5P_mean_zScore];

plot(S5P_mean_all,S2P_mean_all,...
            'ko','MarkerSize',6,...
            'MarkerEdgeColor',[0,0,0],...
            'MarkerFaceColor',[1,1,1])

clustRepeats = 50;
clustMethod = 'kmeans';

numClust_vec = zeros(1,clustRepeats);

for rr = 1:clustRepeats

    clustRepeats-rr

    ClustEval = evalclusters(observMatrix,clustMethod,"gap","KList",2:4);
    numClust = ClustEval.OptimalK;
    if strcmp(clustMethod,'kmeans')
        [clustIdx,clusterCentr,sumdist] = ...
            kmeans(observMatrix,numClust,'Replicates',10);
    elseif strcmp(clustMethod,'gmdistribution')
            mixtModel = fitgmdist(observMatrix,numClust,'Replicates',10);
            clustIdx = cluster(mixtModel,observMatrix);
    end

    numClust_vec(rr) = numClust;

    for cc = 1:numClust

        if sum(cc==clustIdx)>2
            S2P_clust = S2P_mean_all(cc==clustIdx);
            S5P_clust = S5P_mean_all(cc==clustIdx);

            hullInds = convhull(S5P_clust,S2P_clust);

            plot(S5P_clust(hullInds),S2P_clust(hullInds),'k-',...
                'LineWidth',1.5,'Color',[0,0,0,1./clustRepeats])
        end

    end

end

%title('Contact%')
xlabel('Mean S5P Int. (a.u.)')
ylabel('Mean S2P Int. (a.u.)')
c = colorbar;
c.Label.String='10-percentile distance [\mum]';
%colormap(flipud(parula))
colormap(flipud(parula))
set(gca,'Box','on')