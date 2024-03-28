%% --- load analysis results

clear all
load('ConditionSortedResults')

%% -- Interpolation figure

%For interpolation of data
dist_mean_all = [];
S5P_mean_all = [];
S2P_mean_all = [];
%End for interpolation

for cc = 1:numConds

    thisCondInd = cc;

    dist_vals = [sortedDistCell{thisCondInd}];
    OP_S5P_vals = [sortedOPIntCell{2}{thisCondInd}];
    OP_S2P_vals = [sortedOPIntCell{1}{thisCondInd}];
    
    dist_mean_all = [dist_mean_all, mean(dist_vals)];
    S5P_mean_all = [S5P_mean_all,mean(OP_S5P_vals)];
    S2P_mean_all = [S2P_mean_all,mean(OP_S2P_vals)];

end

figure(1)
clf

grid_N = 150; %150
grid_vals = zeros(grid_N,grid_N);
averaging_N = numel(S2P_mean_all); %8
averaging_KK = 0.4;

%Min and max values obtained from figure 2
S2P_min = 0.95;
%S2P_min = floor(min(S2P_mean_all))
S2P_max = 1.3;
S5P_min = 0.9;
S5P_max = 2.0;
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
        grid_vals(mm,nn) = sum(weights.*dist_mean_all)./sum(weights);

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

observMatrix = [S2P_mean_zScore',S5P_mean_zScore'];

plot(S5P_mean_all,S2P_mean_all,...
            'ko','MarkerSize',6,...
            'MarkerEdgeColor',[0,0,0],...
            'MarkerFaceColor',[1,1,1])

% clustRepeats = 50;
% clustMethod = 'gmdistribution';
% 
% for rr = 1:clustRepeats
% 
%     clustRepeats-rr
% 
%     ClustEval = evalclusters(observMatrix,clustMethod,"gap","KList",2:4);
%     numClust = ClustEval.OptimalK;
%     if strcmp(clustMethod,'kmeans')
%         [clustIdx,clusterCentr,sumdist] = ...
%             kmeans(observMatrix,numClust);
%     elseif strcmp(clustMethod,'gmdistribution')
%             mixtModel = fitgmdist(observMatrix,numClust);
%             clustIdx = cluster(mixtModel,observMatrix);
%     end
% 
%     for cc = 1:numClust
% 
%         if sum(cc==clustIdx)>2
%             S2P_clust = S2P_mean_all(cc==clustIdx);
%             S5P_clust = S5P_mean_all(cc==clustIdx);
% 
%             hullInds = convhull(S5P_clust,S2P_clust);
% 
%             plot(S5P_clust(hullInds),S2P_clust(hullInds),'w-',...
%                 'LineWidth',1.5,'Color',[1,1,1,1./clustRepeats])
%         end
% 
%     end
% 
% end

%title('Contact%')
xlabel('Mean S5P Int. (a.u.)')
ylabel('Mean S2P Int. (a.u.)')
c = colorbar;
c.Label.String='Gene-cluster distance [\mum]';
%colormap(flipud(parula))
colormap(parula)
set(gca,'Box','on')