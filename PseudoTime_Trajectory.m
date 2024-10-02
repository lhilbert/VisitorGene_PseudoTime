%% --- load analysis results

clear all
load('ConditionSortedResults')


%% --- Pseudo-time reconstruction

% Choose target gene and stage

target_cond = 2; % 2 - foxd5, 6 - ripply1
all_conds = 1:numel(sortedCondNames);
% all_conds = target_cond;

disp(sortedCondNames{target_cond})

Vol_threshold_pseudotime = 0.0;

dist_vals = [sortedDistCell{target_cond}];
OP_S5P_vals = [sortedOPIntCell{2}{target_cond}];
OP_S2P_vals = [sortedOPIntCell{1}{target_cond}];
OP_OP_vals = [sortedOPIntCell{3}{target_cond}];
Cluster_S5P_vals = [sortedIntCell{2}{target_cond}];
Cluster_S2P_vals = [sortedIntCell{1}{target_cond}];
Cluster_OP_vals = [sortedIntCell{3}{target_cond}];
Vol_vals = [sortedVolCell{target_cond}];
Elo_vals = [sortedEloCell{target_cond}];
Sol_vals = [sortedSolCell{target_cond}];

inclInds = Vol_vals>=Vol_threshold_pseudotime;

dist_vals = dist_vals(inclInds);
OP_S5P_vals = OP_S5P_vals(inclInds);
OP_S2P_vals = OP_S2P_vals(inclInds);
OP_OP_vals = OP_OP_vals(inclInds);
Cluster_S5P_vals = Cluster_S5P_vals(inclInds);
Cluster_S2P_vals = Cluster_S2P_vals(inclInds);
Cluster_OP_vals = Cluster_OP_vals(inclInds);
Vol_vals = Vol_vals(inclInds);
Elo_vals = Elo_vals(inclInds);
Sol_vals = Sol_vals(inclInds);


all_dist_vals = vertcat(sortedDistCell{all_conds});
all_OP_S5P_vals = vertcat(sortedOPIntCell{2}{all_conds});
all_OP_S2P_vals = vertcat(sortedOPIntCell{1}{all_conds});
all_OP_OP_vals = vertcat(sortedOPIntCell{3}{all_conds});
all_Cluster_S5P_vals = vertcat(sortedIntCell{2}{all_conds});
all_Cluster_S2P_vals = vertcat(sortedIntCell{1}{all_conds});
all_Cluster_OP_vals = vertcat(sortedIntCell{3}{all_conds});
all_Vol_vals = vertcat(sortedVolCell{all_conds});
all_Elo_vals = vertcat(sortedEloCell{all_conds});
all_Sol_vals = vertcat(sortedSolCell{all_conds});

inclInds = all_Vol_vals>=Vol_threshold_pseudotime;

all_dist_vals = all_dist_vals(inclInds);
all_OP_S5P_vals = all_OP_S5P_vals(inclInds);
all_OP_S2P_vals = all_OP_S2P_vals(inclInds);
all_OP_OP_vals = all_OP_OP_vals(inclInds);
all_Cluster_S5P_vals = all_Cluster_S5P_vals(inclInds);
all_Cluster_S2P_vals = all_Cluster_S2P_vals(inclInds);
all_Cluster_OP_vals = all_Cluster_OP_vals(inclInds);
all_Vol_vals = all_Vol_vals(inclInds);
all_Elo_vals = all_Elo_vals(inclInds);
all_Sol_vals = all_Sol_vals(inclInds);



% PCA, input: Rows of X are observations, columns the variables

%Create observation matrix with all measured properties
all_observationMatrix = [all_dist_vals,...
    all_OP_S5P_vals,all_OP_S2P_vals,all_OP_OP_vals,...
    all_Cluster_S5P_vals,all_Cluster_S2P_vals,all_Cluster_OP_vals,...
    all_Vol_vals,all_Elo_vals,all_Sol_vals];
target_observationMatrix = [dist_vals,...
    OP_S5P_vals,OP_S2P_vals,OP_OP_vals,...
    Cluster_S5P_vals,Cluster_S2P_vals,Cluster_OP_vals,...
    Vol_vals,Elo_vals,Sol_vals];


%Perform PCA with three components on overall data
[PCA_coeffs,all_PCA_scores,~,~,PCA_percExplained] = ...
    pca(all_observationMatrix,...
    'NumComponents',3);

%Transform target data with same PCA projection
[all_coeff, ~, ~, ~, ~, all_mu] = pca(all_observationMatrix,...
    'NumComponents',3);
% Center the new points
target_Standardized = (target_observationMatrix - all_mu);
% Project the new points onto the PCA axes
target_PCA_scores = target_Standardized * all_coeff;

resortFlag = true;

if resortFlag

    % Get PC that holds most info on Distance changes
    [~,maxIndDist] = max(abs(PCA_coeffs(1,:)));
    % Get PC that holds most info on Cluster Elongation changes
    [~,maxIndElo] = max(abs(PCA_coeffs(9,:)));

    if maxIndElo == maxIndDist
        maxIndElo = mod(maxIndDist+1,2);
    end

    %Set order for PCA (according to max vol cluster)
    PCA_order = [maxIndDist,maxIndElo,...
        setdiff([1,2,3],[maxIndDist,maxIndElo])];
    if PCA_order(1) == PCA_order(2)
        PCA_order = PCA_order(2:end);
    end
else
    PCA_order = 1:3;
end

PCA_coeffs_plot = PCA_coeffs(:,PCA_order);

figure(1);
clf

subplot(2,3,1)
imagesc(PCA_coeffs_plot',[-1,+1])
colorbar

PCA_labels = arrayfun(@(nn) ...
    sprintf('PC %d (%1.1f%%)',...
    PCA_order(nn),PCA_percExplained(PCA_order(nn))),1:3,...
    'UniformOutput',false);

set(gca,'YTick',1:3,'YTickLabel',PCA_labels)
set(gca,'XTick',1:10,'XTickLabel',...
    {'Distance',...
    'OP Ser5P','OP Ser2P','OP OP',...
    'Clus Ser5P','Clus Ser2P','Clus OP',...
    'Clus Vol','Clus Elo','Clus Sol'})

%Order for first two components
PCA_coeffs = PCA_coeffs(:,PCA_order([1,2]));
all_PCA_scores = all_PCA_scores(:,PCA_order([1,2]));
target_PCA_scores = target_PCA_scores(:,PCA_order([1,2]));
PCA_percExplained = PCA_percExplained(PCA_order([1,2]));
PCA_labels = PCA_labels([1,2]);

topCount = 50;
[~,dist_sortInds] = sort(all_dist_vals,'ascend');
dist_topInds = dist_sortInds(1:topCount);
[~,S5P_sortInds] = sort(all_OP_S5P_vals,'ascend');
S5P_topInds = S5P_sortInds(1:topCount);
[~,Vol_sortInds] = sort(all_Vol_vals,'descend');
Vol_topInds = Vol_sortInds(1:topCount);
[~,Elo_sortInds] = sort(all_Elo_vals,'descend');
Elo_topInds = Elo_sortInds(1:topCount);

%Vector for minimal distance
dist_vec = [...
    mean(all_PCA_scores(dist_topInds,1)),...
    mean(all_PCA_scores(dist_topInds,2))];

%Vector for cluster elongation
Elo_vec = [...
    mean(all_PCA_scores(Elo_topInds,1)),...
    mean(all_PCA_scores(Elo_topInds,2))];

%Vector for cluster elongation
S5P_vec = [...
    mean(all_PCA_scores(S5P_topInds,1)),...
    mean(all_PCA_scores(S5P_topInds,2))];


anchor_vec = [-2,-0.2];
%[-1.65,-0.025];[-1.91,-0.14];%[-1.4,-0.245];

subplot(2,3,2)
scatter(all_PCA_scores(:,1),all_PCA_scores(:,2),30,...
    all_dist_vals,'.')
set(gca,'CLim',[0,1.5],'Colormap',flipud(parula))
colorbar
hold on
plot(anchor_vec(1),anchor_vec(2),...
    'ko','MarkerSize',8)
dist_h = plot([anchor_vec(1),dist_vec(1)],...
    [anchor_vec(2),dist_vec(2)],'k-','LineWidth',1);
Elo_h = plot([anchor_vec(1),Elo_vec(1)],...
    [anchor_vec(2),Elo_vec(2)],'r-','LineWidth',1);
% S5P_h = plot([anchor_vec(1),S5P_vec(1)],...
%     [anchor_vec(2),S5P_vec(2)],'b--','LineWidth',1);
xlabel(PCA_labels{1})
ylabel(PCA_labels{2})
%axis equal
set(gca,'Box','on')
%set(gca,'XLim',[-6,3],'YLim',[-2,8])
legend([dist_h,Elo_h],'Max. Cluster Elongation',...
    'Min. Gene-cluster distance')

%% --- Linear transformations of PCA values

dist_vec_shifted = dist_vec - anchor_vec;
Elo_vec_shifted = Elo_vec - anchor_vec;
all_scores_shifted = all_PCA_scores - anchor_vec;
target_scores_shifted = target_PCA_scores - anchor_vec;

if dist_vec_shifted(2)>0
    dist_ortho_vec_shifted = dist_vec_shifted*[0,-1;+1,0];
else
    dist_ortho_vec_shifted = dist_vec_shifted*[0,+1;-1,0];
end

%Perform transformation, such that close contact points East

%Define unit vectors
unit_vec_1 = dist_vec_shifted./norm(dist_vec_shifted);
unit_vec_2 = dist_ortho_vec_shifted./norm(dist_ortho_vec_shifted);

%Define transformation matrix
Trafo_matrix = [unit_vec_1',unit_vec_2'];

%Perform transformation
all_trafo_scores = all_scores_shifted*Trafo_matrix;
target_trafo_scores = target_scores_shifted*Trafo_matrix;
trafo_dist_vec = dist_vec_shifted*Trafo_matrix;
trafo_Elo_vec = Elo_vec_shifted*Trafo_matrix;

subplot(2,3,3)
scatter(all_trafo_scores(:,1),all_trafo_scores(:,2),30,...
    all_dist_vals,'.')
hold on
plot(0,0,'ko','MarkerSize',8)
set(gca,'CLim',[0,1.5],'Colormap',flipud(parula),'Box','on')
colorbar
dist_h = plot([0,trafo_dist_vec(1)],...
    [0,trafo_dist_vec(2)],...
    'k-','LineWidth',1);
Elo_h = plot([0,trafo_Elo_vec(1)],...
    [0,trafo_Elo_vec(2)],...
    'r-','LineWidth',1);
xlabel('Transformed PC_1')
ylabel('Transformed PC_2')
hold off

legend([dist_h,Elo_h],'Max. Cluster Elongation',...
    'Min. Gene-cluster distance')

%% Pseudo-time sorting

figure(1)

subplot(2,3,4)

dd_support = linspace(0,6,250);
dd_edges = dd_support;
dd_centers = (dd_edges(1:end-1)+dd_edges(2:end))./2;

[all_counts,~] = histcounts(all_dist_vals,dd_edges);
[target_counts,~] = histcounts(dist_vals,dd_edges);

all_counts = all_counts./sum(all_counts);
target_counts = target_counts./sum(target_counts);

[all_pp] = ksdensity(all_dist_vals,dd_support,...
    'Support','positive','Bandwidth',0.2);
[target_pp] = ksdensity(dist_vals,dd_support,...
    'Support','positive','Bandwidth',0.2);

% all_h = plot(dd_support,all_pp,'k-','LineWidth',1.5,...
%     'Color',[0.7,0.7,0.1]);
all_h = area(dd_support,all_pp,...
    'FaceColor',[0.7,0.7,0.7],'EdgeColor','none');
hold on
target_h = plot(dd_support,target_pp,'k-','LineWidth',1.5,...
    'Color',[0.3,0.3,0.9]);

prctlDist = prctile(dist_vals,5);
% title(sprintf('%s, 5-percentile dist. %2.2f \\mum',...
%     sortedCondNames{target_cond},prctlDist),...
%     'FontWeight','normal')

prctl_h = plot([1,1].*prctlDist,[0,0.8],'k--','LineWidth',1.5);

legend([target_h,all_h,prctl_h],...
    sortedCondNames{target_cond},...
    'All observations',...
    '5-percentile dist.')

hold off

xlabel('Gene-cluster dist. d [\mum]')
ylabel('P(d)')




subplot(2,3,5)

plot(all_trafo_scores(:,1),all_trafo_scores(:,2),'ko',...
    'MarkerFaceColor',[0.7,0.7,0.7],...
    'MarkerEdgeColor','none',...
    'MarkerSize',2)
hold on
plot(target_trafo_scores(:,1),target_trafo_scores(:,2),'ko',...
    'MarkerFaceColor',[0.3,0.3,0.9],...
    'MarkerEdgeColor','none',...
    'MarkerSize',3)
hold on

dist_h = plot([0,trafo_dist_vec(1)],...
    [0,trafo_dist_vec(2)],...
    'k-','LineWidth',1);
Elo_h = plot([0,trafo_Elo_vec(1)],...
    [0,trafo_Elo_vec(2)],...
    'r-','LineWidth',1);
xlabel('Transformed PC_1')
ylabel('Transformed PC_2')

%set(gca,'XLim',[-6,3],'YLim',[-2,8])

% Obtain pseudo-time sorted data for target condition

%Get angles from x,y coordinate
angles = -atan2(target_trafo_scores(:,2),...
    target_trafo_scores(:,1))./2./pi; % in radians

angles = mod(angles,1);

%Get number of angles
numPoints = numel(angles);

%Define pseudo-time coordinate
coord_s = ((1:numPoints)-1)./numPoints;
register_shift = +0.5;
coord_s = mod(coord_s+register_shift,1);
[angles,sortInds] = sort(angles);

%Calculate moving mean
windowSize = 5;
curveSmooth = @(xx) movmean(...
    padarray(xx,windowSize,'circular','both'),...
    windowSize,'Endpoints','discard');

traj_h = plot(curveSmooth(target_trafo_scores(sortInds,1)),...
    curveSmooth(target_trafo_scores(sortInds,2)),...
    'k-','LineWidth',1.5);

legend([dist_h,Elo_h,traj_h],...
    'Max. Cluster Elongation',...
    'Min. Gene-cluster distance','Reconstructed trajectory')

%Sort all variables
dist_vals = dist_vals(sortInds);
OP_S5P_vals = OP_S5P_vals(sortInds);
OP_S2P_vals = OP_S2P_vals(sortInds);
Cluster_S5P_vals = Cluster_S5P_vals(sortInds);
Cluster_S2P_vals = Cluster_S2P_vals(sortInds);
Vol_vals = Vol_vals(sortInds);
Elo_vals = Elo_vals(sortInds);
Sol_vals = Sol_vals(sortInds);
%coord_s = coord_s(sortInds);

n_boot = 500;
dist_CI = bootci(n_boot,@median,dist_vals);
OP_S5P_CI = bootci(n_boot,@median,OP_S5P_vals);
OP_S2P_CI = bootci(n_boot,@median,OP_S2P_vals);
Cluster_S5P_CI = bootci(n_boot,@median,Cluster_S5P_vals);
Cluster_S2P_CI = bootci(n_boot,@median,Cluster_S2P_vals);
Vol_CI = bootci(n_boot,@median,Vol_vals);
Elo_CI = bootci(n_boot,@median,Elo_vals);
Sol_CI = bootci(n_boot,@median,Sol_vals);

%Bin discretization for plotting of properties vs. pseudo time s

%Define bin borders
windowWidth = 0.07;
numWindows = numPoints;
windowCenters = linspace(0,1,numWindows);
leftEdges = windowCenters-windowWidth./2;
rightEdges = windowCenters+windowWidth./2;

%Empty arrays/cells for each plot
mean_dist = zeros(1,numWindows);
mean_OP_S5P = zeros(1,numWindows);
mean_OP_S2P = zeros(1,numWindows);
mean_Cluster_S5P = zeros(1,numWindows);
mean_Cluster_S2P = zeros(1,numWindows);
mean_Vol = zeros(1,numWindows);
mean_Elo = zeros(1,numWindows);
mean_Sol = zeros(1,numWindows);

%Loop through all windows
for nn = 1:numWindows

    %Get window edges
    thisLeftEdge = leftEdges(nn);
    thisRightEdge = rightEdges(nn);

    %Get window indices
    windowInds = find(coord_s>=thisLeftEdge & coord_s<thisRightEdge);
    if thisRightEdge>1
        windowInds = [windowInds,find(coord_s<(thisRightEdge-1))];
    end
    if thisLeftEdge<0
        windowInds = [windowInds,find(coord_s>(thisLeftEdge+1))];
    end

    %Get mean properties of each window
    mean_dist(nn) = median(dist_vals(windowInds));
    mean_OP_S5P(nn) = median(OP_S5P_vals(windowInds));
    mean_OP_S2P(nn) = median(OP_S2P_vals(windowInds));
    mean_Cluster_S5P(nn) = median(Cluster_S5P_vals(windowInds));
    mean_Cluster_S2P(nn) = median(Cluster_S2P_vals(windowInds));
    mean_Vol(nn) = median(Vol_vals(windowInds));
    mean_Elo(nn) = median(Elo_vals(windowInds));
    mean_Sol(nn) = median(Sol_vals(windowInds));
end
%
% windowCenters = curveSmooth(coord_s');
% mean_dist = curveSmooth(dist_vals);
% mean_OP_S5P = curveSmooth(OP_S5P_vals);
% mean_OP_S2P = curveSmooth(OP_S2P_vals);


%Plot all properties within subplots

figure(4)
clf

square_ind = round(numWindows.*0.465);
square_marker = 'ks';
circle_ind = round(numWindows.*0.5);
circle_marker = 'ko';
diamond_ind = round(numWindows.*0.6);
diamond_marker = 'kd';

subplot(3,3,1)
patch([0,1,1,0],dist_CI([1,1,2,2]),[0.8,0.8,0.8],'EdgeColor','none')
hold on
plot(coord_s,dist_vals,'ko','MarkerFaceColor',[0.4,0.4,0.4],...
    'MarkerEdgeColor','none','MarkerSize',3)
hold on
plot(windowCenters,mean_dist,'k-','LineWidth',1)
legend_array = zeros(1,3);
legend_array(2) = plot(windowCenters(circle_ind),mean_dist(circle_ind),circle_marker,...
    'MarkerFaceColor',[1,0,0]);
legend_array(1) = plot(windowCenters(square_ind),mean_dist(square_ind),square_marker,...
    'MarkerFaceColor',[0,0,1]);
legend_array(3) = plot(windowCenters(diamond_ind),mean_dist(diamond_ind),diamond_marker,...
    'MarkerFaceColor',[0.5,0.5,0.5]);
xlabel('Pseudo-time s')
ylabel('Distance [\mum]')
set(gca,'YLim',[0,4],'Box','on')

legend(legend_array,{'Induced','Associated','Transcribing'})


subplot(3,3,4)
patch([0,1,1,0],OP_S5P_CI([1,1,2,2]),[0.8,0.8,0.8],'EdgeColor','none')
hold on
plot(coord_s,OP_S5P_vals,'ko','MarkerFaceColor',[0.4,0.4,0.4],...
    'MarkerEdgeColor','none','MarkerSize',3)
hold on
plot(windowCenters,mean_OP_S5P,'k-','LineWidth',1)
plot(windowCenters(circle_ind),mean_OP_S5P(circle_ind),circle_marker,...
    'MarkerFaceColor',[1,0,0])
plot(windowCenters(square_ind),mean_OP_S5P(square_ind),square_marker,...
    'MarkerFaceColor',[0,0,1])
plot(windowCenters(diamond_ind),mean_OP_S5P(diamond_ind),diamond_marker,...
    'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Pseudo-time s')
ylabel('Gene Pol II Ser5P')
set(gca,'YLim',[0.5,5],'Box','on')

subplot(3,3,7)
patch([0,1,1,0],OP_S2P_CI([1,1,2,2]),[0.8,0.8,0.8],'EdgeColor','none')
hold on
plot(coord_s,OP_S2P_vals,'ko','MarkerFaceColor',[0.4,0.4,0.4],...
    'MarkerEdgeColor','none','MarkerSize',3)
hold on
plot(windowCenters,mean_OP_S2P,'k-','LineWidth',1)
plot(windowCenters(circle_ind),mean_OP_S2P(circle_ind),circle_marker,...
    'MarkerFaceColor',[1,0,0]);
plot(windowCenters(square_ind),mean_OP_S2P(square_ind),square_marker,...
    'MarkerFaceColor',[0,0,1]);
plot(windowCenters(diamond_ind),mean_OP_S2P(diamond_ind),diamond_marker,...
    'MarkerFaceColor',[0.5,0.5,0.5]);
xlabel('Pseudo-time s')
ylabel('Gene Pol II Ser2P')
set(gca,'YLim',[0.5,2],'Box','on')




subplot(3,3,5)
patch([0,1,1,0],Cluster_S5P_CI([1,1,2,2]),[0.8,0.8,0.8],'EdgeColor','none')
hold on
plot(coord_s,Cluster_S5P_vals,'ko','MarkerFaceColor',[0.4,0.4,0.4],...
    'MarkerEdgeColor','none','MarkerSize',3)
hold on
plot(windowCenters,mean_Cluster_S5P,'k-','LineWidth',1)
plot(windowCenters(circle_ind),mean_Cluster_S5P(circle_ind),circle_marker,...
    'MarkerFaceColor',[1,0,0])
plot(windowCenters(square_ind),mean_Cluster_S5P(square_ind),square_marker,...
    'MarkerFaceColor',[0,0,1])
plot(windowCenters(diamond_ind),mean_Cluster_S5P(diamond_ind),diamond_marker,...
    'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Pseudo-time s')
ylabel('Cluster Pol II Ser5P')
set(gca,'YLim',[0,6],'Box','on')

subplot(3,3,8)
patch([0,1,1,0],Cluster_S2P_CI([1,1,2,2]),[0.8,0.8,0.8],'EdgeColor','none')
hold on
plot(coord_s,Cluster_S2P_vals,'ko','MarkerFaceColor',[0.4,0.4,0.4],...
    'MarkerEdgeColor','none','MarkerSize',3)
hold on
plot(windowCenters,mean_Cluster_S2P,'k-','LineWidth',1)
plot(windowCenters(circle_ind),mean_Cluster_S2P(circle_ind),circle_marker,...
    'MarkerFaceColor',[1,0,0])
plot(windowCenters(square_ind),mean_Cluster_S2P(square_ind),square_marker,...
    'MarkerFaceColor',[0,0,1])
plot(windowCenters(diamond_ind),mean_Cluster_S2P(diamond_ind),diamond_marker,...
    'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Pseudo-time s')
ylabel('Cluster Pol II Ser2P')
set(gca,'YLim',[0,2.5],'Box','on')


subplot(3,3,3)
patch([0,1,1,0],Vol_CI([1,1,2,2]),[0.8,0.8,0.8],'EdgeColor','none')
hold on
plot(coord_s,Vol_vals,'ko','MarkerFaceColor',[0.4,0.4,0.4],...
    'MarkerEdgeColor','none','MarkerSize',3)
hold on
plot(windowCenters,mean_Vol,'k-','LineWidth',1)
legend_array = zeros(1,3);
legend_array(2) = ...
    plot(windowCenters(circle_ind),mean_Vol(circle_ind),circle_marker,...
    'MarkerFaceColor',[1,0,0]);
legend_array(1) = ...
    plot(windowCenters(square_ind),mean_Vol(square_ind),square_marker,...
    'MarkerFaceColor',[0,0,1]);
legend_array(3) = ...
    plot(windowCenters(diamond_ind),mean_Vol(diamond_ind),diamond_marker,...
    'MarkerFaceColor',[0.5,0.5,0.5]);
xlabel('Pseudo-time s')
ylabel('Cluster Volume [\mum^3]')
%set(gca,'YLim',[1,6],'Box','on')

legend(legend_array,{'Induced','Associated','Transcribing'})





subplot(3,3,6)
patch([0,1,1,0],Elo_CI([1,1,2,2]),[0.8,0.8,0.8],'EdgeColor','none')
hold on
plot(coord_s,Elo_vals,'ko','MarkerFaceColor',[0.4,0.4,0.4],...
    'MarkerEdgeColor','none','MarkerSize',3)
hold on
plot(windowCenters,mean_Elo,'k-','LineWidth',1)
legend_array = zeros(1,3);
legend_array(2) = ...
    plot(windowCenters(circle_ind),mean_Elo(circle_ind),circle_marker,...
    'MarkerFaceColor',[1,0,0]);
legend_array(1) = ...
    plot(windowCenters(square_ind),mean_Elo(square_ind),square_marker,...
    'MarkerFaceColor',[0,0,1]);
legend_array(3) = ...
    plot(windowCenters(diamond_ind),mean_Elo(diamond_ind),diamond_marker,...
    'MarkerFaceColor',[0.5,0.5,0.5]);
xlabel('Pseudo-time s')
ylabel('Cluster Elongation')
set(gca,'YLim',[1,6],'Box','on')

%legend(legend_array,{'Induced','Associated','Transcribing'})


subplot(3,3,9)
patch([0,1,1,0],Sol_CI([1,1,2,2]),[0.8,0.8,0.8],'EdgeColor','none')
hold on
plot(coord_s,Sol_vals,'ko','MarkerFaceColor',[0.4,0.4,0.4],...
    'MarkerEdgeColor','none','MarkerSize',3)
hold on
plot(windowCenters,mean_Sol,'k-','LineWidth',1)
legend_array = zeros(1,3);
legend_array(2) = ...
    plot(windowCenters(circle_ind),mean_Sol(circle_ind),circle_marker,...
    'MarkerFaceColor',[1,0,0]);
legend_array(1) = ...
    plot(windowCenters(square_ind),mean_Sol(square_ind),square_marker,...
    'MarkerFaceColor',[0,0,1]);
legend_array(3) = ...
    plot(windowCenters(diamond_ind),mean_Sol(diamond_ind),diamond_marker,...
    'MarkerFaceColor',[0.5,0.5,0.5]);
xlabel('Pseudo-time s')
ylabel('Cluster Solidity')
set(gca,'YLim',[0,1],'Box','on')

%legend(legend_array,{'Induced','Associated','Transcribing'})







subplot(3,3,2)

scatter(OP_S5P_vals,dist_vals,12,OP_S2P_vals,'o','filled')
xlabel('Gene Pol II Ser5P')
ylabel('Distance [\mum]')
colormap(flipud(parula))
colormap(parula)
clim([0.5,1.8])
colorbar
set(gca,'Box','on')









hold on
plot(mean_OP_S5P,mean_dist,'k-','LineWidth',1)
plot(mean_OP_S5P(circle_ind),mean_dist(circle_ind),circle_marker,...
    'MarkerFaceColor',[1,0,0])
plot(mean_OP_S5P(square_ind),mean_dist(square_ind),square_marker,...
    'MarkerFaceColor',[0,0,1])
plot(mean_OP_S5P(diamond_ind),mean_dist(diamond_ind),diamond_marker,...
    'MarkerFaceColor',[0.5,0.5,0.5])



% patch([mean_OP_S5P nan],[mean_OP_S2P nan],[mean_dist nan],...
%     [mean_dist nan], 'edgecolor', 'interp');


%% -- special section on visualization of gene-cluster visit

% To use this section, you first need to use the
% Visualization_Segmentation_FeatureExtraction.m script to extract example
% images for a desired condition and gene. This will be used here to source
% the image data needed for plotting.

if target_cond == 2

    figure(5)
    clf

    imgCell = load(...
        'Visualization_ConditionSortedResults.mat',...
        'sortedOPDisplCell','sortedImgCell',...
        'sortedPixelSize','centralSliceExtension');

    visImgCell = imgCell.sortedImgCell{1}(sortInds);
    visDispls = imgCell.sortedOPDisplCell{1};
    visDispls = visDispls(sortInds,:);
    visPixelSize = imgCell.sortedPixelSize(1);
    extension = imgCell.centralSliceExtension;

    OP_z_Displ = visDispls(:,3);

    img_pixel_count = cellfun(...
        @(xx)prod(size(xx{1})),visImgCell);

    target_ss = [0.03,0.07,0.15,0.41,0.5,...
        0.51,0.545,0.61,0.7,0.77];

    numExamples = numel(target_ss);

    for ee = 1:numExamples

        coord_dist = abs(coord_s - target_ss(ee));
        coord_dist(abs(OP_z_Displ)>0.5) = Inf;
        coord_dist(img_pixel_count<(2.*extension+1).^2) = Inf;

        [~,this_ind] = min(coord_dist);

        OP_img = visImgCell{this_ind}{3};
        S5P_img = visImgCell{this_ind}{2};
        S2P_img = visImgCell{this_ind}{1};

        OP_lims = [1.05,1.4];
        OP_img = OP_img - OP_lims(1);
        OP_img(OP_img<0) = 0;
        OP_img = OP_img./diff(OP_lims);
        OP_img(OP_img>1) = 1;
        OP_img = 1-OP_img;

        S5P_lims = [0.95,2.4];
        S5P_img = S5P_img./median(S5P_img(:));
        S5P_img = S5P_img - S5P_lims(1);
        S5P_img(S5P_img<0) = 0;
        S5P_img = S5P_img./diff(S5P_lims);
        S5P_img(S5P_img>1) = 1;
        S5P_img = 1-S5P_img;

        S2P_lims = [0.7,2.7];
        S2P_img = S2P_img - S2P_lims(1);
        S2P_img(S2P_img<0) = 0;
        S2P_img = S2P_img./diff(S2P_lims);
        S2P_img(S2P_img>1) = 1;

        OP_Displ = visDispls(this_ind,:);

        displAngle = atan2d(OP_Displ(2),OP_Displ(1))...
            +180.*(coord_s(this_ind)<0.5);
        OP_Displ(2) = norm(OP_Displ(1:2));
        OP_Displ(1) = 0;
        if coord_s(this_ind)<0.5
            OP_Displ(2) = OP_Displ(2).*(-1);
        end
%         OP_Displ(1) = OP_Displ(1)./sqrt(2);
%         OP_Displ(2) = OP_Displ(1);
        S2P_img = imrotate(S2P_img,displAngle-90,'crop');
        S5P_img = imrotate(S5P_img,displAngle-90,'crop');
        OP_img = imrotate(OP_img,displAngle-90,'crop');

        %     	subplot(numExamples,4,(ee-1).*4+1)
        %     	imagesc(...
        %             [-extension,extension].*visPixelSize,...
        %             [-extension,extension].*visPixelSize,...
        %             OP_img)
        %         hold on
        %     	axis equal tight
        %     	colormap(gray)
        %     	if ee == 1
        %     		plot(-extension.*visPixelSize+[0.2,1.2],...
        %                 -extension.*visPixelSize+[0.2,0.2],'w-','LineWidth',2)
        %     		title('Oligopaint','FontWeight','normal')
        %         end
        %         plot(OP_Displ(1),OP_Displ(2),'wo',...
        %             'MarkerFaceColor',[1,1,1],...
        %             'MarkerEdgeColor',[0,0,0])
        %     	set(gca,'XLim',[-extension,extension].*visPixelSize,...
        %     		'YLim',[-extension,extension].*visPixelSize,...
        %     		'XTick',[],'YTick',[])
        %
        %     	subplot(numExamples,4,(ee-1).*4+2)
        %     	imagesc(...
        %             [-extension,extension].*visPixelSize,...
        %             [-extension,extension].*visPixelSize,...
        %             S5P_img)
        %     	axis equal tight
        %     	colormap(gray)
        %     	if ee == 1
        %     		hold on
        %     		title('Pol II Ser5P','FontWeight','normal')
        %     	end
        %     	set(gca,'XLim',[-extension,extension].*visPixelSize,...
        %     		'YLim',[-extension,extension].*visPixelSize,...
        %     		'XTick',[],'YTick',[])
        %
        %
        %     	subplot(numExamples,4,(ee-1).*4+3)
        %     	imagesc(...
        %             [-extension,extension].*visPixelSize,...
        %             [-extension,extension].*visPixelSize,...
        %             S2P_img)
        %     	axis equal tight
        %     	colormap(gray)
        %     	if ee == 1
        %     		hold on
        %     		title('Pol II Ser2P','FontWeight','normal')
        %     	end
        %     	set(gca,'XLim',[-extension,extension].*visPixelSize,...
        %     		'YLim',[-extension,extension].*visPixelSize,...
        %     		'XTick',[],'YTick',[])
        %
        %
        %     	subplot(numExamples,4,(ee-1).*4+4)
        subplot(4,numExamples,ee)
        greenChannel = S5P_img;
        greenChannel = greenChannel./max(greenChannel(:));
        magentaChannel = OP_img;
        magentaChannel = magentaChannel./max(magentaChannel(:));
        rgb_img = zeros(size(OP_img,1),size(OP_img,2),3);
        rgb_img(:,:,1) = magentaChannel;
        rgb_img(:,:,3) = greenChannel;
        rgb_img(:,:,2) = min(magentaChannel,greenChannel);
        image(...
            [-extension,extension].*visPixelSize,...
            [-extension,extension].*visPixelSize,...
            rgb_img)
        hold on
        axis equal tight
        if ee == 1
            title('Merge','FontWeight','normal')
            plot(...
                -0.15.*extension.*visPixelSize-[0.4,0.4],...
                -extension.*visPixelSize+[0.6,2.6],...
                'k-','LineWidth',2.5)

        end
        plot([0,OP_Displ(1)],[0,OP_Displ(2)],'k-')
        plot(OP_Displ(1),OP_Displ(2),'wo',...
            'MarkerFaceColor',[1,1,1],...
            'MarkerEdgeColor',[0,0,0],...
            'MarkerSize',4)
        plot(0,0,'ws',...
            'MarkerFaceColor',[1,1,1],...
            'MarkerEdgeColor',[0,0,0],...
            'MarkerSize',4)

        plot([-1,-1,1,1,-1].*0.4,...
            [-1,1,1,-1,-1].*0.4+OP_Displ(2),...
            'k-','LineWidth',1)

%         plot(OP_Displ(1),OP_Displ(2),'wo',...
%             'MarkerFaceColor','none',...
%             'MarkerEdgeColor',[0,0,0],...
%             'MarkerSize',6)
%         plot(0,0,'ws',...
%             'MarkerFaceColor','none',...
%             'MarkerEdgeColor',[0,0,0],...
%             'MarkerSize',6)

        set(gca,'XLim',0.4.*[-extension,extension].*visPixelSize,...
            'YLim',0.95.*[-extension,extension].*visPixelSize,...
            'XTick',[],'YTick',-2:1:+2,'YDir','normal')

        if ee>1
            set(gca,'YTickLabel',[])
        end

        title(sprintf('s=%2.2f, d=%2.2f \\mum',...
            coord_s(this_ind),dist_vals(this_ind)),...
            'FontWeight','normal')

        subplot(4,numExamples,numExamples.*1+ee)
        imagesc(...
            [-extension,extension].*visPixelSize,...
            [-extension,extension].*visPixelSize,...
            1-OP_img,[0,1])

        title(sprintf('d=%2.2f\\mum',dist_vals(this_ind)),...
            'FontWeight','normal')

        colormap(inferno)

        axis equal
        axis tight

        set(gca,'XLim',[-0.4,0.4],...
            'YLim',OP_Displ(2)+[-0.4,0.4],...
            'XTick',[],'YTick',[],'YDir','normal')

        subplot(4,numExamples,numExamples.*2+ee)
        imagesc(...
            [-extension,extension].*visPixelSize,...
            [-extension,extension].*visPixelSize,...
            1-S5P_img,[0,1])

        title(sprintf('I_{Ser5P}=%2.2f',OP_S5P_vals(this_ind)),...
            'FontWeight','normal')

        colormap(inferno)

        axis equal
        axis tight

        set(gca,'XLim',[-0.4,0.4],...
            'YLim',OP_Displ(2)+[-0.4,0.4],...
            'XTick',[],'YTick',[],'YDir','normal')

        subplot(4,numExamples,numExamples.*3+ee)
        
        yellowChannel = S2P_img;
        blueChannel = OP_img;
        rgb_img = zeros(size(OP_img,1),size(OP_img,2),3);
        rgb_img(:,:,1) = yellowChannel;
        rgb_img(:,:,2) = yellowChannel;
        rgb_img(:,:,3) = blueChannel;
        image(...
            [-extension,extension].*visPixelSize,...
            [-extension,extension].*visPixelSize,...
            rgb_img)

        title(sprintf('I_{Ser2P}=%2.2f',OP_S2P_vals(this_ind)),...
            'FontWeight','normal')

        colormap(inferno)

        axis equal
        axis tight

        set(gca,'XLim',[-0.4,0.4],...
            'YLim',OP_Displ(2)+[-0.4,0.4],...
            'XTick',[],'YTick',[],'YDir','normal')


        


    end

end