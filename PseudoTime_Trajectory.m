%% --- load analysis results

clear all
load('ConditionSortedResults')


%% --- Pseudo-time reconstruction

% Choose target gene and stage

cc = 1;

disp(sortedCondNames{cc})
geneName = 'zgc';
Vol_threshold_pseudotime = 0.0;0.03;
inclInds = sortedVolCell{cc}>=Vol_threshold_pseudotime;

dist_vals = [sortedDistCell{cc}(inclInds)];
OP_S5P_vals = [sortedOPIntCell{2}{cc}(inclInds)];
OP_S2P_vals = [sortedOPIntCell{1}{cc}(inclInds)];
%Add for third channel
Cluster_S5P_vals = [sortedIntCell{2}{cc}(inclInds)];
Cluster_S2P_vals = [sortedIntCell{1}{cc}(inclInds)];
Vol_vals = [sortedVolCell{cc}(inclInds)];
Elo_vals = [sortedEloCell{cc}(inclInds)];
Sol_vals = [sortedSolCell{cc}(inclInds)];
Central_slices = [sortedMaskCell{cc}(inclInds)];

% PCA, input: Rows of X are observations, columns the variables

%Create observation matrix with all measured properties
observationMatrix = [dist_vals,...
    OP_S5P_vals,OP_S2P_vals, ...
    Cluster_S5P_vals,Cluster_S2P_vals,...
    Vol_vals,Elo_vals,Sol_vals];

%Perform PCA with three components
[PCA_coeffs,PCA_scores,~,~,PCA_percExplained] = ...
    pca(observationMatrix,...
    'NumComponents',3);

%Get maximal values
[maxVal,maxInd] = max(abs(PCA_coeffs),[],2);
[~,maxIndCluster] = max(PCA_coeffs(6,:));

maxValOP = max(abs(PCA_coeffs(1:3,:)),[],1);
[~,maxIndOP] = max(maxValOP);
maxIndOP = maxIndOP + (maxIndOP==maxIndCluster);

%Set order for PCA (according to max vol cluster)
PCA_order = [maxIndCluster,maxIndOP,...
    setdiff([1,2,3],[maxIndCluster,maxIndOP])];
if PCA_order(1) == PCA_order(2)
    PCA_order = PCA_order(2:end);
end
PCA_coeffs_plot = PCA_coeffs(:,PCA_order);

figure(1);
clf

subplot(1,4,1)
imagesc(PCA_coeffs_plot',[-1,+1])
%

PCA_labels = arrayfun(@(nn) ...
    sprintf('PC %d (%1.1f%%)',...
    PCA_order(nn),PCA_percExplained(PCA_order(nn))),1:3,...
    'UniformOutput',false);

set(gca,'YTick',1:3,'YTickLabel',PCA_labels)
set(gca,'XTick',1:8,'XTickLabel',...
    {'Distance',...
    'OP Ser5P','OP Ser2P','Clus Ser5P','Clus Ser2P',...
    'Clus Vol','Clus Elo','Clus Sol'})

% title(sprintf('f(d<%d nm)=%1.1f%%',...
%     dist_threshold.*1000,frac_close.*100),...
%     'FontWeight','normal')

% % 
% % 	
% %     %Define colormap matrix
% % 	colormap_matrix = ones(65,3);
% % 	colormap_matrix(1:33,2) = linspace(0,1,33);
% % 	colormap_matrix(1:33,3) = linspace(0,1,33);
% % 	colormap_matrix(33:65,1) = linspace(1,0,33);
% % 	colormap_matrix(33:65,2) = linspace(1,0,33);
% % 	colormap(colormap_matrix);
% % 	

%Order for first two components
PCA_coeffs = PCA_coeffs(:,PCA_order([1,2]));
PCA_scores = PCA_scores(:,PCA_order([1,2]));
PCA_percExplained = PCA_percExplained(PCA_order([1,2]));
PCA_labels = PCA_labels([1,2]);

topCount = 5;
[~,S5P_sortInds] = sort(OP_S5P_vals,'descend');
S5P_topInds = S5P_sortInds(1:topCount);
[~,Vol_sortInds] = sort(Vol_vals,'descend');
Vol_topInds = Vol_sortInds(1:topCount);

%Vector for S5P intensity
S5P_vec = [...
    mean(PCA_scores(S5P_topInds,1)),...
    mean(PCA_scores(S5P_topInds,2))];

%Vector for cluster volume
Vol_vec = [...
    mean(PCA_scores(Vol_topInds,1)),...
    mean(PCA_scores(Vol_topInds,2))];

subplot(1,4,2)
plot(PCA_scores(:,1),PCA_scores(:,2),'k.',...
    'MarkerEdgeColor',[0.6,0.6,0.6])
hold on
plot([0,S5P_vec(1)],[0,S5P_vec(2)],'m-','LineWidth',1)
plot([0,Vol_vec(1)],[0,Vol_vec(2)],'b-','LineWidth',1)
xlabel(PCA_labels{1})
ylabel(PCA_labels{2})
axis equal
%set(gca,'XLim',[-6,3],'YLim',[-2,8])

if S5P_vec(2)>0
    Vol_ortho_vec = Vol_vec*[0,-1;+1,0];
else
    Vol_ortho_vec = Vol_vec*[0,+1;-1,0];
end

%Perform transformation, such that volume points to north

%Define unit vectors
unit_vec_2 = Vol_vec./norm(Vol_vec);
unit_vec_1 = Vol_ortho_vec./norm(Vol_ortho_vec);

%Define transformation matrix
Trafo_matrix = [unit_vec_1',unit_vec_2'];

%Perform transformation
trafo_scores = PCA_scores*Trafo_matrix;
trafo_Vol_vec = Vol_vec*Trafo_matrix;
trafo_S5P_vec = S5P_vec*Trafo_matrix;


subplot(1,4,3)
plot(trafo_scores(:,1),trafo_scores(:,2),'k.',...
    'MarkerEdgeColor',[0.6,0.6,0.6])
hold on
plot([0,trafo_S5P_vec(1)],[0,trafo_S5P_vec(2)],'m-','LineWidth',1)
plot([0,trafo_Vol_vec(1)],[0,trafo_Vol_vec(2)],'b-','LineWidth',1)
xlabel('Transformed 1')
ylabel('Transformed 2')
axis equal
%set(gca,'XLim',[-6,3],'YLim',[-2,8])

angle_shift = -0.132;
%Get angles from x,y coordinate
angles = -atan2(trafo_scores(:,2),trafo_scores(:,1))./2./pi; % in radians
%apply modulo operator
angles = mod(angles,1);
%Perform angle shift
angles = mod(angles+angle_shift,1);
%Get number of angles
numPoints = numel(angles);
%Define pseudo-time coordinate
coord_s = ((1:numPoints)-1)./numPoints;

[angles,sortInds] = sort(angles);

%Calculate moving mean
windowSize = 20;
curveSmooth = @(xx) movmean(...
    padarray(xx,windowSize,'circular','both'),...
    windowSize,'Endpoints','discard');

plot(curveSmooth(trafo_scores(sortInds,1)),...
    curveSmooth(trafo_scores(sortInds,2)),'k-','LineWidth',1)

%Sort all variables
dist_vals = dist_vals(sortInds);
OP_S5P_vals = OP_S5P_vals(sortInds);
OP_S2P_vals = OP_S2P_vals(sortInds);
Cluster_S5P_vals = Cluster_S5P_vals(sortInds);
Cluster_S2P_vals = Cluster_S2P_vals(sortInds);
Vol_vals = Vol_vals(sortInds);
Elo_vals = Elo_vals(sortInds);
Sol_vals = Sol_vals(sortInds);
pseudoTime_central_slices{cc} = Central_slices(sortInds);

%Bin discretization for plotting of properties vs. pseudo time s

%Define bin borders
windowWidth = 0.15;
numWindows = 100;
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

square_ind = round(numWindows.*0.35);
square_marker = 'ks';
circle_ind = round(numWindows.*0.51);
circle_marker = 'ko';
diamond_ind = round(numWindows.*0.6);
diamond_marker = 'kd';

subplot(2,2,3)
plot(coord_s,dist_vals,'ko','MarkerFaceColor',[0.6,0.6,0.6],...
    'MarkerEdgeColor','none','MarkerSize',3)
hold on
plot(windowCenters,mean_dist,'k-','LineWidth',1)
plot(windowCenters(circle_ind),mean_dist(circle_ind),circle_marker,...
    'MarkerFaceColor',[1,0,0])
plot(windowCenters(square_ind),mean_dist(square_ind),square_marker,...
    'MarkerFaceColor',[0,0,1])
plot(windowCenters(diamond_ind),mean_dist(diamond_ind),diamond_marker,...
    'MarkerFaceColor',[0.3,0.3,0.3])
xlabel('Pseudo-time s')
ylabel('Distance [\mum]')
set(gca,'YLim',[0,4])

subplot(2,2,2)
plot(coord_s,OP_S5P_vals,'ko','MarkerFaceColor',[0.6,0.6,0.6],...
    'MarkerEdgeColor','none','MarkerSize',3)
hold on
plot(windowCenters,mean_OP_S5P,'k-','LineWidth',1)
plot(windowCenters(circle_ind),mean_OP_S5P(circle_ind),circle_marker,...
    'MarkerFaceColor',[1,0,0])
plot(windowCenters(square_ind),mean_OP_S5P(square_ind),square_marker,...
    'MarkerFaceColor',[0,0,1])
plot(windowCenters(diamond_ind),mean_OP_S5P(diamond_ind),diamond_marker,...
    'MarkerFaceColor',[0.3,0.3,0.3])
xlabel('Pseudo-time s')
ylabel('Gene Pol II Ser5P')
set(gca,'YLim',[0,6])

subplot(2,2,1)
plot(coord_s,OP_S2P_vals,'ko','MarkerFaceColor',[0.6,0.6,0.6],...
    'MarkerEdgeColor','none','MarkerSize',3)
hold on
plot(windowCenters,mean_OP_S2P,'k-','LineWidth',1)
legend_array = zeros(1,3);
legend_array(2) = ...
    plot(windowCenters(circle_ind),mean_OP_S2P(circle_ind),circle_marker,...
    'MarkerFaceColor',[1,0,0]);
legend_array(1) = ...
    plot(windowCenters(square_ind),mean_OP_S2P(square_ind),square_marker,...
    'MarkerFaceColor',[0,0,1]);
legend_array(3) = ...
    plot(windowCenters(diamond_ind),mean_OP_S2P(diamond_ind),diamond_marker,...
    'MarkerFaceColor',[0.3,0.3,0.3]);
xlabel('Pseudo-time s')
ylabel('Gene Pol II Ser2P')
set(gca,'YLim',[0.5,2.5])

legend(legend_array,{'Induced','Associated','Transcribing'})

% title(sprintf('%s, f(d<%d nm)=%1.1f%%',...
%     geneName,dist_threshold.*1000,frac_close.*100),...
%     'FontWeight','normal')

subplot(2,2,4)

scatter(OP_S5P_vals,dist_vals,12,OP_S2P_vals,'o','filled')
xlabel('Gene Pol II Ser5P')
ylabel('Distance [\mum]')
colormap(flipud(parula))
colormap(parula)
clim([0,1.5])
colorbar
set(gca,'Box','on')







hold on
plot(mean_OP_S5P,mean_dist,'k-','LineWidth',1)
plot(mean_OP_S5P(circle_ind),mean_dist(circle_ind),circle_marker,...
    'MarkerFaceColor',[1,0,0])
plot(mean_OP_S5P(square_ind),mean_dist(square_ind),square_marker,...
    'MarkerFaceColor',[0,0,1])
plot(mean_OP_S5P(diamond_ind),mean_dist(diamond_ind),diamond_marker,...
    'MarkerFaceColor',[0.3,0.3,0.3])



% patch([mean_OP_S5P nan],[mean_OP_S2P nan],[mean_dist nan],...
%     [mean_dist nan], 'edgecolor', 'interp');