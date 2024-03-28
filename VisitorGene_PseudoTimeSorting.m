
%% Principal Component based sorting of observations

load('ConditionSortedResults')

figure(4)
clf

sorted_central_slices = cell(1,numConds);

radius_median = zeros(1,numConds);
radius_CI = zeros(2,numConds);
in_range_perc = zeros(1,numConds);
radii_cell = cell(1,numConds);
crossCorr_vals = zeros(1,numConds);
crossCorr_CI = zeros(2,numConds);

Vol_threshold = 0.2;
angle_shift = -0.1;


for cc = 1:numConds
	
% 	figure(1)
	
	inclInds = ...
		sortedDistCell{cc}<=Inf & sortedVolCell{cc}>=Vol_threshold;
		
	dist_vals = [sortedDistCell{cc}(inclInds)];
	OP_S5P_vals = [sortedOPIntCell{1}{cc}(inclInds)];
	OP_S2P_vals = [sortedOPIntCell{2}{cc}(inclInds)];
	Cluster_S5P_vals = [sortedIntCell{1}{cc}(inclInds)];
	Cluster_S2P_vals = [sortedIntCell{2}{cc}(inclInds)];
	Vol_vals = [sortedVolCell{cc}(inclInds)];
	Elo_vals = [sortedEloCell{cc}(inclInds)];
	Sol_vals = [sortedSolCell{cc}(inclInds)];
	Central_slices = [sortedMaskCell{cc}(inclInds)];
	
	frac_close = sum(dist_vals<=dist_threshold)./numel(dist_vals);
	in_range_perc(cc) = frac_close.*100;

	[ff,xx] = ksdensity(dist_vals,'BandWidth',0.07,'support','unbounded');
	
	subplot(4,numConds,0.*numConds+cc)
	plot(xx,ff,'k-','LineWidth',1)
	hold on
	plot([1,1].*dist_threshold,[0,1],'k-','LineWidth',1)
	xlabel('Gene-cluster distance')
	ylabel('Probability')
	title(sprintf('%s (n=%d)',...
		sortedCondNames{cc},numel(sortedDistCell{cc})))
	set(gca,'XLim',[-0.05,4.5])
		
	
	% PCA, input: Rows of X are observations, columns to variables
		
	observationMatrix = [dist_vals,...
		OP_S5P_vals,OP_S2P_vals, ...
		Cluster_S5P_vals,Cluster_S2P_vals,...
		Vol_vals,Elo_vals,Sol_vals];
	
	[PCA_coeffs,PCA_scores,~,~,PCA_percExplained] = ...
		pca(observationMatrix,...
		'NumComponents',3);
	
	[maxVal,maxInd] = max(abs(PCA_coeffs),[],2);	
	[~,maxIndCluster] = max(PCA_coeffs(6,:));
	
	maxValOP = max(abs(PCA_coeffs(1:3,:)),[],1);
	[~,maxIndOP] = max(maxValOP);
	maxIndOP = maxIndOP + (maxIndOP==maxIndCluster);
	
	PCA_order = [maxIndCluster,maxIndOP,...
		setdiff([1,2,3],[maxIndCluster,maxIndOP])];
	if PCA_order(1) == PCA_order(2)
		PCA_order = PCA_order(2:end);
	end
	PCA_coeffs_plot = PCA_coeffs(:,PCA_order);
	
	subplot(4,numConds,1.*numConds+cc)
	imagesc(PCA_coeffs_plot',[-1,+1])
	
	PCA_labels = arrayfun(@(nn) ...
		sprintf('PC %d (%1.1f%%)',...
		PCA_order(nn),PCA_percExplained(PCA_order(nn))),1:3,...
		'UniformOutput',false);
	set(gca,'YTick',1:3,'YTickLabel',PCA_labels)
	
	title(sprintf('f(d<%d nm)=%1.1f%%',...
		dist_threshold.*1000,frac_close.*100),...
		'FontWeight','normal')

	
	colormap_matrix = ones(65,3);
	colormap_matrix(1:33,2) = linspace(0,1,33);
	colormap_matrix(1:33,3) = linspace(0,1,33);
	colormap_matrix(33:65,1) = linspace(1,0,33);
	colormap_matrix(33:65,2) = linspace(1,0,33);
	colormap(colormap_matrix);
	
	PCA_coeffs = PCA_coeffs(:,PCA_order([1,2]));
	PCA_scores = PCA_scores(:,PCA_order([1,2]));
	PCA_percExplained = PCA_percExplained(PCA_order([1,2]));
	PCA_labels = PCA_labels([1,2]);
	
	topCount = 10;
	[~,S5P_sortInds] = sort(OP_S5P_vals,'descend');
	S5P_topInds = S5P_sortInds(1:topCount);
	[~,Vol_sortInds] = sort(Vol_vals,'descend');
	Vol_topInds = Vol_sortInds(1:topCount);
	
	S5P_vec = [...
		mean(PCA_scores(S5P_topInds,1)),...
		mean(PCA_scores(S5P_topInds,2))];
	Vol_vec = [...
		mean(PCA_scores(Vol_topInds,1)),...
		mean(PCA_scores(Vol_topInds,2))];
	
	subplot(4,numConds,2.*numConds+cc)
	plot(PCA_scores(:,1),PCA_scores(:,2),'k.',...
		'MarkerEdgeColor',[0.6,0.6,0.6])
	hold on
	plot([0,S5P_vec(1)],[0,S5P_vec(2)],'m-','LineWidth',1)
	plot([0,Vol_vec(1)],[0,Vol_vec(2)],'b-','LineWidth',1)
	xlabel(PCA_labels{1})
	ylabel(PCA_labels{2})
	axis equal
% 	set(gca,'XLim',[-6,3],'YLim',[-2,8])

	if S5P_vec(2)>0
		Vol_ortho_vec = Vol_vec*[0,-1;+1,0];
	else
		Vol_ortho_vec = Vol_vec*[0,+1;-1,0];
	end
		
% 	plot([0,Vol_ortho_vec(1)],[0,Vol_ortho_vec(2)],'k-','LineWidth',1)
	
	unit_vec_2 = Vol_vec./norm(Vol_vec);
	unit_vec_1 = Vol_ortho_vec./norm(Vol_ortho_vec);

	Trafo_matrix = [unit_vec_1',unit_vec_2'];
	
	trafo_scores = PCA_scores*Trafo_matrix;
	trafo_Vol_vec = Vol_vec*Trafo_matrix;
	trafo_S5P_vec = S5P_vec*Trafo_matrix;	
	
	
	subplot(4,numConds,3.*numConds+cc)
	plot(trafo_scores(:,1),trafo_scores(:,2),'k.',...
		'MarkerEdgeColor',[0.6,0.6,0.6])
	hold on
	plot([0,trafo_S5P_vec(1)],[0,trafo_S5P_vec(2)],'m-','LineWidth',1)
	plot([0,trafo_Vol_vec(1)],[0,trafo_Vol_vec(2)],'b-','LineWidth',1)
	xlabel('Transformed 1')
	ylabel('Transformed 2')
	axis equal
	set(gca,'XLim',[-6,3],'YLim',[-2,8])
	
	angles = -atan2(trafo_scores(:,2),trafo_scores(:,1))./2./pi; % in radians
	angles = mod(angles,1);
	angles = mod(angles+angle_shift,1);
	numPoints = numel(angles);
	coord_s = ((1:numPoints)-1)./numPoints;
	
	[angles,sortInds] = sort(angles);
	
	windowSize = 20;
	curveSmooth = @(xx) movmean(...
		padarray(xx,windowSize,'circular','both'),...
		windowSize,'Endpoints','discard');
% 	plot(curveSmooth(trafo_scores(sortInds,1)),...
% 		curveSmooth(trafo_scores(sortInds,2)),'k-','LineWidth',1)
	
	radii = sqrt(sum(trafo_scores.^2,2));
	radius_median(cc) = median(radii);
	radius_CI(:,cc) = bootci(1000,@median,radii);
	radii_cell{cc} = radii;
		
	dist_vals = dist_vals(sortInds);
	OP_S5P_vals = OP_S5P_vals(sortInds);
	OP_S2P_vals = OP_S2P_vals(sortInds);
	Cluster_S5P_vals = Cluster_S5P_vals(sortInds);
	Cluster_S2P_vals = Cluster_S2P_vals(sortInds);
	Vol_vals = Vol_vals(sortInds);
	Elo_vals = Elo_vals(sortInds);
	Sol_vals = Sol_vals(sortInds);
	sorted_central_slices{cc} = Central_slices(sortInds);
	
	crossCorr_vals(cc) = corr(...
		OP_S5P_vals,circshift(Elo_vals,-ceil(numPoints.*0.2)));
	crossCorr_CI(:,cc) = bootci(1000,@corr,...
		OP_S5P_vals,circshift(Elo_vals,-ceil(numPoints.*0.2)));
	
%     crossCorr_vals(cc) = corr(...
% 		OP_S5P_vals,circshift(Cluster_S2P_vals,-ceil(numPoints.*0.25)));
% 	crossCorr_CI(:,cc) = bootci(1000,@corr,...
% 		OP_S5P_vals,circshift(Cluster_S2P_vals,-ceil(numPoints.*0.25)));
	
% 	title(sprintf('\\rho=%2.2f [%2.2f,%2.2f]',crossCorr_vals(cc),...
% 		crossCorr_CI(1,cc),crossCorr_CI(2,cc)),...
% 		'FontWeight','normal')
% 
% 	
% 	figure(2)
% 	
% 	% --- bin discretization
% 	
% 	windowWidth = 0.2;
% 	numWindows = 100;
% 	windowCenters = linspace(0,1,numWindows);
% 	leftEdges = windowCenters-windowWidth./2;
% 	rightEdges = windowCenters+windowWidth./2;
% 
% 		
% 	mean_dist = zeros(1,numWindows);
% 	mean_OP_S5P = zeros(1,numWindows);
% 	mean_OP_S2P = zeros(1,numWindows);
% 	mean_Cluster_S5P = zeros(1,numWindows);
% 	mean_Cluster_S2P = zeros(1,numWindows);
% 	mean_Vol = zeros(1,numWindows);
% 	mean_Elo = zeros(1,numWindows);
% 	mean_Sol = zeros(1,numWindows);
% 
% 	for nn = 1:numWindows
% 		
% 		thisLeftEdge = leftEdges(nn);
% 		thisRightEdge = rightEdges(nn);
% 		
% 		windowInds = find(coord_s>=thisLeftEdge & coord_s<thisRightEdge);
% 		if thisRightEdge>1
% 			windowInds = [windowInds,find(coord_s<(thisRightEdge-1))];
% 		end
% 		if thisLeftEdge<0
% 			windowInds = [windowInds,find(coord_s>(thisLeftEdge+1))];
% 		end
% 		
% 		mean_dist(nn) = median(dist_vals(windowInds));
% 		mean_OP_S5P(nn) = median(OP_S5P_vals(windowInds));
% 		mean_OP_S2P(nn) = median(OP_S2P_vals(windowInds));
% 		mean_Cluster_S5P(nn) = median(Cluster_S5P_vals(windowInds));
% 		mean_Cluster_S2P(nn) = median(Cluster_S2P_vals(windowInds));
% 		mean_Vol(nn) = median(Vol_vals(windowInds));
% 		mean_Elo(nn) = median(Elo_vals(windowInds));
% 		mean_Sol(nn) = median(Sol_vals(windowInds));
% 	end
% 	
% 	subplot(11,numConds,0.*numConds+cc)
% 	plot(coord_s,dist_vals,'ko','MarkerFaceColor',[0.6,0.6,0.6],...
% 		'MarkerEdgeColor','none','MarkerSize',3)
% 	hold on
% 	plot(windowCenters,mean_dist,'k-','LineWidth',1)
% 	title(sprintf('Gene: %s',sortedCondNames{cc}))
% 	if cc == 1
% 		ylabel('Distance [\mum]')
% 	end
% 	set(gca,'XTickLabel',[],'YLim',[0,4])
% 
% 	subplot(11,numConds,1.*numConds+cc)
% 	plot(coord_s,OP_S5P_vals,'ko','MarkerFaceColor',[0.6,0.6,0.6],...
% 		'MarkerEdgeColor','none','MarkerSize',3)
% 	hold on
% 	plot(windowCenters,mean_OP_S5P,'k-','LineWidth',1)
% 	if cc == 1
% 		ylabel('Gene S5P')
% 	end
% 	set(gca,'XTickLabel',[],'YLim',[0,6])
% 
% 	subplot(11,numConds,2.*numConds+cc)
% 	plot(coord_s,OP_S2P_vals,'ko','MarkerFaceColor',[0.6,0.6,0.6],...
% 		'MarkerEdgeColor','none','MarkerSize',3)
% 	hold on
% 	plot(windowCenters,mean_OP_S2P,'k-','LineWidth',1)
% 	if cc == 1
% 		ylabel('Gene S2P')
% 	end
% 	set(gca,'XTickLabel',[],'YLim',[0.5,2.5])
% 	
% 	subplot(11,numConds,3.*numConds+cc)
% 	plot(coord_s,Cluster_S5P_vals,'ko','MarkerFaceColor',[0.6,0.6,0.6],...
% 		'MarkerEdgeColor','none','MarkerSize',3)
% 	hold on
% 	plot(windowCenters,mean_Cluster_S5P,'k-','LineWidth',1)
% 	if cc == 1
% 		ylabel('Cluster S5P')
% 	end
% 	set(gca,'XTickLabel',[])
% 
% 	subplot(11,numConds,4.*numConds+cc)
% 	plot(coord_s,Cluster_S2P_vals,'ko','MarkerFaceColor',[0.6,0.6,0.6],...
% 		'MarkerEdgeColor','none','MarkerSize',3)
% 	hold on
% 	plot(windowCenters,mean_Cluster_S2P,'k-','LineWidth',1)
% 	if cc == 1
% 		ylabel('Cluster S2P')
% 	end
% 	set(gca,'XTickLabel',[],'YLim',[0.5,2.5])
% 	
% 	subplot(11,numConds,5.*numConds+cc)
% 	plot(coord_s,Vol_vals,'ko','MarkerFaceColor',[0.6,0.6,0.6],...
% 		'MarkerEdgeColor','none','MarkerSize',3)
% 	hold on
% 	plot(windowCenters,mean_Vol,'k-','LineWidth',1)
% 	if cc == 1
% 		ylabel('Volume [\mum^3]')
% 	end
% 	set(gca,'XTickLabel',[],'YLim',[0,10])
% 	
% 	subplot(11,numConds,6.*numConds+cc)
% 	plot(coord_s,Elo_vals,'ko','MarkerFaceColor',[0.6,0.6,0.6],...
% 		'MarkerEdgeColor','none','MarkerSize',3)
% 	hold on
% 	plot(windowCenters,mean_Elo,'k-','LineWidth',1)
% 	if cc == 1
% 		ylabel('Elongation')
% 	end
% 	set(gca,'XTickLabel',[],'YLim',[0,7])
% 	
% 	subplot(11,numConds,7.*numConds+cc)
% 	plot(coord_s,Sol_vals,'ko','MarkerFaceColor',[0.6,0.6,0.6],...
% 		'MarkerEdgeColor','none','MarkerSize',3)
% 	hold on
% 	plot(windowCenters,mean_Sol,'k-','LineWidth',1)
% 	if cc == 1
% 		ylabel('Solidity')
% 	end
% 	xlabel('Progress s')
% 	set(gca,'YLim',[0,1])
% 	
% 	% Cross-correlation analysis
% 	
% 	shiftRange = numPoints;
% 	shiftVals = (1:5:shiftRange)-1-ceil(numPoints./2);
% 	numShifts = numel(shiftVals);
% 	padRange = numPoints;
% 	n_boot = 10;
% 
% 	
% 	S5P_vals = padarray(Cluster_S5P_vals,padRange,'circular','both');
% 	S2P_vals = padarray(Cluster_S2P_vals,padRange,'circular','both');
% 	Sol_vals = padarray(Sol_vals,padRange,'circular','both');
% 	Elo_vals = padarray(Elo_vals,padRange,'circular','both');
% 	
% 	corr_S5P = zeros(1,numShifts);
% 	corr_S5P_CI = zeros(2,numShifts);
% 	corr_S2P = zeros(1,numShifts);
% 	corr_S2P_CI = zeros(2,numShifts);
% 	corr_Sol = zeros(1,numShifts);
% 	corr_Sol_CI = zeros(2,numShifts);
% 	shiftDist = zeros(1,numShifts);
% 
% 	for nn = 1:numShifts
% 		
% 		shiftInds = numPoints+(1:numPoints)+shiftVals(nn);
% 		
% 		corr_S5P(nn) = corr(...
% 			Elo_vals(numPoints+(1:numPoints)),...
% 			S5P_vals(shiftInds));
% 		corr_S5P_CI(:,nn) = bootci(n_boot,@corr,...
% 			Elo_vals(numPoints+(1:numPoints)),...
% 			S5P_vals(shiftInds));
% 		corr_S2P(nn) = corr(...
% 			Elo_vals(numPoints+(1:numPoints)),...
% 			S2P_vals(shiftInds));
% 		corr_S2P_CI(:,nn) = bootci(n_boot,@corr,...
% 			Elo_vals(numPoints+(1:numPoints)),...
% 			S2P_vals(shiftInds));
% 		corr_Sol(nn) = corr(...
% 			Elo_vals(numPoints+(1:numPoints)),...
% 			Sol_vals(shiftInds));
% 		corr_Sol_CI(:,nn) = bootci(n_boot,@corr,...
% 			Elo_vals(numPoints+(1:numPoints)),...
% 			Sol_vals(shiftInds));
% 		
% 	end
% 	
% 	windowSize = 1;
% 	curveSmooth = @(xx) movmean(...
% 		padarray(xx,windowSize,'circular','both'),...
% 		windowSize,'Endpoints','discard');
% 	
% 	subplot(11,numConds,8.*numConds+cc)
% 	plot(shiftVals([1,end]),[0,0],'k-','LineWidth',1,'Color',[0.5,0.5,0.5])
% 	hold on
% 	plot([0,0],[-1,+1],'k-','LineWidth',1,'Color',[0.5,0.5,0.5])
% 	plot(shiftVals,curveSmooth(corr_S5P),'k-','Linewidth',1)
% 	plot(shiftVals,curveSmooth(corr_S5P_CI(1,:)),'k-','Linewidth',0.7,...
% 		'Color',[0.6,0.6,0.6])
% 	plot(shiftVals,curveSmooth(corr_S5P_CI(2,:)),'k-','Linewidth',0.7,...
% 		'Color',[0.6,0.6,0.6])
% 	set(gca,'XLim',[-Inf,Inf],'YLim',[-1,+1],'XTickLabels',[])
% 	ylabel('Corr S5P-Elon.')
% 	
% 	subplot(11,numConds,9.*numConds+cc)
% 	plot(shiftVals([1,end]),[0,0],'k-','LineWidth',1,'Color',[0.5,0.5,0.5])
% 	hold on
% 	plot([0,0],[-1,+1],'k-','LineWidth',1,'Color',[0.5,0.5,0.5])
% 	plot(shiftVals,curveSmooth(corr_S2P),'k-','Linewidth',1)
% 	plot(shiftVals,curveSmooth(corr_S2P_CI(1,:)),'k-','Linewidth',0.7,...
% 		'Color',[0.6,0.6,0.6])
% 	plot(shiftVals,curveSmooth(corr_S2P_CI(2,:)),'k-','Linewidth',0.7,...
% 		'Color',[0.6,0.6,0.6])
% 	set(gca,'XLim',[-Inf,Inf],'YLim',[-1,+1],'XTickLabels',[])
% 	ylabel('Corr S2P-Elon.')
% 	
% 	subplot(11,numConds,10.*numConds+cc)
% 	plot(shiftVals([1,end]),[0,0],'k-','LineWidth',1,'Color',[0.5,0.5,0.5])
% 	hold on
% 	plot([0,0],[-1,+1],'k-','LineWidth',1,'Color',[0.5,0.5,0.5])
% 	plot(shiftVals,curveSmooth(corr_Sol),'k-','Linewidth',1)
% 	plot(shiftVals,curveSmooth(corr_Sol_CI(1,:)),'k-','Linewidth',0.7,...
% 		'Color',[0.6,0.6,0.6])
% 	plot(shiftVals,curveSmooth(corr_Sol_CI(2,:)),'k-','Linewidth',0.7,...
% 		'Color',[0.6,0.6,0.6])
% 	ylabel('Corr Sol.-Elon.')
% 	xlabel('Register shift')
% 	set(gca,'XLim',[-Inf,Inf],'YLim',[-1,+1])
	
end

figure(5)
clf

errorbar(in_range_perc,crossCorr_vals,...
	+crossCorr_CI(2,:)-crossCorr_vals,...
	-crossCorr_CI(1,:)+crossCorr_vals,'ko');
xlabel(sprintf('f(d<%d nm) [percent]',1000.*dist_threshold))
ylabel('1/4 shift correlation')
set(gca,'YLim',[-1,1])


%% Overview plot of changes with developmental stages

figure(1)
clf

figure(2)
clf

plotSets = {...
	{'Oblong_foxd5','Sphere_foxd5','Dome_foxd5','Epi30_foxd5','Epi50_foxd5'},...
	{'Oblong_ripply1','Sphere_ripply1','Dome_ripply1','Epi30_ripply1','Epi50_ripply1'},...
	{'Oblong_vamp2','Sphere_vamp2','Dome_vamp2','Epi30_vamp2','Epi50_vamp2'},...
	{'Oblong_drll2','Sphere_drll2','Dome_drll2','Epi30_drll2','Epi50_drll2'},...
	{'Oblong_klf2b','Sphere_klf2b','Dome_klf2b','Epi30_klf2b','Epi50_klf2b'},...
	{'Oblong_zgc:64022','Sphere_zgc:64022','Dome_zgc:64022','Epi30_zgc:64022','Epi50_zgc:64022'},...
	{'Oblong_gadd45ga','Sphere_gadd45ga','Dome_gadd45ga','Epi30_gadd45ga','Epi50_gadd45ga'},...
	{'Oblong_iscub','Sphere_iscub','Dome_iscub','Epi30_iscub','Epi50_iscub'}...
	};
setNames = {'foxd5','ripply1','vamp2','drll.2',...
    'klf2b','zgc:64022','gadd45ga','iscub'};
stageNames = {'Obl','Sph','Dome','Epi30','Epi50'};
setColors = {[0,0,1],[1,0,0],[0.5,0.5,0.5]};
numSets = numel(plotSets);

sorted_central_slices = cell(1,numConds);

for nn_set = 1:numSets
	
	condNamesInSet = plotSets{nn_set};
	
	numCondsInSet = numel(condNamesInSet);
	
	in_range_perc = zeros(1,numCondsInSet);
	in_range_perc_CI = zeros(2,numCondsInSet);
	dist_median = zeros(1,numCondsInSet);
	dist_CI = zeros(2,numCondsInSet);
	S5P_median = zeros(1,numCondsInSet);
	S5P_CI = zeros(2,numCondsInSet);
	S2P_median = zeros(1,numCondsInSet);
	S2P_CI = zeros(2,numCondsInSet);
    S2P_top = zeros(1,numCondsInSet);
    S2P_top_CI = zeros(2,numCondsInSet);
    S2P_perc = zeros(1,numCondsInSet);
    S2P_perc_CI = zeros(2,numCondsInSet);
    S5P_perc = zeros(1,numCondsInSet);
    S5P_perc_CI = zeros(2,numCondsInSet);

    setCrossCorr = zeros(1,numCondsInSet);
    setCrossCorr_CI = zeros(2,numCondsInSet);

	dist_value_array = [];
	S5P_value_array = [];
	S2P_value_array = [];
	grouping_array = [];

	
	for cc = 1:numCondsInSet
		
		thisCondInd = find(...
			cellfun(@(elmt)strcmp(elmt,plotSets{nn_set}{cc}),sortedCondNames));
		
        setCrossCorr(cc) = crossCorr_vals(thisCondInd);
        setCrossCorr_CI(:,cc) = crossCorr_CI(:,thisCondInd);

        top_prctile = 5;
		dist_threshold = 0.3;
		Vol_threshold = 0.02;0.4;0.02;
		inclInds = ...
			sortedDistCell{thisCondInd}<=Inf ...
			& sortedVolCell{thisCondInd}>=Vol_threshold;
		
		dist_vals = [sortedDistCell{thisCondInd}(inclInds)];
		
		OP_S5P_vals = [sortedOPIntCell{1}{thisCondInd}(inclInds)];
		OP_S2P_vals = [sortedOPIntCell{2}{thisCondInd}(inclInds)];
		Cluster_S5P_vals = [sortedIntCell{1}{thisCondInd}(inclInds)];
		Cluster_S2P_vals = [sortedIntCell{2}{thisCondInd}(inclInds)];
		Vol_vals = [sortedVolCell{thisCondInd}(inclInds)];
		Elo_vals = [sortedEloCell{thisCondInd}(inclInds)];
		Sol_vals = [sortedSolCell{thisCondInd}(inclInds)];
		% 	Central_slices = [sortedCentralSliceCell{thisCondInd}(inclInds)];
		
        figure(2)
        
        subplot(1,numSets,nn_set)
        
%         plot3(dist_vals,OP_S5P_vals,OP_S2P_vals,'k.')
%         hold on
%         set(gca,'XLim',[0,2],'YLim',[0,7],'ZLim',[0,3])
%         ylabel('Pol II Ser5P')
%         xlabel('Distance [\mum]')
%         zlabel('Pol II Ser2P')
        
        plot(dist_vals,OP_S5P_vals,'k.')
        hold on
        
        set(gca,'XLim',[0,5],'YLim',[0,6])
        
        
        
        ylabel('Pol II Ser5P')
        xlabel('Distance [\mum]')
        
        title(setNames{nn_set})

        
		n_boot = 1000;
		dist_median(cc) = median(dist_vals);
		dist_CI(:,cc) = bootci(n_boot,@median,dist_vals);
		S5P_median(cc) = median(OP_S5P_vals);
		S5P_CI(:,cc) = bootci(n_boot,@median,OP_S5P_vals);
		S2P_median(cc) = median(OP_S2P_vals);
		S2P_CI(:,cc) = bootci(n_boot,@median,OP_S2P_vals);
		
        S2P_top(cc) = prctile(OP_S2P_vals,100-top_prctile);
        S2P_top_CI(:,cc) = bootci(n_boot,...
            @(xx)prctile(xx,100-top_prctile),OP_S2P_vals);
        
% 		numPoints = numel(dist_vals);
% 		dist_value_array = [dist_value_array,dist_vals'];
% 		S5P_value_array = [S5P_value_array,OP_S5P_vals'];
% 		S2P_value_array = [S2P_value_array,OP_S2P_vals'];
% 		grouping_array = [grouping_array,ones(1,numPoints).*cc];
		
		frac_close = sum(dist_vals<=dist_threshold)./numel(dist_vals);
		in_range_perc(cc) = frac_close.*100;
		
		in_range_perc(cc) = 100.*...
			sum(dist_vals<=dist_threshold)./numel(dist_vals);
		in_range_perc_CI(:,cc) = 100./numel(dist_vals).*...
			bootci(n_boot,@(dd)sum(dd<=dist_threshold),dist_vals);
		
		S2P_perc(cc) = 100.*...
			sum(OP_S2P_vals>=1.5)./numel(OP_S2P_vals);
		S2P_perc_CI(:,cc) = 100./numel(OP_S2P_vals).*...
			bootci(n_boot,@(xx)sum(xx>=1.5),OP_S2P_vals);
		
        S5P_perc(cc) = 100.*...
			sum(OP_S5P_vals>=2.0)./numel(OP_S5P_vals);
		S5P_perc_CI(:,cc) = 100./numel(OP_S5P_vals).*...
			bootci(n_boot,@(xx)sum(xx>=2.0),OP_S5P_vals);
		
        
	end
	
	figure(1)
	
    subplot(4,numSets,0.*numSets+nn_set)
% 	errorbar(1:numCondsInSet,dist_median,...
% 		dist_CI(1,:)-dist_median,...
% 		dist_median-dist_CI(2,:),...
% 		'k-o')
% 	ylabel('OP-Cl d [\mum]')

	errorbar(1:numCondsInSet,in_range_perc,...
		in_range_perc_CI(1,:)-in_range_perc,...
		in_range_perc-in_range_perc_CI(2,:),...
		'k-o')
    if nn_set == 1
        ylabel('Contact %')
    else
        set(gca,'YTickLabel',[])
    end
    
	set(gca,'XTick',1:numCondsInSet,'XTickLabel',stageNames,...
		'XLim',([1,numCondsInSet])+[-0.5,+0.5],...
		'YLim',[0,25])

    title(setNames{nn_set})

    
	subplot(4,numSets,numSets+nn_set)
% 	errorbar(1:numCondsInSet,S5P_median,...
% 		S5P_CI(1,:)-S5P_median,...
% 		S5P_median-S5P_CI(2,:),...
% 		'k-o')
    errorbar(1:numCondsInSet,S5P_perc,...
		S5P_perc_CI(1,:)-S5P_perc,...
		S5P_perc-S5P_perc_CI(2,:),...
		'k-o')
	if nn_set == 1
        ylabel('S5P %')
    else
        set(gca,'YTickLabel',[])
    end
    
	set(gca,'XTick',1:numCondsInSet,'XTickLabel',stageNames,...
		'XLim',([1,numCondsInSet])+[-0.5,+0.5],...
		'YLim',[0,30])
	

	subplot(4,numSets,2.*numSets+nn_set)
% 	errorbar(1:numCondsInSet,S2P_median,...
% 		S2P_CI(1,:)-S2P_median,...
% 		S2P_median-S2P_CI(2,:),...
% 		'k-o')
    errorbar(1:numCondsInSet,S2P_perc,...
		S2P_perc_CI(1,:)-S2P_perc,...
		S2P_perc-S2P_perc_CI(2,:),...
		'k-o')
    if nn_set == 1
        ylabel('S2P %')
    else
        set(gca,'YTickLabel',[])
    end
	
	set(gca,'XTick',1:numCondsInSet,'XTickLabel',stageNames,...
		'XLim',([1,numCondsInSet])+[-0.5,+0.5],...
		'YLim',[0,25])

    
    subplot(4,numSets,3.*numSets+nn_set)
% 	errorbar(1:numCondsInSet,S2P_median,...
% 		S2P_CI(1,:)-S2P_median,...
% 		S2P_median-S2P_CI(2,:),...
% 		'k-o')
    errorbar(1:numCondsInSet,setCrossCorr,...
		setCrossCorr_CI(1,:)-setCrossCorr,...
		setCrossCorr-setCrossCorr_CI(2,:),...
		'k-o')
    if nn_set == 1
        ylabel('Clock radius')
    else
        set(gca,'YTickLabel',[])
    end
	
	set(gca,'XTick',1:numCondsInSet,'XTickLabel',stageNames,...
		'XLim',([1,numCondsInSet])+[-0.5,+0.5],...
		'YLim',[0,2])


	
% 	figure(2)
%     
%     subplot(1,2,2)
% 	plot(in_range_perc,S2P_perc,'o-')
% 	hold on
% 	
% 	xlabel('Contact %')
% 	ylabel('Pol II Ser2P %')
				
end

% figure(2)
% 
% subplot(1,2,1)
% legend(setNames)
% 
% subplot(1,2,2)
% legend(setNames)