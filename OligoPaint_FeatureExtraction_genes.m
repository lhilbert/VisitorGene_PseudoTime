clear all

sourceDirectory = './ExtractedStacks/**/';

% Channels for segmentation
NucSegChannel = 3;
OPSegChannel = 1;

QuantChannels = [3,2];

segBlurSigma_small = 1.0; % in microns
segBlurSigma_large = 10; % in microns
OPsegBlurSigma_small = 0.1; % in microns
OPsegBlurSigma_large = 5; % in microns
OPseg_numStdDev = 6; % number of standard deviations in roobust threshold

Nuc_min_vol = 40; % cubic microns
Nuc_min_sol = 0.7; % to ensure round nuclei 
OP_minVol = 0.05; % cubic microns

listing = rdir([sourceDirectory,'*Image*.mat']);

numFiles = numel(listing);

% Condition index retrieval
condInds = [];
for ff = 1:numFiles
	thisFilePath = listing(ff).name;
	thisCondInd = load(thisFilePath,'condInd');
	thisCondInd = thisCondInd.condInd;
	condInds = [condInds,thisCondInd];
end
uniqueCondInds = unique(condInds);
numConds = numel(uniqueCondInds);
numFiles_perCond = arrayfun(@(nn)sum(condInds==nn),uniqueCondInds);

% analyze image stacks one by one

numNuclei_vec = zeros(1,numFiles);

OP_volCell = cell(1,numFiles);
OP_intCell = cell(1,numFiles);
OP_solCell = cell(1,numFiles);

numQuantChannels = numel(QuantChannels);
OP_quantCell = cell(numQuantChannels,numFiles);

parfor ff = 1:numFiles
	
	fprintf('Processing file %d of %d\n',ff,numFiles)
	
	thisCondInd = condInds(ff);	
	thisFilePath = listing(ff).name;
	
	loadStruct = load(thisFilePath,...
		'imgStack','imgSize','pixelSize','zStepSize');
	imgStack = loadStruct.imgStack;
	imgSize = loadStruct.imgSize;
	pixelSize = loadStruct.pixelSize;
	zStepSize = loadStruct.zStepSize;

	% Nuclei segmentation
	segImg = imgStack{NucSegChannel};
	segImg = ...
		+ imgaussfilt(segImg,segBlurSigma_small./pixelSize) ...
		- imgaussfilt(segImg,segBlurSigma_large./pixelSize);

	[bin_counts,bin_centers] = hist(segImg(:),1000);
	[nuc_seg_thresh,~] = otsuLimit(bin_centers,bin_counts,[0,Inf]);
	NucSegMask = segImg>1.0.*nuc_seg_thresh;
	
	subplot(1,3,1)
	imagesc(squeeze(imgStack{NucSegChannel}(:,:,ceil(imgSize(3)./2))))
	axis tight equal
	
	subplot(1,3,2)
	imagesc(squeeze(segImg(:,:,ceil(imgSize(3)./2))))
	axis tight equal

	subplot(1,3,3)
	imagesc(squeeze(NucSegMask(:,:,ceil(imgSize(3)./2))))
	axis tight equal
	
% 	waitforbuttonpress
	
	% --- Connected component segmentation of nuclei
	comps = bwconncomp(NucSegMask,18);
	numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList);
	minPixels = Nuc_min_vol./(pixelSize.^2)./zStepSize;
	comps.NumObjects = sum(numPxls>=minPixels);
	comps.PixelIdxList = comps.PixelIdxList(numPxls>=minPixels);
	numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList);
	
	props = regionprops3(comps,imgStack{NucSegChannel},...
		'Solidity');
	
	Solidity_array = [props.Solidity];
	
	comps.NumObjects = sum(Solidity_array>=Nuc_min_sol);
	comps.PixelIdxList = comps.PixelIdxList(Solidity_array>=Nuc_min_sol);
	numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList);
	
	props = regionprops3(comps,imgStack{NucSegChannel},...
		'Volume','VoxelValues','Solidity','VoxelIdxList',...
		'BoundingBox');
	
	Volume_array = [props.Volume].*pixelSize.^2.*zStepSize;
	Intensity_array = cellfun(@(vals)median(vals),props.VoxelValues);
	Solidity_array = [props.Solidity];
	
	
	numNuclei = numel(Intensity_array);
	
	numNuclei_vec(ff) = numNuclei;
	
	% --- For each nucleus, get oligopaint features
	
	OP_img = imgStack{OPSegChannel};
	OP_volume = cell(1,numNuclei);
	OP_solidity = cell(1,numNuclei);
	OP_intensity = cell(numQuantChannels,numNuclei);
	
	for nn = 1:numNuclei
		
		boxArray = props.BoundingBox(nn,:);
		Nuc_subImage = imgStack{NucSegChannel}(...
			boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
			boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
			boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
		OP_subImage = OP_img(...
			boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
			boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
			boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
		OP_subImage = imgaussfilt(OP_subImage,...
			OPsegBlurSigma_small./pixelSize);
		OP_subImage = OP_subImage - imgaussfilt(OP_subImage,...
			OPsegBlurSigma_large./pixelSize);
		NucMask = false(imgSize);
		NucMask(props.VoxelIdxList{nn}) = true;
		NucMask_subImage = NucMask(...
			boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
			boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
			boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
		
		
		OP_intensities = OP_subImage(NucMask_subImage);
		OP_mean = mean(OP_intensities);
		OP_std = std(OP_intensities);
		OP_mask = (OP_subImage.*NucMask_subImage)...
			>(OP_mean+OPseg_numStdDev.*OP_std);
		
		subImgSize = size(Nuc_subImage);
		if numel(subImgSize)==2
			subImgSize(3)=1;
		else
			subplot(1,3,1)
			imagesc(squeeze(max(Nuc_subImage,[],3)))
			axis tight equal
			
			subplot(1,3,2)
			imagesc(squeeze(max(OP_subImage,[],3)))
			axis tight equal
			
			subplot(1,3,3)
			imagesc(squeeze(max(OP_mask,[],3)))
			axis tight equal

		end
				
% 		waitforbuttonpress
		
		OP_comps = bwconncomp(OP_mask,18);
		OP_numPxls = cellfun(@(elmt)numel(elmt),OP_comps.PixelIdxList);
		minPixels = OP_minVol./(pixelSize.^2)./zStepSize;
		OP_comps.NumObjects = sum(OP_numPxls>minPixels);
		OP_comps.PixelIdxList = OP_comps.PixelIdxList(OP_numPxls>minPixels);
		OP_numPxls = cellfun(@(elmt)numel(elmt),OP_comps.PixelIdxList);
		
		if OP_comps.NumObjects>0
		
			OP_props = regionprops3(OP_comps,Nuc_subImage,...
				'Volume','VoxelValues','Solidity');
			
			OP_Volume_array = [OP_props.Volume].*pixelSize.^2.*zStepSize;
			OP_Median_array = cellfun(@(vals)median(vals),OP_props.VoxelValues);
			OP_Solidity_array = [OP_props.Solidity];
									
			OP_Median_array = ...
				OP_Median_array./Intensity_array(nn);
			
			OP_volume{nn} = OP_Volume_array;
			OP_solidity{nn} = OP_Solidity_array;

			% --- quantification for all target channels
			for qq = 1:numQuantChannels
				channelInd = QuantChannels(qq);
				quant_subImage = imgStack{channelInd}(...
					boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
					boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
					boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
		
				Quant_nucleusMedian = ...
					median(quant_subImage(NucMask_subImage));
				OP_quant_props = regionprops3(...
					OP_comps,quant_subImage,'VoxelValues');
				Quant_OPMedian = cellfun(...
					@(vals)median(vals),OP_quant_props.VoxelValues);
				OP_intensity{qq,nn} = ...
					Quant_OPMedian./Quant_nucleusMedian;

			end
			
		else
			
			OP_volume{nn} = [];
			OP_solidity{nn} = [];
			for qq = 1:numQuantChannels
				OP_intensity{qq,nn} = [];
			end
			
		end
		
	end
	
	OP_volCell{ff} = vertcat(OP_volume{:});
	OP_solCell{ff} = vertcat(OP_solidity{:});
	OP_intCell{ff} = cell(1,numQuantChannels);
	for qq = 1:numQuantChannels
		OP_intCell{ff}{qq} = vertcat(OP_intensity{qq,:});
	end
	
end

%% Sort into conditions and plot histograms

figure(1)
clf

sortedNumNuclei = zeros(1,numConds);

sortedVolCell = cell(1,numConds);
sortedSolCell = cell(1,numConds);
sortedIntCell = cell(1,numQuantChannels);
for qq = 1:numQuantChannels
	sortedIntCell{qq} = cell(1,numConds);
end

for cc = 1:numConds
	
	OP_vols = vertcat(OP_volCell{condInds==uniqueCondInds(cc)});
	OP_sols = vertcat(OP_solCell{condInds==uniqueCondInds(cc)});
	OP_ints = OP_intCell(condInds==uniqueCondInds(cc));
	
	sortedVolCell{cc} = OP_vols;
	sortedSolCell{cc} = OP_sols;
	
	sortedNumNuclei(cc) = ...
		sum(numNuclei_vec(condInds==uniqueCondInds(cc)));

	numSubplots = 2+numQuantChannels;
	
	subplot(1,numSubplots,1)
	[counts,locs] = hist(OP_vols,12);
	plot(locs,counts,'-')
	hold on
	xlabel('Foci Volumue [\mum^2]')
	ylabel('Count')
		
	subplot(1,numSubplots,2)
	[counts,locs] = hist(OP_sols,12);
	plot(locs,counts,'-')
	hold on
	xlabel('Foci solidity')
	ylabel('Count')

	for qq = 1:numQuantChannels
		
		subplot(1,numSubplots,2+qq)
		int_vals = vertcat(OP_ints{qq});
		int_vals = cellfun(@(elmt)elmt{qq},OP_ints,'UniformOutput',false);
		int_vals = vertcat(int_vals{:});
	
		sortedIntCell{qq}{cc} = int_vals;
				
		[pp,xx] = ksdensity(int_vals);
		plot(xx,pp,'-');
		
		hold on
		xlabel('Intensity at foci')
		ylabel('Count')
		
		title(sprintf('Channel %d',qq))
		
	end
	
end

subplot(1,numSubplots,1)

legend(arrayfun(@(xx)sprintf('Condition %d',xx),1:numConds,...
	'UniformOutput',false))


%% Boxplots

figure(2)
clf

compareCells = {sortedVolCell,sortedSolCell,...
	sortedIntCell{1},sortedIntCell{2}}; 
value_name_cell = {...
	'Area [\mum^2]','Solidity',...
	'Pol II Ser5P','Pol II Ser2P'};
title_name_cell = {...
	'','','',''};
naming_cell = {...
	'cdc25b','celf1','crsp7','drll2','foxd5',...
	'gadd45ga','iscub','klf2b','lft1',...
	'ripply1','rnf19a','vamp2','zgc'};
numComparisons = numel(compareCells);

n_cell = cell(1,numComparisons);
diff_cell = cell(1,numComparisons);
p_values_cell = cell(1,numComparisons);

diff_cell_mean = cell(1,numComparisons);
p_values_cell_mean = cell(1,numComparisons);

group_median = zeros(numConds,numComparisons);
group_CI_low = zeros(numConds,numComparisons);
group_CI_high = zeros(numConds,numComparisons);

n_permute = 10000;

bonf_Factor = 1;

for cc = 1:numComparisons
	
	thisCell = compareCells{cc};
	numConds = numel(thisCell);
	
	n_cell{cc} = zeros(1,numConds);
	
	diff_cell{cc} = zeros(numConds,numConds);
	p_values_cell{cc} = zeros(numConds,numConds);
	
	diff_cell_mean{cc} = zeros(numConds,numConds);
	p_values_cell_mean{cc} = zeros(numConds,numConds);
	
	% for boxplots
	value_array = [];
	grouping_array = [];
	
	% for error bar plots
	median_vec = zeros(1,numConds);
	CI_vec = zeros(2,numConds);
	
	for kk = 1:numConds
		
		allAVals = thisCell{kk};
		AInclInds = isfinite(allAVals);
		AVals = allAVals(AInclInds);
		n_A = numel(AVals);
		n_cell{cc}(kk) = n_A;
		
		value_array = [value_array,AVals'];
		grouping_array = [grouping_array,ones(1,n_A).*kk];
		
		median_vec(kk) = median(AVals(:));
		CI_vec(:,kk) = bootci(n_permute,@median,AVals(:));

		group_median(kk,cc) = median_vec(kk);
		group_CI_low(kk,cc) = CI_vec(2,kk)-median_vec(kk);
		group_CI_high(kk,cc) = CI_vec(1,kk)-median_vec(kk);

% 		for ll = 1:numConds
% 			
% 			allBVals = thisCell{ll};
% 			BInclInds = isfinite(allBVals);
% 			BVals = allBVals(BInclInds);
% 			n_B = numel(BVals);
% 			
% 			n_joint = n_A+n_B;
% 			jointVals = [AVals;BVals];
% 			diffMedians = median(BVals)-median(AVals);
% 			diffMeans = mean(BVals)-mean(AVals);
% 			
% 			HnullMedians = zeros(1,n_permute);
% 			HnullMeans = zeros(1,n_permute);
% 			
% 			for pp = 1:n_permute
% 				
% 				permutedVals = jointVals(randperm(n_joint));
% 				HnullMedians(pp) = median(permutedVals(1:n_A)) ...
% 					- median(permutedVals((n_A+1):end));
% 				HnullMeans(pp) = mean(permutedVals(1:n_A)) ...
% 					- mean(permutedVals((n_A+1):end));
% 				
% 			end
% 			
% 			% two-tailed p value
% 			p_value = (1+sum(HnullMedians>=max(diffMedians,-diffMedians)) ...
% 				+sum(HnullMedians<=min(diffMedians,-diffMedians)))./(1+n_permute);
% 			
% 			p_value_mean = (1+sum(HnullMeans>=max(diffMeans,-diffMeans)) ...
% 				+sum(HnullMeans<=min(diffMeans,-diffMeans)))./(1+n_permute);
% 			
% 			diff_cell{cc}(kk,ll) = diffMedians;
% 			p_values_cell{cc}(kk,ll) = p_value.*bonf_Factor;
% 			
% 			diff_cell_mean{cc}(kk,ll) = diffMeans;
% 			p_values_cell_mean{cc}(kk,ll) = p_value_mean.*bonf_Factor;
% 			
% 		end
		
	end
	
	subplot(2,numComparisons,cc)
	plot([0.5,numConds+0.5],[1,1],'k-')
	hold on
	boxplot(value_array,grouping_array,'labels',naming_cell,...
		'symbol','o','notch','on')
	xlabel('')
	ylabel(value_name_cell{cc})
	
	title(title_name_cell{cc})

	subplot(2,numComparisons,numComparisons+cc)
	plot([0.5,numConds+0.5],[1,1],'k-')
	hold on
	errorbar(1:numConds,median_vec,...
		CI_vec(2,:)-median_vec,CI_vec(1,:)-median_vec,'ko')
	xlabel('')
	ylabel(value_name_cell{cc})
	set(gca,'XLim',[0.5,numConds+0.5])
	
	set(gca,'XTick',[1:numConds],'XTickLabel',naming_cell)
	
	title(title_name_cell{cc})

	
end

%% Make 2D plot

figure(3)

clf

plot([0.9,1.2],[1,1],'k-','LineWidth',1,'Color',[0.6,0.6,0.6])
hold on
plot([1,1],[0.9,1.2],'k-','LineWidth',1,'Color',[0.6,0.6,0.6])

plotInds = [5,8,4,13];%[5,8,4,13];
active_genes = errorbar(group_median(plotInds,4),group_median(plotInds,3),...
	group_CI_low(plotInds,3),group_CI_high(plotInds,3),...
	group_CI_low(plotInds,4),group_CI_high(plotInds,4),...
	'bd','Linewidth',1,'MarkerFaceColor',[0,0,1],...
	'MarkerEdgeColor','none','Color',[0.7,0.7,1]);
plotInds = [6,7,10,12];
inactive_genes = errorbar(group_median(plotInds,4),group_median(plotInds,3),...
	group_CI_low(plotInds,3),group_CI_high(plotInds,3),...
	group_CI_low(plotInds,4),group_CI_high(plotInds,4),...
	'rd','Linewidth',1,'MarkerFaceColor',[1,0,0],...
	'MarkerEdgeColor','none','Color',[1,0.7,0.7]);

% low_SE = errorbar(group_median(7:9,4),group_median(7:9,3),...
% 	group_CI_low(7:9,3),group_CI_high(7:9,3),...
% 	group_CI_low(7:9,4),group_CI_high(7:9,4),...
% 	'rd','Linewidth',1,'MarkerFaceColor',[1,0,0],...
% 	'MarkerEdgeColor','none','Color',[1,0.7,0.7]);
% Epi_SE = errorbar(group_median(10:12,4),group_median(10:12,3),...
% 	group_CI_low(10:12,3),group_CI_high(10:12,3),...
% 	group_CI_low(10:12,4),group_CI_high(10:12,4),...
% 	'ro','Linewidth',1,'MarkerFaceColor','none',...
% 	'MarkerEdgeColor',[1,0,0],'Color',[1,0.7,0.7]);

xlabel('Pol II Ser2P')
ylabel('Pol II Ser5P')

legend([active_genes,inactive_genes],{'Active genes','Inactive genes'})

figure(4)
clf

for cc = 1:numConds
	
	subplot(2,ceil(numConds./2),cc)
	plot([0,4],[1,1],'k-','LineWidth',1,'Color',[0.6,0.6,0.6])
	hold on
	plot([1,1],[0,6],'k-','LineWidth',1,'Color',[0.6,0.6,0.6])
	plot(sortedIntCell{2}{cc},sortedIntCell{1}{cc},'k.')
	xlabel('Pol II Ser2P')
	ylabel('Pol II Ser5P')
	set(gca,'XLim',[0,4],'YLim',[0,6])
	title(sprintf('%s (n=%d,N=%d)',naming_cell{cc},...
		numel(sortedIntCell{2}{cc}),sortedNumNuclei(cc)),...
		'FontWeight','normal')
	
end