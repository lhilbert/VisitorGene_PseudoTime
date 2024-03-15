clear all

sourceDirectory = './ExtractedStacks_Stages/Epi_30/Cond_25/';
sourceDirectory = './ExtractedStacks_Stages/**/';
sourceDirectory = './ExtractedStacks_Stages/Sphere/';

% Channels for segmentation
NucSegChannel = 3;
ClusterSegChannel = 3;
OPSegChannel = 1;

ImgSquareSize = 24; % pixels for cut-out images

quantChannels = [3,2];
quantBlurSigma = [0,0];

storeImgChannels = [1,3,2];
numStoreChannels = numel(storeImgChannels);

segBlurSigma_small = 1.0; % in microns
segBlurSigma_large = 10; % in microns

clusterSegBlurSigma_large = 3.0;
clusterSeg_numStdDev = 2.0;%5;

OPsegBlurSigma_small = 0.1; % in microns
OPsegBlurSigma_large = 5; % in microns
OPseg_numStdDev = 6; % number of standard deviations in robust threshold

Nuc_min_vol = 40; % cubic microns
Nuc_min_sol = 0.7; % to ensure round nuclei
Cluster_minVol = 0.03;
OP_minVol = 0.05; % cubic microns

listing = rdir([sourceDirectory,'*Image*.mat']);

numFiles = numel(listing);

% Condition index retrieval
condInds = [];
condNames = {};
for ff = 1:numFiles
	thisFilePath = listing(ff).name;
	thisCondInd = load(thisFilePath,'condInd');
	thisCondInd = thisCondInd.condInd;
	condInds = [condInds,thisCondInd];
	thisCondName = load(thisFilePath,'condName');
	thisCondName = thisCondName.condName;
	condNames = [condNames,thisCondName];
end

% analyze image stacks one by one

numNuclei_vec = zeros(1,numFiles);

Cluster_distCell = cell(1,numFiles);
Cluster_volCell = cell(1,numFiles);
Cluster_solCell = cell(1,numFiles);
Cluster_eloCell = cell(1,numFiles);
CentralSlicesCell = cell(1,numFiles);
Cluster_centCell = cell(1,numFiles);
Cluster_intCell = cell(1,numFiles);
OP_intCell = cell(1,numFiles);
OP_displCell = cell(1,numFiles);

numQuantChannels = numel(quantChannels);

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
	
	% --- For each nucleus, get clusters and oligopaints
	
	Cluster_img = imgStack{ClusterSegChannel};
	OP_img = imgStack{OPSegChannel};
	
	Cluster_OP_dist = cell(1,numNuclei);
	Cluster_volume = cell(1,numNuclei);
	Cluster_solidity = cell(1,numNuclei);
	Cluster_elongation = cell(1,numNuclei);
	CentralSlices_store = cell(1,numNuclei);
	Cluster_intensity = cell(numQuantChannels,numNuclei);
	ClusterCent_store = cell(1,numNuclei);
	OP_intensity = cell(numQuantChannels,numNuclei);
	
	for nn = 1:numNuclei
		
		boxArray = props.BoundingBox(nn,:);
		Cluster_subImage = imgStack{NucSegChannel}(...
			boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
			boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
			boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
		OP_subImage = OP_img(...
			boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
			boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
			boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
		
		store_subImages = cell(1,numStoreChannels);
		for color = 1:numStoreChannels
			store_subImages{color} = imgStack{storeImgChannels(color)}(...
				boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
				boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
				boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
		end
		
		Cluster_subImage = Cluster_subImage ...
			- imgaussfilt(Cluster_subImage,...
			clusterSegBlurSigma_large./pixelSize);		
		
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
		
		Cluster_intensities = Cluster_subImage(NucMask_subImage);
		Cluster_mean = mean(Cluster_intensities);
		Cluster_std = std(Cluster_intensities);
		Cluster_mask = (Cluster_subImage.*NucMask_subImage)...
			>(Cluster_mean+clusterSeg_numStdDev.*Cluster_std);
		
		OP_intensities = OP_subImage(NucMask_subImage);
		OP_mean = mean(OP_intensities);
		OP_std = std(OP_intensities);
		OP_mask = (OP_subImage.*NucMask_subImage)...
			>(OP_mean+OPseg_numStdDev.*OP_std);
		
		subImgSize = size(Cluster_subImage);
		if numel(subImgSize)==2
			subImgSize(3)=1;
		else
			subplot(2,2,1)
			imagesc(squeeze(max(Cluster_subImage,[],3)))
			axis tight equal
			
			subplot(2,2,2)
			imagesc(squeeze(max(OP_subImage,[],3)))
			axis tight equal
			
			subplot(2,2,3)
			imagesc(squeeze(max(Cluster_mask,[],3)))
			axis tight equal

			subplot(2,2,4)
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
		
		Cluster_comps = bwconncomp(Cluster_mask,18);
		Cluster_numPxls = cellfun(@(elmt)numel(elmt),Cluster_comps.PixelIdxList);
		minPixels = Cluster_minVol./(pixelSize.^2)./zStepSize;
		Cluster_comps.NumObjects = sum(Cluster_numPxls>minPixels);
		Cluster_comps.PixelIdxList = Cluster_comps.PixelIdxList(Cluster_numPxls>minPixels);
		Cluster_numPxls = cellfun(@(elmt)numel(elmt),Cluster_comps.PixelIdxList);
		
		
		if OP_comps.NumObjects>0 && Cluster_comps.NumObjects>0
		
			OP_props = regionprops3(OP_comps,Cluster_subImage,...
				'Centroid','VoxelIdxList');
			
			Cluster_props = regionprops3(Cluster_comps,Cluster_subImage,...
				'Volume','MeanIntensity','Solidity',...
				'Centroid','Image','BoundingBox');
			
			pwDist = pdist2(...
				OP_props.Centroid.*[pixelSize,pixelSize,zStepSize],...
				Cluster_props.Centroid.*[pixelSize,pixelSize,zStepSize]);
			[minDist,minInds] = min(pwDist,[],2);
			
			Cluster_OP_minDist = minDist;
			Cluster_OP_minInd = minInds;
			
			Cluster_props = Cluster_props(Cluster_OP_minInd,:);
			Cluster_comps.NumObjects = numel(Cluster_OP_minDist);
			Cluster_comps.PixelIdxList = ...
				Cluster_comps.PixelIdxList(Cluster_OP_minInd);
			
			Cluster_Volume_array = ...
				[Cluster_props.Volume].*pixelSize.^2.*zStepSize;
			Cluster_Solidity_array = [Cluster_props.Solidity];
			Cluster_centroid_array = ...
				Cluster_props.Centroid.*[pixelSize,pixelSize,zStepSize];
			
			% --- get cluster central plane and elongation in-plane
			Cluster_Elongation_array = ...
				zeros(size(Cluster_Solidity_array));
			Slices_cell = ...
				cell(size(Cluster_Solidity_array));
			for cl = 1:numel(Cluster_Volume_array)
				boundingBox = Cluster_props.BoundingBox(cl,:);
				thisImage = squeeze(Cluster_subImage(...
					boundingBox(2)+0.5:boundingBox(2)+boundingBox(5)-0.5,...
					boundingBox(1)+0.5:boundingBox(1)+boundingBox(4)-0.5,...
					ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));
				thisImage = thisImage-min(thisImage(:));
				thisImage = thisImage./max(thisImage(:));
				
				thisMask = Cluster_props.Image{cl};
				centerInd = ceil(size(thisMask,3)./2);
				thisMask = squeeze(thisMask(:,:,centerInd));
				thisImage((bwperim(thisMask))) = 0;
				Slices_cell{cl} = cell(1,numStoreChannels+2);
				for color = 1:numStoreChannels
					Slices_cell{cl}{color} = ...
						squeeze(store_subImages{color}(:,:,...
						ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));
				end
				Slices_cell{cl}{numStoreChannels+1} = ...
					squeeze(Cluster_mask(:,:,...
					ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));
				Slices_cell{cl}{numStoreChannels+2} = ...
					squeeze(OP_mask(:,:,...
					ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));

				thisProps = regionprops(uint8(thisMask),...
					'MajorAxisLength','MinorAxisLength');
				Cluster_Elongation_array(cl) = ...
					thisProps.MajorAxisLength./thisProps.MinorAxisLength;
			end
			Cluster_OP_dist{nn} = Cluster_OP_minDist;
			Cluster_volume{nn} = Cluster_Volume_array;
			Cluster_solidity{nn} = Cluster_Solidity_array;
			Cluster_elongation{nn} = Cluster_Elongation_array;
			CentralSlices_store{nn} = Slices_cell;
			ClusterCent_store{nn} = Cluster_centroid_array;

			% --- quantification for all target channels
			for qq = 1:numQuantChannels
				
				channelInd = quantChannels(qq);
				quant_subImage = imgStack{channelInd}(...
					boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
					boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
					boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
				
				if quantBlurSigma(qq)>0
					quant_subImage = imgaussfilt(quant_subImage,...
						quantBlurSigma(qq)./pixelSize);
				end
				
				Quant_nucleusMedian = ...
					median(quant_subImage(NucMask_subImage));
				
				Cluster_quant_props = regionprops3(...
					Cluster_comps,quant_subImage,'VoxelValues');
				Quant_ClusterMedian = cellfun(...
					@(vals)median(vals),Cluster_quant_props.VoxelValues);
				Cluster_intensity{qq,nn} = ...
					Quant_ClusterMedian./Quant_nucleusMedian;
				
				OP_quant_props = regionprops3(...
					OP_comps,quant_subImage,'VoxelValues');
				Quant_OPMedian = cellfun(...
					@(vals)median(vals),OP_quant_props.VoxelValues);
				OP_intensity{qq,nn} = ...
					Quant_OPMedian./Quant_nucleusMedian;

			end
			
		else
			Cluster_OP_dist{nn} = [];
			Cluster_volume{nn} = [];
			Cluster_solidity{nn} = [];
			Cluster_elongation{nn} = [];
			CentralSlices_store{nn} = {};
			ClusterCent_store{nn} = [];
			for qq = 1:numQuantChannels
				Cluster_intensity{qq,nn} = [];
				OP_intensity{qq,nn} = [];
			end
			
		end
		
	end
	
	Cluster_distCell{ff} = vertcat(Cluster_OP_dist{:});
	Cluster_volCell{ff} = vertcat(Cluster_volume{:});
	Cluster_solCell{ff} = vertcat(Cluster_solidity{:});
	Cluster_eloCell{ff} = vertcat(Cluster_elongation{:});
	CentralSlicesCell{ff} = vertcat(CentralSlices_store{:});
	Cluster_centCell{ff} = vertcat(ClusterCent_store{:});
	Cluster_intCell{ff} = cell(1,numQuantChannels);
	for qq = 1:numQuantChannels
		Cluster_intCell{ff}{qq} = vertcat(Cluster_intensity{qq,:});
	end
	OP_intCell{ff} = cell(1,numQuantChannels);
	for qq = 1:numQuantChannels
		OP_intCell{ff}{qq} = vertcat(OP_intensity{qq,:});
	end
	
end

%% Sort into conditions

uniqueCondNames = unique(condNames);
numConds = numel(uniqueCondNames);
fileIndsCell = cell(1,numConds);
numFiles_perCond = zeros(1,numConds);
for cc = 1:numConds
	fileIndsCell{cc} = cellfun(...
		@(elmt)strcmp(elmt,uniqueCondNames{cc}),condNames);
	numFiles_perCond(cc) = sum(fileIndsCell{cc});
end

sortedCondNames = cell(1,numConds);
sortedNumNuclei = zeros(1,numConds);
sortedNumFiles = zeros(1,numConds);
sortedDistCell = cell(1,numConds);
sortedVolCell = cell(1,numConds);
sortedSolCell = cell(1,numConds);
sortedEloCell = cell(1,numConds);
sortedCentralSliceCell = cell(1,numConds);
sortedCentroidsCell = cell(1,numConds);
sortedIntCell = cell(1,numQuantChannels);
sortedOPIntCell = cell(1,numQuantChannels);
for qq = 1:numQuantChannels
	sortedIntCell{qq} = cell(1,numConds);
	sortedOPIntCell{qq} = cell(1,numConds);
end

for cc = 1:numConds
	
	sortedCondNames{cc} = ...
		condNames(fileIndsCell{cc});
	sortedCondNames{cc} = sortedCondNames{cc}{1};
	
	sortedNumNuclei(cc) = ...
		sum(numNuclei_vec(fileIndsCell{cc}));
	sortedNumFiles(cc) = sum(fileIndsCell{cc});
	
	Cluster_dists = vertcat(Cluster_distCell{fileIndsCell{cc}});
	Cluster_vols = vertcat(Cluster_volCell{fileIndsCell{cc}});
	Cluster_sols = vertcat(Cluster_solCell{fileIndsCell{cc}});
	Cluster_elos = vertcat(Cluster_eloCell{fileIndsCell{cc}});
	Central_slices = vertcat(CentralSlicesCell{fileIndsCell{cc}});
	Cluster_centroids = vertcat(Cluster_centCell{fileIndsCell{cc}});

	Cluster_ints = Cluster_intCell(fileIndsCell{cc});
	OP_ints = OP_intCell(fileIndsCell{cc});

	sortedDistCell{cc} = Cluster_dists;
	sortedVolCell{cc} = Cluster_vols;
	sortedSolCell{cc} = Cluster_sols;
	sortedEloCell{cc} = Cluster_elos;
	sortedCentralSliceCell{cc} = Central_slices;
	sortedCentroidsCell{cc} = Cluster_centroids;

	for qq = 1:numQuantChannels
		
		sortedOPIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)OP_intCell{ind}{qq},....
			find(fileIndsCell{cc}),...
			'UniformOutput',false));
		sortedOPIntCell{qq}{cc} = vertcat(sortedOPIntCell{qq}{cc}{:});
		sortedIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)Cluster_intCell{ind}{qq},....
			find(fileIndsCell{cc}),...
			'UniformOutput',false));
		sortedIntCell{qq}{cc} = vertcat(sortedIntCell{qq}{cc}{:});

	end
	
end


%% Principal Component based sorting of observations

figure(1)
clf

in_range_perc = zeros(1,numConds);
in_range_perc_CI = zeros(2,numConds);
dist_median = zeros(1,numConds);
dist_CI = zeros(2,numConds);
S5P_median = zeros(1,numConds);
S5P_CI = zeros(2,numConds);
S2P_median = zeros(1,numConds);
S2P_CI = zeros(2,numConds);

dist_value_array = [];
S5P_value_array = [];
S2P_value_array = [];
grouping_array = [];

for cc = 1:numConds
	
	dist_threshold = 0.75;
	Vol_threshold = 0.02;
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
	
	numPoints = numel(dist_vals);
	median_dist(cc) = median(dist_vals);
	
	n_boot = 1000;
	dist_median(cc) = median(dist_vals);
	dist_CI(:,cc) = bootci(n_boot,@median,dist_vals);
	S5P_median(cc) = median(OP_S5P_vals);
	S5P_CI(:,cc) = bootci(n_boot,@median,OP_S5P_vals);
	S2P_median(cc) = median(OP_S2P_vals);
	S2P_CI(:,cc) = bootci(n_boot,@median,OP_S2P_vals);
		
	frac_close = sum(dist_vals<=dist_threshold)./numel(dist_vals);
	in_range_perc(cc) = frac_close.*100;
	
	in_range_perc(cc) = 100.*...
		sum(dist_vals<=dist_threshold)./numel(dist_vals);
	in_range_perc_CI(:,cc) = 100./numel(dist_vals).*...
		bootci(n_boot,@(dd)sum(dd<=dist_threshold),dist_vals);
				
end

plotInds = {[1:numConds]};
numPlotSets = numel(plotInds);
setSymbols = {'ko','ro','ro','ks'};
setFaceColors = {[0,0,0],[1,0,0],'none','none'};
setNames = {'All'};

for pp = 1:numPlotSets
	
	subplot(4,1,1)
	errorbar(plotInds{pp},S5P_median(plotInds{pp}),...
		S5P_CI(1,(plotInds{pp}))-S5P_median((plotInds{pp})),...
		S5P_median((plotInds{pp}))-S5P_CI(2,(plotInds{pp})),...
		setSymbols{pp},'MarkerFaceColor',setFaceColors{pp})
	hold on
	
	set(gca,'XTick',1:numConds,'XTickLabels','')
	xlabel('')
	ylabel('Pol II Ser5P')
	
	
	subplot(4,1,2)
	errorbar(plotInds{pp},S2P_median(plotInds{pp}),...
		S2P_CI(1,(plotInds{pp}))-S2P_median((plotInds{pp})),...
		S2P_median((plotInds{pp}))-S2P_CI(2,(plotInds{pp})),...
		setSymbols{pp},'MarkerFaceColor',setFaceColors{pp})
	hold on
	set(gca,'XTick',1:numConds,'XTickLabels','')
	xlabel('')
	ylabel('Pol II Ser2P')
	
	
	
	subplot(4,1,3)
	errorbar(plotInds{pp},in_range_perc(plotInds{pp}),...
		in_range_perc_CI(1,plotInds{pp})-in_range_perc(plotInds{pp}),...
		in_range_perc(plotInds{pp})-in_range_perc_CI(2,plotInds{pp}),...
		setSymbols{pp},'MarkerFaceColor',setFaceColors{pp})
	hold on
	set(gca,'XTick',1:numConds,'XTickLabels','')
	xlabel('')
	ylabel(sprintf('%% with d<%d nm',dist_threshold.*1000))
	
	
	subplot(4,1,4)
	errorbar(plotInds{pp},dist_median(plotInds{pp}),...
		dist_CI(1,plotInds{pp})-dist_median(plotInds{pp}),...
		dist_median(plotInds{pp})-dist_CI(2,plotInds{pp}),...
		setSymbols{pp},'MarkerFaceColor',setFaceColors{pp})
	hold on
	set(gca,'XTick',1:numConds,'XTickLabels',sortedCondNames)
	xlabel('')
	ylabel('Distance [\mum]')
	xtickangle(45)

end

subplot(4,1,1)

legend(setNames)