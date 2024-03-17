clear all

sourceDirectory = './ExtractedStacks/**/';

% Channels for segmentation
NucSegChannel = 2;
ClusterSegChannel = 2;
OPSegChannel = 3;

quantChannels = [1,2,3];
quantBlurSigma = [0,0,0.07];

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

centralSliceExtension = 0; % pixels from centroid

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
uniqueCondInds = unique(condInds);
numConds = numel(uniqueCondInds);
numFiles_perCond = arrayfun(@(nn)sum(condInds==nn),uniqueCondInds);

% analyze image stacks one by one

numNuclei_vec = zeros(1,numFiles);

Cluster_distCell = cell(1,numFiles);
Cluster_volCell = cell(1,numFiles);
Cluster_solCell = cell(1,numFiles);
Cluster_eloCell = cell(1,numFiles);
Cluster_maskCell = cell(1,numFiles);
Cluster_intCell = cell(1,numFiles);
OP_intCell = cell(1,numFiles);
OP_displCell = cell(1,numFiles);
pixelSize_array = zeros(1,numFiles);
zStep_array = zeros(1,numFiles);

validFileFlags = false(1,numFiles);

numQuantChannels = numel(quantChannels);

for ff = 1:numFiles
%parfor ff = 1:numFiles
	
	fprintf('Processing file %d of %d\n',ff,numFiles)
	
	thisCondInd = condInds(ff);	
	thisFilePath = listing(ff).name;
	
	loadStruct = load(thisFilePath,...
		'imgStack','imgSize','pixelSize','zStepSize');
	imgStack = loadStruct.imgStack;
	imgSize = loadStruct.imgSize;
	pixelSize = loadStruct.pixelSize;
	zStepSize = loadStruct.zStepSize;

	pixelSize_array(ff) = pixelSize;
	zStep_array(ff) = zStepSize;
	
	% Nuclei segmentation
	segImg = imgStack{NucSegChannel};
	segImg = ...
		+ imgaussfilt(segImg,segBlurSigma_small./pixelSize) ...
		- imgaussfilt(segImg,segBlurSigma_large./pixelSize);

	[bin_counts,bin_centers] = hist(segImg(:),1000);
	[nuc_seg_thresh,~] = otsuLimit(bin_centers,bin_counts,[0,Inf]);
	NucSegMask = segImg>1.0.*nuc_seg_thresh;
	
    figure(1)
    clf

	subplot(2,3,1)
	imagesc(squeeze(imgStack{NucSegChannel}(:,:,ceil(imgSize(3)./2))))
	axis tight equal
	
	subplot(2,3,2)
	imagesc(squeeze(segImg(:,:,ceil(imgSize(3)./2))))
	axis tight equal

	subplot(2,3,3)
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

    if comps.NumObjects == 0
        % If no valid nuclei are deteced, record empty vectors

        Volume_array = [];
    	Intensity_array = [];
    	Solidity_array = [];
    	numNuclei = 0;

    else
        % If no valid nuclei are detected, record their properties

    	comps.PixelIdxList = comps.PixelIdxList(Solidity_array>=Nuc_min_sol);
    	numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList);

    	props = regionprops3(comps,imgStack{NucSegChannel},...
    		'Volume','VoxelValues','Solidity','VoxelIdxList',...
    		'BoundingBox');

    	Volume_array = [props.Volume].*pixelSize.^2.*zStepSize;
    	Intensity_array = cellfun(@(vals)median(vals),props.VoxelValues);
    	Solidity_array = [props.Solidity];
    	numNuclei = numel(Intensity_array);

    end
	
	numNuclei_vec(ff) = numNuclei;
	
	% --- For each nucleus, get clusters and oligopaints
    % The case of no nuclei detected is automatically covered by the
    % following code. If there are no nuclei, veftors and cells of length 0
    % will be saved, and no iteration will be carried out, so that the case
    % is treated accurately.
	
	Cluster_img = imgStack{ClusterSegChannel};
	OP_img = imgStack{OPSegChannel};
	
	Cluster_OP_dist = cell(1,numNuclei);
	Cluster_volume = cell(1,numNuclei);
	Cluster_solidity = cell(1,numNuclei);
	Cluster_elongation = cell(1,numNuclei);
	Cluster_mask_store = cell(1,numNuclei);
	Cluster_intensity = cell(numQuantChannels,numNuclei);
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

            figure(2)
            clf
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
			
			% --- get cluster central plane and elongation in-plane
			Cluster_Elongation_array = ...
				zeros(size(Cluster_Solidity_array));
			for cl = 1:numel(Cluster_Volume_array)
				
				thisMask = Cluster_props.Image{cl};
				centerInd = ceil(size(thisMask,3)./2);
				thisMask = squeeze(thisMask(:,:,centerInd));
				thisProps = regionprops(uint8(thisMask),...
					'MajorAxisLength','MinorAxisLength');
				Cluster_Elongation_array(cl) = ...
					thisProps.MajorAxisLength./thisProps.MinorAxisLength;
				
			end
			Cluster_OP_dist{nn} = Cluster_OP_minDist;
			Cluster_volume{nn} = Cluster_Volume_array;
			Cluster_solidity{nn} = Cluster_Solidity_array;
			Cluster_elongation{nn} = Cluster_Elongation_array;

			% --- quantification for all target channels
			Cluster_Mask_cell = ...
				cell(size(Cluster_Solidity_array));
			for cl = 1:numel(Cluster_Volume_array)
				Cluster_Mask_cell{cl} = cell(1,numQuantChannels);
			end
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
				
				if centralSliceExtension==0
					
					for cl = 1:numel(Cluster_Volume_array)
						Cluster_Mask_cell{cl}{qq} = [];
					end
					
				else
					
					for cl = 1:numel(Cluster_Volume_array)
						
						centroid_1 = round(Cluster_props.Centroid(cl,2));
						centroid_2 = round(Cluster_props.Centroid(cl,1));
						center_z = round(Cluster_props.Centroid(cl,3));
						% 					boundingBox = Cluster_props.BoundingBox(cl,:);
						% 					center_z = ...
						% 						ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1));
						img_limits = [...
							centroid_1-centralSliceExtension,...
							centroid_1+centralSliceExtension,...
							centroid_2-centralSliceExtension,...
							centroid_2+centralSliceExtension];
						img_limits = [...
							max(1,img_limits(1)),...
							min(subImgSize(1),img_limits(2)),...
							max(1,img_limits(3)),...
							min(subImgSize(2),img_limits(4))];
						
						thisImage = squeeze(quant_subImage(...
							img_limits(1):img_limits(2),...
							img_limits(3):img_limits(4),...
							center_z))./Quant_nucleusMedian;
						
						% thisImage = thisImage-min(thisImage(:));
						% thisImage = thisImage./max(thisImage(:));
						
						% Where the example image gets stored:
						Cluster_Mask_cell{cl}{qq} = thisImage;
						
					end
				end
			end
			Cluster_mask_store{nn} = Cluster_Mask_cell;
						
		else
			Cluster_OP_dist{nn} = [];
			Cluster_volume{nn} = [];
			Cluster_solidity{nn} = [];
			Cluster_elongation{nn} = [];
			Cluster_mask_store{nn} = {};
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
	Cluster_maskCell{ff} = vertcat(Cluster_mask_store{:});
	Cluster_intCell{ff} = cell(1,numQuantChannels);
	for qq = 1:numQuantChannels
		Cluster_intCell{ff}{qq} = vertcat(Cluster_intensity{qq,:});
	end
	OP_intCell{ff} = cell(1,numQuantChannels);
	for qq = 1:numQuantChannels
		OP_intCell{ff}{qq} = vertcat(OP_intensity{qq,:});
    end

    figure(1)



    
    subplot(2,3,4)
    % Count of OP detections per nucleus

    OP_countPerNucleus = ...
        cellfun(@(vals)numel(vals),Cluster_OP_dist);
    binCounts = histcounts(OP_countPerNucleus,[-0.5:1:12.5]);

    fractionAboveFour = ...
        sum(OP_countPerNucleus>4)...
        ./numel(OP_countPerNucleus);

    bar(0:12,binCounts,'FaceColor',[0.5,0.5,0.5],'EdgeColor',[0,0,0],...
        'LineWidth',1)
    set(gca,'XTick',0:12)
    xlabel('FISH spots per nucleus')
    ylabel('Count')
    title(sprintf('Fraction above 4 spots: %2.2f (%s)',...
        fractionAboveFour,condNames{ff}))
   




    subplot(2,3,5)
    % Cross talk assessment within OP foci
    dist_threshold = 0.25;
    inContact_inds = Cluster_distCell{ff}<=dist_threshold;
    plot(OP_intCell{ff}{ClusterSegChannel}(inContact_inds),...
        OP_intCell{ff}{OPSegChannel}(inContact_inds),'ko',...
        'MarkerEdgeColor','none','MarkerFaceColor',[1,0,0])
    hold on
    
    plot(OP_intCell{ff}{ClusterSegChannel},...
        OP_intCell{ff}{OPSegChannel},...
        'ko')

    hold on
    
    xlabel('Pol II Ser5P Intensity (a.u.)')
    ylabel('OP Intensity (a.u.)')
    title(sprintf('%d OP foci from %d nuclei',...
        numel(Cluster_distCell{ff}),numNuclei_vec(ff)))
    %set(gca,'XLim',[0,3],'YLim',[0,3])
    hold on
    plot([0,3],[0,3],'k-','LineWidth',1.5)
    hold off

    subplot(2,3,6)
    % Cross talk assessment within Pol II Ser5P clusters
    plot(Cluster_intCell{ff}{ClusterSegChannel}(inContact_inds),...
        Cluster_intCell{ff}{OPSegChannel}(inContact_inds),'ko',...
        'MarkerEdgeColor','none','MarkerFaceColor',[1,0,0])
    hold on
    
    plot(Cluster_intCell{ff}{ClusterSegChannel},...
        Cluster_intCell{ff}{OPSegChannel},...
        'ko')

    hold on
    
    xlabel('Pol II Ser5P Intensity (a.u.)')
    ylabel('OP Intensity (a.u.)')
    title(sprintf('Nearest neighbor Pol II clusters',...
        numel(Cluster_distCell{ff}),numNuclei_vec(ff)))
    %set(gca,'XLim',[0,3],'YLim',[0,3])
    hold on
    plot([0,3],[0,3],'k-','LineWidth',1.5)
    hold off


    if numNuclei_vec(ff)>0
        waitforbuttonpress
    end

    if fractionAboveFour > 0.15
        validFileFlags(ff) = false
    else
        validFileFlags(ff) = true;
    end
	
end

%% Sort into conditions

sortedCondNames = cell(1,numConds);
sortedNumNuclei = zeros(1,numConds);
sortedNumFiles = zeros(1,numConds);
sortedDistCell = cell(1,numConds);
sortedVolCell = cell(1,numConds);
sortedSolCell = cell(1,numConds);
sortedEloCell = cell(1,numConds);
sortedMaskCell = cell(1,numConds);
sortedIntCell = cell(1,numQuantChannels);
sortedOPIntCell = cell(1,numQuantChannels);
for qq = 1:numQuantChannels
	sortedIntCell{qq} = cell(1,numConds);
	sortedOPIntCell{qq} = cell(1,numConds);
end

sortedPixelSize = zeros(1,numConds);
sortedZStep = zeros(1,numConds);

for cc = 1:numConds
	
	sortedCondNames{cc} = ...
		condNames(condInds==uniqueCondInds(cc));
	sortedCondNames{cc} = sortedCondNames{cc}{1};
	
	sortedNumNuclei(cc) = ...
		sum(numNuclei_vec(condInds==uniqueCondInds(cc)));
	sortedNumFiles(cc) = sum(condInds==uniqueCondInds(cc));
	
	Cluster_dists = vertcat(Cluster_distCell{condInds==uniqueCondInds(cc)});
	Cluster_vols = vertcat(Cluster_volCell{condInds==uniqueCondInds(cc)});
	Cluster_sols = vertcat(Cluster_solCell{condInds==uniqueCondInds(cc)});
	Cluster_elos = vertcat(Cluster_eloCell{condInds==uniqueCondInds(cc)});
	Cluster_masks = vertcat(Cluster_maskCell{condInds==uniqueCondInds(cc)});
	Cluster_ints = Cluster_intCell(condInds==uniqueCondInds(cc));
	OP_ints = OP_intCell(condInds==uniqueCondInds(cc));

	sortedDistCell{cc} = Cluster_dists;
	sortedVolCell{cc} = Cluster_vols;
	sortedSolCell{cc} = Cluster_sols;
	sortedEloCell{cc} = Cluster_elos;
	sortedMaskCell{cc} = Cluster_masks;

	pixelSizes = pixelSize_array(condInds==uniqueCondInds(cc));
	zStepSizes = zStep_array(condInds==uniqueCondInds(cc));
	
	sortedPixelSize(cc) = pixelSizes(1);
	sortedZStep(cc) = zStepSizes(1);

	for qq = 1:numQuantChannels
		
		sortedOPIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)OP_intCell{ind}{qq},....
			find(condInds==uniqueCondInds(cc)),...
			'UniformOutput',false));
		sortedOPIntCell{qq}{cc} = vertcat(sortedOPIntCell{qq}{cc}{:});
		sortedIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)Cluster_intCell{ind}{qq},....
			find(condInds==uniqueCondInds(cc)),...
			'UniformOutput',false));
		sortedIntCell{qq}{cc} = vertcat(sortedIntCell{qq}{cc}{:});

	end
	
end


%% Principal Component based sorting of observations

% figure(1)
% clf
% 
% figure(2)
% clf

sorted_central_slices = cell(1,numConds);

radius_median = zeros(1,numConds);
radius_CI = zeros(2,numConds);
in_range_perc = zeros(1,numConds);
radii_cell = cell(1,numConds);
crossCorr_vals = zeros(1,numConds);
crossCorr_CI = zeros(2,numConds);

dist_threshold = 0.25;
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
	
% 	subplot(4,numConds,0.*numConds+cc)
% 	plot(xx,ff,'k-','LineWidth',1)
% 	hold on
% 	plot([1,1].*dist_threshold,[0,1],'k-','LineWidth',1)
% 	xlabel('Gene-cluster distance')
% 	ylabel('Probability')
% 	title(sprintf('%s (n=%d)',...
% 		sortedCondNames{cc},numel(sortedDistCell{cc})))
% 	set(gca,'XLim',[-0.05,4.5])
		
	
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
	
% 	subplot(4,numConds,1.*numConds+cc)
% 	imagesc(PCA_coeffs_plot',[-1,+1])
% 	
% 	PCA_labels = arrayfun(@(nn) ...
% 		sprintf('PC %d (%1.1f%%)',...
% 		PCA_order(nn),PCA_percExplained(PCA_order(nn))),1:3,...
% 		'UniformOutput',false);
% 	set(gca,'YTick',1:3,'YTickLabel',PCA_labels)
% 	
% 	title(sprintf('f(d<%d nm)=%1.1f%%',...
% 		dist_threshold.*1000,frac_close.*100),...
% 		'FontWeight','normal')

	
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
	
% 	subplot(4,numConds,2.*numConds+cc)
% 	plot(PCA_scores(:,1),PCA_scores(:,2),'k.',...
% 		'MarkerEdgeColor',[0.6,0.6,0.6])
% 	hold on
% 	plot([0,S5P_vec(1)],[0,S5P_vec(2)],'m-','LineWidth',1)
% 	plot([0,Vol_vec(1)],[0,Vol_vec(2)],'b-','LineWidth',1)
% 	xlabel(PCA_labels{1})
% 	ylabel(PCA_labels{2})
% 	axis equal
% % 	set(gca,'XLim',[-6,3],'YLim',[-2,8])

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
	
	
% 	subplot(4,numConds,3.*numConds+cc)
% 	plot(trafo_scores(:,1),trafo_scores(:,2),'k.',...
% 		'MarkerEdgeColor',[0.6,0.6,0.6])
% 	hold on
% 	plot([0,trafo_S5P_vec(1)],[0,trafo_S5P_vec(2)],'m-','LineWidth',1)
% 	plot([0,trafo_Vol_vec(1)],[0,trafo_Vol_vec(2)],'b-','LineWidth',1)
% 	xlabel('Transformed 1')
% 	ylabel('Transformed 2')
% 	axis equal
% 	set(gca,'XLim',[-6,3],'YLim',[-2,8])
	
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

figure(3)
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