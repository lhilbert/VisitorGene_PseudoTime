clear all

sourceDirectory = './ExtractedStacks_Stages/**/';

% Channels for segmentation
NucSegChannel = 2;
ClusterSegChannel = 2;
OPSegChannel = 3;

% The channels must be assigned so that the they represent the following
% labels in the biological sample:
% Channel 1: Elongating Pol II, Pol II Ser2Phos
% Channel 2: Recruited Pol II, Pol II Ser5Phos
% Channel 3: Oligopaint DNA-FISH labeling of the gene of interest
quantChannels = [1,2,3];
quantBlurSigma = [0,0,0.07];

segBlurSigma_small = 1.0; % in microns
segBlurSigma_large = 10; % in microns

clusterSegBlurSigma_large = 0.5;
clusterSeg_numStdDev = 3.0;%5;

OPsegBlurSigma_small = 0.07; % in microns
OPsegBlurSigma_large = 0.3; % in microns
OPseg_numStdDev = 5; % number of standard deviations in robust threshold

Nuc_min_vol = 40; % cubic microns
Nuc_min_sol = 0.7; % to ensure round nuclei
Cluster_minVol = 0.08; % to only include large clusters
OP_minVol = 0.05; % cubic microns

dist_threshold = 0.5;0.25; % contact distance in micrometers

centralSliceExtension = 0; % pixels from centroid


% ------end of analysis parameters

% Get number of parallel workers

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end



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

%for ff = 1:numFiles
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
        % If valid nuclei are detected, record their properties

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
        y1=boxArray(2).*pixelSize;
        y2=(boxArray(2)+boxArray(5)).*pixelSize;
        x1=boxArray(1).*pixelSize;
        x2=(boxArray(1)+boxArray(4)).*pixelSize;

        showZSlice = round(boxArray(3)+0.5.*boxArray(6));
        if showZSlice==0
            showZSlice=1;
        end
        
        figure(1)
        clf

    	subplot(3,3,1)
    	imagesc(...
            [0,imgSize(2)].*pixelSize,...
            [0,imgSize(1)].*pixelSize,...
            squeeze(imgStack{NucSegChannel}(:,:,showZSlice)))
        axis tight equal
        xlabel('x [\mum]')
        ylabel('y [\mum]')
        title('Segmentation channel')
        hold on
        plot([x1,x2,x2,x1,x1],...
            [y1,y1,y2,y2,y1],'w-')

        
    	subplot(3,3,2)
    	imagesc(...
            [0,imgSize(2)].*pixelSize,...
            [0,imgSize(1)].*pixelSize,...
            squeeze(segImg(:,:,showZSlice)))
    	axis tight equal
        xlabel('x [\mum]')
        ylabel('y [\mum]')
        title('Preprocessed')
        hold on
        plot([x1,x2,x2,x1,x1],...
            [y1,y1,y2,y2,y1],'w-')


    	subplot(3,3,3)
    	imagesc(...
            [0,imgSize(2)].*pixelSize,...
            [0,imgSize(1)].*pixelSize,...
            squeeze(NucSegMask(:,:,showZSlice)))
    	axis tight equal
    	xlabel('x [\mum]')
        ylabel('y [\mum]')
        title('Segmentation mask')
        hold on
        plot([x1,x2,x2,x1,x1],...
            [y1,y1,y2,y2,y1],'w-')


		Cluster_subImage = imgStack{NucSegChannel}(...
			boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
			boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
			boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
		OP_subImage = OP_img(...
			boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
			boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
			boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);

        subImgSize = size(Cluster_subImage);
        if numel(subImgSize)==2
            % This is here in case we are facing 2D data
            subImgSize(3)=1;
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
		
       %Plotting of z-projected data with OP detection and segmentation
       %masks. If you want this displayed, deactivate parallel processing
       %and also uncomment the waitforbuttonpress at the end of this
       %plotting section.

       figure(1)

       subplot(3,3,4)
       
       imagesc(...
           [0,subImgSize(2)].*pixelSize,...
           [0,subImgSize(1)].*pixelSize,...
           squeeze(max(NucMask_subImage,[],3)))
       axis tight equal
       xlabel('x [\mum]')
       ylabel('y [\mum]')
       title('Nuclear segmentation')

       subplot(3,3,5)
       imagesc(...
           [0,subImgSize(2)].*pixelSize,...
           [0,subImgSize(1)].*pixelSize,...
           squeeze(max(Cluster_subImage,[],3)))
       axis tight equal
       xlabel('x [\mum]')
       ylabel('y [\mum]')
       title('Pol II Ser5P cluster')

       subplot(3,3,6)
       imagesc(...
           [0,subImgSize(2)].*pixelSize,...
           [0,subImgSize(1)].*pixelSize,...
           squeeze(max(OP_subImage,[],3)))
       axis tight equal
       xlabel('x [\mum]')
       ylabel('y [\mum]')
       title(sprintf('Oligopaint (%s)',...
           condNames{ff}))

       subplot(3,3,8)
       imagesc(...
           [0,subImgSize(2)].*pixelSize,...
           [0,subImgSize(1)].*pixelSize,...
           squeeze(max(Cluster_mask,[],3)))
       xlabel('x [\mum]')
       ylabel('y [\mum]')
       axis tight equal

       subplot(3,3,9)
       imagesc(...
           [0,subImgSize(2)].*pixelSize,...
           [0,subImgSize(1)].*pixelSize,...
           squeeze(max(OP_mask,[],3)))
       xlabel('x [\mum]')
       ylabel('y [\mum]')
       axis tight equal
       title(sprintf('Oligopaint (%d detected)',...
           OP_comps.NumObjects))


       if OP_comps.NumObjects>0

           markerSize = 12;

           OP_centroids = ...
               OP_props.Centroid.*[pixelSize,pixelSize,zStepSize];
           
           subplot(3,3,5)
           hold on
           plot(OP_centroids(:,1),OP_centroids(:,2),'k+',...
               'MarkerSize',markerSize)
           plot(OP_centroids(:,1),OP_centroids(:,2),'ws',...
               'MarkerSize',markerSize)
           hold off

           subplot(3,3,8)
           hold on
           plot(OP_centroids(:,1),OP_centroids(:,2),'k+',...
               'MarkerSize',markerSize)
           plot(OP_centroids(:,1),OP_centroids(:,2),'ws',...
               'MarkerSize',markerSize)
           hold off
                  
           subplot(3,3,6)
           hold on
           plot(OP_centroids(:,1),OP_centroids(:,2),'k+',...
               'MarkerSize',markerSize)
           plot(OP_centroids(:,1),OP_centroids(:,2),'ws',...
               'MarkerSize',markerSize)
           hold off

           subplot(3,3,9)
           hold on
           plot(OP_centroids(:,1),OP_centroids(:,2),'k+',...
               'MarkerSize',markerSize)
           plot(OP_centroids(:,1),OP_centroids(:,2),'ws',...
               'MarkerSize',markerSize)
           hold off
       
       end

       %waitforbuttonpress % uncomment this if you want image outputs

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

    figure(2)
    clf

    subplot(1,3,1)
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
   




    subplot(1,3,2)
    % Cross talk assessment within OP foci
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

    subplot(1,3,3)
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


    if numNuclei_vec(ff)>0 && poolsize == 0
        %waitforbuttonpress
    end

    if fractionAboveFour > 0.15
        validFileFlags(ff) = false;
    else
        validFileFlags(ff) = true;
    end
	
end

%% Report conditions from which data was rejected

invalidInds = find(~validFileFlags);
numInvalid = numel(invalidInds);

if numInvalid>0

    disp('There were files rejected during oligopaint quality control:')

    for ff = 1:numInvalid

        disp(' ')
        disp('Condition:')
        disp(condNames{invalidInds(ff)})
        disp('File name:')
        disp(listing(ff).name)

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
	
    validCondInds = ...
        condInds==uniqueCondInds(cc) ...
        & validFileFlags;

	sortedCondNames{cc} = ...
		condNames(validCondInds);
	sortedCondNames{cc} = sortedCondNames{cc}{1};
	
	sortedNumNuclei(cc) = ...
		sum(numNuclei_vec(validCondInds));
	sortedNumFiles(cc) = sum(validCondInds);
	
	Cluster_dists = vertcat(Cluster_distCell{validCondInds});
	Cluster_vols = vertcat(Cluster_volCell{validCondInds});
	Cluster_sols = vertcat(Cluster_solCell{validCondInds});
	Cluster_elos = vertcat(Cluster_eloCell{validCondInds});
	Cluster_masks = vertcat(Cluster_maskCell{validCondInds});
	Cluster_ints = Cluster_intCell(validCondInds);
	OP_ints = OP_intCell(validCondInds);

	sortedDistCell{cc} = Cluster_dists;
	sortedVolCell{cc} = Cluster_vols;
	sortedSolCell{cc} = Cluster_sols;
	sortedEloCell{cc} = Cluster_elos;
	sortedMaskCell{cc} = Cluster_masks;

	pixelSizes = pixelSize_array(validCondInds);
	zStepSizes = zStep_array(validCondInds);
	
	sortedPixelSize(cc) = pixelSizes(1);
	sortedZStep(cc) = zStepSizes(1);

	for qq = 1:numQuantChannels
		
		sortedOPIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)OP_intCell{ind}{qq},....
			find(validCondInds),...
			'UniformOutput',false));
		sortedOPIntCell{qq}{cc} = vertcat(sortedOPIntCell{qq}{cc}{:});
		sortedIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)Cluster_intCell{ind}{qq},....
			find(validCondInds),...
			'UniformOutput',false));
		sortedIntCell{qq}{cc} = vertcat(sortedIntCell{qq}{cc}{:});

	end
	
end

%% Saving of results

save('ConditionSortedResults')

%% Overview plots for all analyzed conditions

figure(4)
clf

prct_dist = cellfun(@(xx)prctile(xx,10),sortedDistCell); % in micrometers
prct_Ser2P = cellfun(@(xx)prctile(xx,90),sortedIntCell{1});

contact_freq = cellfun(@(xx)mean(xx<=0.50),sortedDistCell); % in micrometers);
avg_dist = cellfun(@mean,sortedDistCell); % in micrometers
avg_Ser5P = cellfun(@median,sortedIntCell{2});
avg_Ser2P = cellfun(@mean,sortedIntCell{1});

scatter(avg_Ser5P,avg_Ser2P,100,prct_dist,'filled')
xlabel('Mean Pol II Ser5P')
ylabel('Mean Pol II Ser2P')

colormap(parula)
colorbar
set(gca,'Box','on')
title('Mean gene-cluster distance [\mum]','FontWeight','normal')

