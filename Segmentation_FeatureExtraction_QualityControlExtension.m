clear all

visualizationFlag = false;

if visualizationFlag

    % Script modified to only extract for a single gene and stage, but with
    % image data to be plotted as examples

    sourceDirectory = './ExtractedStacks_Stages/Oblong/Cond_2/';
    centralSliceExtension = 43;40; % pixels from centroid

else

    % Script configuration to process all image data

    sourceDirectory = './ExtractedStacks_Stages/Oblong/**/';
    %Actin overexpression data
    sourceDirectory = '/Volumes/Seagate Bas/ActinPerturbation_OPExtractedData/ExtractedStacks/**/';
    centralSliceExtension = 0; % pixels from centroid

end

% Channels for segmentation
Nuc_SegChannel = 2; % Channel used to detect nuclei
S5P_SegChannel = 2; % Channel used to detect Pol II S5P clusters
S2P_SegChannel = 1; % Channel used to detect Pol II S2P clusters
OP_SegChannel = 3; % Channel used to detect oligopaint spots

% Save images of the clusters
ImgSquareExtension = 0; % pixels for cut-out image extension, set 0 for no images
% Which image channels to store in example images
storeImgChannels = [];
numStoreChannels = numel(storeImgChannels);

% The channels must be assigned so that the they represent the following
% labels in the biological sample:
% Channel 1: Elongating Pol II, Pol II Ser2Phos
% Channel 2: Recruited Pol II, Pol II Ser5Phos
% Channel 3: Oligopaint DNA-FISH labeling of the gene of interest
quantChannels = [1,2,3];
quantBlurSigma = [0,0,0.07];

nuc_segBlurSigma_small = 1.0; % in microns
nuc_segBlurSigma_large = 10; % in microns
nuc_segErosion = 0.5; % range of erosion (in microns) to avoid margin effects
% Use topological operation to fill holes in the nuclei segmentation masks?
% Default: 3D hole-filing, set flag to value 1
% 2D hole-filling, usefil if the stack cuts nuclei on top or bottom, so
% that 3D hole-filling does not work, set flag value to 2
% To deactivate hole-filling, set flag to any other number
fillHolesFlag = 1;

% Minimum volume of nuclei, typical ranges for a full nucieus 10-100 cubic
% microns of volume, so set a cut-off oof 10 or 30 or so
Nuc_min_vol = 40; % cubic microns
Nuc_min_sol = 0.7; % to ensure round nuclei
Nuc_min_CoV = 0.0; % to ensure transcriptionally active foci

% Inner and outer extension of a nuclear masks to cover cytoplasm
cytoMask_extension = 1.5; % in microns
cytoMask_distance = 1.0; % in microns

S5P_segBlurSigma_small = 0.0001; % in microns
S5P_segBlurSigma_large = 0.5; % in microns
S5P_seg_numStdDev = 3.0;
% Cluster connection range:
S5P_DBSCAN_epsilon = 0.5; % in microns, choose 0 for no clustering

S2P_segBlurSigma_small = 0.03; % in microns
S2P_segBlurSigma_large = 0.1; % in microns
S2P_seg_numStdDev = 2.25; % number of standard deviations in robust threshold
% Cluster connection range:
S2P_DBSCAN_epsilon = 0.0; % in microns, choose 0 for no clustering

OP_segBlurSigma_small = 0.07; % in microns
OP_segBlurSigma_large = 0.3; % in microns
OP_seg_numStdDev = 5; % number of standard deviations in robust threshold
% Cluster connection range:
OP_DBSCAN_epsilon = 0.2; % in microns, choose 0 for no clustering


% Minimum volumes for objects inside the nuclei
S5P_minVol = 0.005; % cubic microns
S2P_minVol = 0.005; % cubic microns
OP_minVol = 0.05; % cubic microns

paired_cluster_minVol = 0.08; % to only include large clusters for paired objects

dist_threshold = 0.5; % contact distance in micrometers


% ------end of analysis parameters

% Get number of parallel workers

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end

numQuantChannels = numel(quantChannels);

% --- Get all files and that should be analyzed
listing = rdir([sourceDirectory,'*Image*.mat']);
numFiles = numel(listing);
validFileFlag = false(1,numFiles);

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
% ---


% --- preallocate variables to store analysis results from the files
% Variables to store properties of nuclei
numNuclei_vec = zeros(1,numFiles);
nuc_volCell = cell(1,numFiles);
nuc_intCell = cell(1,numFiles);
cyto_intCell = cell(1,numFiles);
nuc_stdCell = cell(1,numFiles);
nuc_medianVolCell = cell(1,numFiles);
perNuc_countCell = cell(1,numFiles);
perNuc_volCell = cell(1,numFiles);

% Variables to store the voxel size
voxelSize_array_XY = zeros(1,numFiles);
voxelSize_array_Z = zeros(1,numFiles);

% Variables to store properties of objects inside nuclei
S5P_volCell = cell(1,numFiles);
S5P_solCell = cell(1,numFiles);
S5P_eloCell = cell(1,numFiles);
S5P_intCell = cell(1,numFiles);
S5P_centCell = cell(1,numFiles);
S5P_imgCell = cell(1,numFiles);
S5P_nucIntCell = cell(1,numFiles);
S5P_nucVolCell = cell(1,numFiles);
S5P_nucClustVolCell = cell(1,numFiles);

S2P_volCell = cell(1,numFiles);
S2P_solCell = cell(1,numFiles);
S2P_eloCell = cell(1,numFiles);
S2P_intCell = cell(1,numFiles);
S2P_centCell = cell(1,numFiles);
S2P_imgCell = cell(1,numFiles);
S2P_nucIntCell = cell(1,numFiles);
S2P_nucVolCell = cell(1,numFiles);
S2P_nucClustVolCell = cell(1,numFiles);

OP_volCell = cell(1,numFiles);
OP_solCell = cell(1,numFiles);
OP_eloCell = cell(1,numFiles);
OP_intCell = cell(1,numFiles);
OP_centCell = cell(1,numFiles);
OP_imgCell = cell(1,numFiles);
OP_nucIntCell = cell(1,numFiles);
OP_nucVolCell = cell(1,numFiles);
OP_nucClustVolCell = cell(1,numFiles);


% Variables to store properties of oligopaint-cluster pairs
LargeS5P_NNDistCell = cell(1,numFiles);

paired_Cluster_distCell = cell(1,numFiles);
paired_Cluster_volCell = cell(1,numFiles);
paired_Cluster_solCell = cell(1,numFiles);
paired_Cluster_eloCell = cell(1,numFiles);
paired_Cropped_ImgCell = cell(1,numFiles);
paired_Cluster_intCell = cell(1,numFiles);
paired_OP_intCell = cell(1,numFiles);
paired_OP_displCell = cell(1,numFiles);
paired_nucVolCell = cell(1,numFiles);
paired_nucClustVolumeCell = cell(1,numFiles);
paired_nucNNDistCell = cell(1,numFiles);
paired_nucIntCell = cell(1,numFiles);

% ---


% --- Analyze image stacks one by one

% If you want to switch between parallel and single thread processing, the
% easiest is to comment out the for or the parfor statement below. To not
% overload common computers, control the number of threads that are used
% with the parpool(n_threads) command before you run this script. If you
% run out of memory, reduce the number of threads in the parpool() command.

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

    voxelSize_array_XY(ff) = pixelSize;
    voxelSize_array_Z(ff) = zStepSize;

    % Nuclei segmentation
    segImg = imgStack{Nuc_SegChannel};
    if nuc_segBlurSigma_small>0
        segImg = ...
            + imgaussfilt(segImg,nuc_segBlurSigma_small./pixelSize) ...
            - imgaussfilt(segImg,nuc_segBlurSigma_large./pixelSize);
    else
        segImg = ...
            + segImg ...
            - imgaussfilt(segImg,nuc_segBlurSigma_large./pixelSize);
    end

    [bin_counts,bin_centers] = hist(segImg(:),1000);
	[nuc_seg_thresh,~] = otsuLimit(bin_centers,bin_counts,[0,Inf]);
	NucSegMask = segImg>1.0.*nuc_seg_thresh;
    if fillHolesFlag == 1
        % 3D hole-filling, default
        NucSegMask = imfill(NucSegMask,18,'holes');
    elseif fillHolesFlag == 2
        % 2D hole-filling, useful if the stack cuts most nuclei on top or
        % bottom
        NucSegMask = imfill(NucSegMask,8,'holes');
    end
    if nuc_segErosion>0
        se = strel('disk', ...
            round(nuc_segErosion./pixelSize)); % Two-dimensional erosion disk
        NucSegMask = imerode(NucSegMask,se);
    end

    subplot(1,3,1)
	imagesc(squeeze(imgStack{Nuc_SegChannel}(:,:,ceil(imgSize(3)./2))))
	axis tight equal
	
	subplot(1,3,2)
	imagesc(squeeze(segImg(:,:,ceil(imgSize(3)./2))))
	axis tight equal

	subplot(1,3,3)
	imagesc(squeeze(NucSegMask(:,:,ceil(imgSize(3)./2))))
	axis tight equal
	
    % Uncomment the following two lines, and remove the par in parfor above
    % if you want to check the extracted images one by one
    %   fprintf('File name: %s\n',thisFilePath)
    % 	waitforbuttonpress
    

    % --- Connected component segmentation of nuclei
	comps = bwconncomp(NucSegMask,18);
	numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList);
	minPixels = Nuc_min_vol./(pixelSize.^2)./zStepSize;
	comps.NumObjects = sum(numPxls>=minPixels);
	comps.PixelIdxList = comps.PixelIdxList(numPxls>=minPixels);
	numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList);
	
	props = regionprops3(comps,imgStack{Nuc_SegChannel},...
		'Solidity','VoxelValues');
	
    if comps.NumObjects>0
        Solidity_array = [props.Solidity];
        CoV_array = ...
            cellfun(@(vals)std(vals(:))./mean(vals(:)),...
            props.VoxelValues);
        inclNucInds = Solidity_array>=Nuc_min_sol ...
            & CoV_array>=Nuc_min_CoV;
        comps.NumObjects = sum(Solidity_array>=Nuc_min_sol);
        comps.PixelIdxList = comps.PixelIdxList(Solidity_array>=Nuc_min_sol);
        numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList);
    end
    
	numNuclei = comps.NumObjects;
	numNuclei_vec(ff) = numNuclei;

    if comps.NumObjects>0
        
        validFileFlag(ff) = true;
        
        nuc_intCell{ff} = cell(1,numQuantChannels);
        cyto_intCell{ff} = cell(1,numQuantChannels);
        nuc_stdCell{ff} = cell(1,numQuantChannels);
        
        for qq = 1:numQuantChannels
            quantImg = imgStack{quantChannels(qq)};
            quantProps = regionprops3(comps,quantImg,...
                'MeanIntensity','VoxelIdxList','VoxelValues');
            nuc_intCell{ff}{qq} = [quantProps.MeanIntensity];
            cyto_intCell{ff}{qq} = zeros(numNuclei,1);
            nuc_stdCell{ff}{qq} = cellfun(...
                @std,quantProps.VoxelValues);
            
            for nn = 1:numNuclei
                cytoMask = false(size(quantImg));
                cytoMask(quantProps.VoxelIdxList{nn}) = true;
                coreMask = cytoMask;
                se = strel('disk',round(cytoMask_extension./pixelSize));
                cytoMask = imdilate(cytoMask,se);
                se = strel('disk',round(cytoMask_distance./pixelSize));
                coreMask = imdilate(coreMask,se);
                cytoMask(coreMask) = false;
                cytoMask = cytoMask & ~NucSegMask;
                cyto_intCell{ff}{qq}(nn) = mean(quantImg(cytoMask));
            end
        end
        
        props = regionprops3(comps,imgStack{Nuc_SegChannel},...
            'Volume','VoxelValues','Solidity','VoxelIdxList',...
            'BoundingBox');
        
        Volume_array = [props.Volume].*pixelSize.^2.*zStepSize;
        Intensity_array = cellfun(@(vals)median(vals),props.VoxelValues);
        Solidity_array = [props.Solidity];

        nuc_volCell{ff} = Volume_array;

 
        
        % --- For each nucleus, get objects from the different channels
        
        nuc_medianVolCell{ff} = cell(1,3);
        perNuc_countCell{ff} = cell(1,3);
        perNuc_volCell{ff} = cell(1,3);
        for qq = 1:3
            nuc_medianVolCell{ff}{qq} = zeros(numNuclei,1);
            perNuc_countCell{ff}{qq} = zeros(numNuclei,1);
            perNuc_volCell{ff}{qq} = cell(numNuclei,1);
        end
        
        S5P_volume = cell(1,numNuclei);
        S5P_solidity = cell(1,numNuclei);
        S5P_elongation = cell(1,numNuclei);
        S5P_centralSlices_store = cell(1,numNuclei);
        S5P_intensity = cell(numQuantChannels,numNuclei);
        S5P_cent_store = cell(1,numNuclei);
        S5P_nucIntensity = cell(numQuantChannels,numNuclei);
        S5P_nucVolume = cell(1,numNuclei);
        S5P_nucClustVolume = cell(1,numNuclei);
        
        S2P_volume = cell(1,numNuclei);
        S2P_solidity = cell(1,numNuclei);
        S2P_elongation = cell(1,numNuclei);
        S2P_centralSlices_store = cell(1,numNuclei);
        S2P_intensity = cell(numQuantChannels,numNuclei);
        S2P_cent_store = cell(1,numNuclei);
        S2P_nucIntensity = cell(numQuantChannels,numNuclei);
        S2P_nucVolume = cell(1,numNuclei);
        S2P_nucClustVolume = cell(1,numNuclei);

        OP_volume = cell(1,numNuclei);
        OP_solidity = cell(1,numNuclei);
        OP_elongation = cell(1,numNuclei);
        OP_centralSlices_store = cell(1,numNuclei);
        OP_intensity = cell(numQuantChannels,numNuclei);
        OP_cent_store = cell(1,numNuclei);
        OP_nucIntensity = cell(numQuantChannels,numNuclei);
        OP_nucVolume = cell(1,numNuclei);
        OP_nucClustVolume = cell(1,numNuclei);


        % --- paired object results
        LargeS5P_NNDist = cell(1,numNuclei);

        paired_Cluster_OP_dist = cell(1,numNuclei);
        paired_Cluster_volume = cell(1,numNuclei);
        paired_Cluster_solidity = cell(1,numNuclei);
        paired_Cluster_elongation = cell(1,numNuclei);
        paired_Cropped_Img_store = cell(1,numNuclei);
        paired_Cluster_intensity = cell(numQuantChannels,numNuclei);
        paired_OP_intensity = cell(numQuantChannels,numNuclei);
        paired_OP_Cluster_disp = cell(1,numNuclei);
        paired_nucVolume = cell(1,numNuclei);
        paired_nucClustVolume = cell(1,numNuclei);
        paired_nucClustNNDist = cell(1,numNuclei);
        paired_nucIntensity = cell(numQuantChannels,numNuclei);



        
        for nn = 1:numNuclei
            
            boxArray = props.BoundingBox(nn,:);
            
            S5P_subImage = imgStack{S5P_SegChannel}(...
                boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
                boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
                boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
            
            S2P_subImage = imgStack{S2P_SegChannel}(...
                boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
                boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
                boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);

            OP_subImage = imgStack{OP_SegChannel}(...
                boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
                boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
                boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);

            
            S5P_subImage = ...
                + imgaussfilt(S5P_subImage,S5P_segBlurSigma_small./pixelSize) ...
                - imgaussfilt(S5P_subImage,S5P_segBlurSigma_large./pixelSize);
            S2P_subImage = ...
                + imgaussfilt(S2P_subImage,S2P_segBlurSigma_small./pixelSize) ...
                - imgaussfilt(S2P_subImage,S2P_segBlurSigma_large./pixelSize);
            OP_subImage = ...
                + imgaussfilt(OP_subImage,OP_segBlurSigma_small./pixelSize) ...
                - imgaussfilt(OP_subImage,OP_segBlurSigma_large./pixelSize);
            

            NucMask = false(imgSize);
            NucMask(props.VoxelIdxList{nn}) = true;
            NucMask_subImage = NucMask(...
                boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
                boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
                boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
            
            
            seg_intensities = S5P_subImage(NucMask_subImage);
            seg_mean = mean(seg_intensities);
            seg_std = std(seg_intensities);
            S5P_mask = (S5P_subImage.*NucMask_subImage)...
                >(seg_mean+S5P_seg_numStdDev.*seg_std);
            
            seg_intensities = S2P_subImage(NucMask_subImage);
            seg_mean = mean(seg_intensities);
            seg_std = std(seg_intensities);
            S2P_mask = (S2P_subImage.*NucMask_subImage)...
                >(seg_mean+S2P_seg_numStdDev.*seg_std);

            seg_intensities = OP_subImage(NucMask_subImage);
            seg_mean = mean(seg_intensities);
            seg_std = std(seg_intensities);
            OP_mask = (OP_subImage.*NucMask_subImage)...
                >(seg_mean+OP_seg_numStdDev.*seg_std);
            
            subImgSize = size(S5P_subImage);
            if numel(subImgSize)==2
                subImgSize(3)=1;
            end
            
            % For storage of example images
            if numStoreChannels > 0
                store_subImages = cell(1,numStoreChannels);
                for color = 1:numStoreChannels
                    store_subImages{color} = imgStack{storeImgChannels(color)}(...
                        boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
                        boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
                        boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
                    nuc_intensities = ...
                        store_subImages{color}(NucMask_subImage);
                    store_subImages{color} = ...
                        store_subImages{color}./median(nuc_intensities(:));
                end
            end
            
            S5P_comps = bwconncomp(S5P_mask,18);
            S5P_numPxls = cellfun(@(elmt)numel(elmt),S5P_comps.PixelIdxList);
            minPixels = S5P_minVol./(pixelSize.^2)./zStepSize;
            S5P_comps.NumObjects = sum(S5P_numPxls>minPixels);
            S5P_comps.PixelIdxList = S5P_comps.PixelIdxList(S5P_numPxls>minPixels);
            S5P_numPxls = cellfun(@(elmt)numel(elmt),S5P_comps.PixelIdxList);
            
            S2P_comps = bwconncomp(S2P_mask,18);
            S2P_numPxls = cellfun(@(elmt)numel(elmt),S2P_comps.PixelIdxList);
            minPixels = S2P_minVol./(pixelSize.^2)./zStepSize;
            S2P_comps.NumObjects = sum(S2P_numPxls>minPixels);
            S2P_comps.PixelIdxList = S2P_comps.PixelIdxList(S2P_numPxls>minPixels);
            S2P_numPxls = cellfun(@(elmt)numel(elmt),S2P_comps.PixelIdxList);

            OP_comps = bwconncomp(OP_mask,18);
            OP_numPxls = cellfun(@(elmt)numel(elmt),OP_comps.PixelIdxList);
            minPixels = OP_minVol./(pixelSize.^2)./zStepSize;
            OP_comps.NumObjects = sum(OP_numPxls>minPixels);
            OP_comps.PixelIdxList = OP_comps.PixelIdxList(OP_numPxls>minPixels);
            OP_numPxls = cellfun(@(elmt)numel(elmt),OP_comps.PixelIdxList);
            
            if S5P_comps.NumObjects>0
                
                % DBSCAN clustering of labeled regions
                if S5P_DBSCAN_epsilon > 0
                    S5P_props = regionprops3(S5P_comps,S5P_subImage,...
                        'Centroid');
                    centroid_coords = ...
                        S5P_props.Centroid.*[pixelSize,pixelSize,zStepSize];
                    dbscan_inds = ...
                        dbscan(centroid_coords,S5P_DBSCAN_epsilon,1);
                    
                    unique_inds = unique(dbscan_inds);
                    num_inds = numel(unique_inds);
                    updated_comps = S5P_comps;
                    updated_comps.NumObjects = num_inds;
                    updated_comps.PixelIdxList = cell(1,num_inds);
                    for ii = 1:num_inds
                        updated_comps.PixelIdxList{ii} = ...
                            sort(vertcat(S5P_comps.PixelIdxList{...
                            dbscan_inds==unique_inds(ii)} ...
                            ));
                    end
                    S5P_comps = updated_comps;
                end
                
                LL = labelmatrix(S5P_comps);
                
                subplot(2,2,1)
                centerPlaneInd = round(boxArray(6).*0.5);
                % 			imagesc(squeeze(max(S5P_subImage,[],3)))
                imagesc(squeeze(S5P_subImage(:,:,centerPlaneInd)))
                axis tight equal
                set(gca,'Colormap',gray)
                
                subplot(2,2,2)
                % 			imagesc(squeeze(max(S2P_subImage,[],3)))
                imagesc(squeeze(S2P_subImage(:,:,centerPlaneInd)))
                axis tight equal
                set(gca,'Colormap',gray)
                
                subplot(2,2,3)
                % 			imagesc(squeeze(max(S5P_mask,[],3)))
                imagesc(squeeze(LL(:,:,centerPlaneInd)))
                axis tight equal
                set(gca,'Colormap',lines)
                
                subplot(2,2,4)
                % 			imagesc(squeeze(max(S2P_mask,[],3)))
                imagesc(squeeze(S2P_mask(:,:,centerPlaneInd)))
                axis tight equal
                set(gca,'Colormap',gray)
                
                % Uncomment this waitforbuttonpress command to see the
                % segmentation results for the two types of foci
                % waitforbuttonpress
                
                
                
                
                S5P_props = regionprops3(S5P_comps,S5P_subImage,...
                    'Volume','Solidity',...
                    'Centroid','Image','BoundingBox');
                
                S5P_Volume_array = ...
                    [S5P_props.Volume].*pixelSize.^2.*zStepSize;
                S5P_Solidity_array = [S5P_props.Solidity];
                S5P_Centroid_array = ...
                    S5P_props.Centroid.*[pixelSize,pixelSize,zStepSize];
                                
                perNuc_countCell{ff}{1}(nn) = numel(S5P_Volume_array);
                perNuc_volCell{ff}{1}{nn} = S5P_Volume_array;
                nuc_medianVolCell{ff}{1}(nn) = median(S5P_Volume_array);
                
                
                % --- get cluster central plane and elongation in-plane
                S5P_Elongation_array = ...
                    zeros(size(S5P_Solidity_array));
                S5P_Slices_cell = ...
                    cell(size(S5P_Solidity_array));
                
                for object_nn = 1:numel(S5P_Volume_array)
                    
                    boundingBox = S5P_props.BoundingBox(object_nn,:);
                    thisImage = squeeze(S5P_subImage(...
                        boundingBox(2)+0.5:boundingBox(2)+boundingBox(5)-0.5,...
                        boundingBox(1)+0.5:boundingBox(1)+boundingBox(4)-0.5,...
                        ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));
                    thisImage = thisImage-min(thisImage(:));
                    thisImage = thisImage./max(thisImage(:));
                    
                    thisMask = S5P_props.Image(object_nn);
                    if iscell(thisMask)
                        % There is some weird inconsistent behavior here,
                        % sometimes .Image(object_nn) results in a cell
                        % output, sometimes in a matrix output. It is not
                        % clear if there is a system to it.
                        thisMask = thisMask{1};
                    else
                        thisMask = S5P_props.Image;
                    end
                    centerInd = ceil(size(thisMask,3)./2);
                    thisMask = squeeze(thisMask(:,:,centerInd));
                    thisImage((bwperim(thisMask))) = 0;
                    
                    S5P_Slices_cell{object_nn} = cell(1,numStoreChannels+2);
                    centroid_1 = round(S5P_props.Centroid(object_nn,2));
                    centroid_2 = round(S5P_props.Centroid(object_nn,1));
                    center_z = round(S5P_props.Centroid(object_nn,3));
                    img_limits = [...
                        centroid_1-ImgSquareExtension,...
                        centroid_1+ImgSquareExtension,...
                        centroid_2-ImgSquareExtension,...
                        centroid_2+ImgSquareExtension];
                    img_limits = [...
                        max(1,img_limits(1)),...
                        min(subImgSize(1),img_limits(2)),...
                        max(1,img_limits(3)),...
                        min(subImgSize(2),img_limits(4))];
                    
                    for color = 1:numStoreChannels
                        % 						S5P_Slices_cell{object_nn}{color} = ...
                        % 							squeeze(store_subImages{color}(:,:,...
                        % 							ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));
                        
                        S5P_Slices_cell{object_nn}{color} = ...
                            squeeze(store_subImages{color}(...
                            img_limits(1):img_limits(2),...
                            img_limits(3):img_limits(4),...
                            center_z));
                        
                    end
                    S5P_Slices_cell{object_nn}{numStoreChannels+1} = ...
                        squeeze(LL(...
                        img_limits(1):img_limits(2),...
                        img_limits(3):img_limits(4),...
                        center_z));
                    S5P_Slices_cell{object_nn}{numStoreChannels+2} = ...
                        squeeze(S2P_mask(...
                        img_limits(1):img_limits(2),...
                        img_limits(3):img_limits(4),...
                        center_z));
                    
                    if sum(sum(uint8(thisMask)))>0
                        thisProps = regionprops(uint8(thisMask),...
                            'MajorAxisLength','MinorAxisLength');
                        S5P_Elongation_array(object_nn) = ...
                            thisProps.MajorAxisLength./thisProps.MinorAxisLength;
                    else
                        S5P_Elongation_array(object_nn) = NaN;
                    end
                end
                S5P_volume{nn} = S5P_Volume_array;
                S5P_solidity{nn} = S5P_Solidity_array;
                S5P_elongation{nn} = S5P_Elongation_array;
                S5P_cent_store{nn} = S5P_Centroid_array;
                S5P_centralSlices_store{nn} = S5P_Slices_cell;
                
                S5P_nucVolume{nn} = ...
                        ones(size(S5P_Volume_array)) ...
                        .*Volume_array(nn);
                S5P_nucClustVolume{nn} = ...
                    cell(size(S5P_Volume_array));
                for object_nn = 1:numel(S5P_Volume_array)
                    S5P_nucClustVolume{nn}{object_nn} = ...
                        S5P_Volume_array;
                end
                                                
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
                    
                    S5P_quant_props = regionprops3(...
                        S5P_comps,quant_subImage,'VoxelValues');
                    Quant_ClusterMedian = cellfun(...
                        @(vals)median(vals),S5P_quant_props.VoxelValues);
                    S5P_intensity{qq,nn} = ...
                        Quant_ClusterMedian./Quant_nucleusMedian;
                    S5P_nucIntensity{qq,nn} = ...
                        ones(size(S5P_intensity{qq,nn})) ...
                        .*nuc_intCell{ff}{qq}(nn);
                    
                end
                
            else
                
                S5P_volume{nn} = [];
                S5P_solidity{nn} = [];
                S5P_elongation{nn} = [];
                S5P_centralSlices_store{nn} = {};
                S5P_cent_store{nn} = [];
                S5P_nucVolume{nn} = [];
                S5P_nucClustVolume{nn} = {};
                
                for qq = 1:numQuantChannels
                    S5P_intensity{qq,nn} = [];
                    S5P_nucIntensity{qq,nn} = [];
                end
                
            end

            if S2P_comps.NumObjects>0

                % DBSCAN clustering of labeled regions
                if S2P_DBSCAN_epsilon > 0
                    S2P_props = regionprops3(S2P_comps,S2P_subImage,...
                        'Centroid');
                    centroid_coords = ...
                        S2P_props.Centroid.*[pixelSize,pixelSize,zStepSize];
                    dbscan_inds = ...
                        dbscan(centroid_coords,S2P_DBSCAN_epsilon,1);

                    unique_inds = unique(dbscan_inds);
                    num_inds = numel(unique_inds);
                    updated_comps = S2P_comps;
                    updated_comps.NumObjects = num_inds;
                    updated_comps.PixelIdxList = cell(1,num_inds);
                    for ii = 1:num_inds
                        updated_comps.PixelIdxList{ii} = ...
                            sort(vertcat(S2P_comps.PixelIdxList{...
                            dbscan_inds==unique_inds(ii)} ...
                            ));
                    end
                    S2P_comps = updated_comps;
                end

                LL = labelmatrix(S2P_comps);

                subplot(2,2,1)
                centerPlaneInd = round(boxArray(6).*0.5);
                % 			imagesc(squeeze(max(S2P_subImage,[],3)))
                imagesc(squeeze(S2P_subImage(:,:,centerPlaneInd)))
                axis tight equal
                set(gca,'Colormap',gray)

                subplot(2,2,2)
                % 			imagesc(squeeze(max(S2P_subImage,[],3)))
                imagesc(squeeze(S2P_subImage(:,:,centerPlaneInd)))
                axis tight equal
                set(gca,'Colormap',gray)

                subplot(2,2,3)
                % 			imagesc(squeeze(max(S2P_mask,[],3)))
                imagesc(squeeze(LL(:,:,centerPlaneInd)))
                axis tight equal
                set(gca,'Colormap',lines)

                subplot(2,2,4)
                % 			imagesc(squeeze(max(S2P_mask,[],3)))
                imagesc(squeeze(S2P_mask(:,:,centerPlaneInd)))
                axis tight equal
                set(gca,'Colormap',gray)

                % Uncomment this waitforbuttonpress command to see the
                % segmentation results for the two types of foci
                % waitforbuttonpress

                S2P_props = regionprops3(S2P_comps,S2P_subImage,...
                    'Volume','Solidity',...
                    'Centroid','Image','BoundingBox');

                S2P_Volume_array = ...
                    [S2P_props.Volume].*pixelSize.^2.*zStepSize;
                S2P_Solidity_array = [S2P_props.Solidity];
                S2P_Centroid_array = ...
                    S2P_props.Centroid.*[pixelSize,pixelSize,zStepSize];

                perNuc_countCell{ff}{2}(nn) = numel(S2P_Volume_array);

                perNuc_volCell{ff}{2}{nn} = S2P_Volume_array;

                nuc_medianVolCell{ff}{2}(nn) = median(S2P_Volume_array);

                S2P_Elongation_array = ...
                    zeros(size(S2P_Solidity_array));
                S2P_Slices_cell = ...
                    cell(size(S2P_Solidity_array));

                for object_nn = 1:numel(S2P_Volume_array)

                    boundingBox = S2P_props.BoundingBox(object_nn,:);
                    thisImage = squeeze(S2P_subImage(...
                        boundingBox(2)+0.5:boundingBox(2)+boundingBox(5)-0.5,...
                        boundingBox(1)+0.5:boundingBox(1)+boundingBox(4)-0.5,...
                        ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));
                    thisImage = thisImage-min(thisImage(:));
                    thisImage = thisImage./max(thisImage(:));

                    thisMask = S2P_props.Image(object_nn);
                    if iscell(thisMask)
                        % There is some weird inconsistent behavior here,
                        % sometimes .Image(object_nn) results in a cell
                        % output, sometimes in a matrix output. It is not
                        % clear if there is a system to it.
                        thisMask = thisMask{1};
                    else
                        thisMask = S2P_props.Image;
                    end
                    centerInd = ceil(size(thisMask,3)./2);
                    thisMask = squeeze(thisMask(:,:,centerInd));
                    thisImage((bwperim(thisMask))) = 0;

                    centroid_1 = round(S2P_props.Centroid(object_nn,2));
                    centroid_2 = round(S2P_props.Centroid(object_nn,1));
                    center_z = round(S2P_props.Centroid(object_nn,3));
                    img_limits = [...
                        centroid_1-ImgSquareExtension,...
                        centroid_1+ImgSquareExtension,...
                        centroid_2-ImgSquareExtension,...
                        centroid_2+ImgSquareExtension];
                    img_limits = [...
                        max(1,img_limits(1)),...
                        min(subImgSize(1),img_limits(2)),...
                        max(1,img_limits(3)),...
                        min(subImgSize(2),img_limits(4))];

                    for color = 1:numStoreChannels
                        S2P_Slices_cell{object_nn}{color} = ...
                            squeeze(store_subImages{color}(...
                            img_limits(1):img_limits(2),...
                            img_limits(3):img_limits(4),...
                            center_z));
                    end
                    S2P_Slices_cell{object_nn}{numStoreChannels+1} = ...
                        squeeze(S5P_mask(...
                        img_limits(1):img_limits(2),...
                        img_limits(3):img_limits(4),...
                        center_z));
                    S2P_Slices_cell{object_nn}{numStoreChannels+2} = ...
                        squeeze(S2P_mask(...
                        img_limits(1):img_limits(2),...
                        img_limits(3):img_limits(4),...
                        center_z));

                    if sum(sum(uint8(thisMask)))>0
                        thisProps = regionprops(uint8(thisMask),...
                            'MajorAxisLength','MinorAxisLength');
                        S2P_Elongation_array(object_nn) = ...
                            thisProps.MajorAxisLength./thisProps.MinorAxisLength;
                    else
                        S2P_Elongation_array(object_nn) = NaN;
                    end

                end
                S2P_volume{nn} = S2P_Volume_array;
                S2P_solidity{nn} = S2P_Solidity_array;
                S2P_elongation{nn} = S2P_Elongation_array;
                S2P_cent_store{nn} = S2P_Centroid_array;
                S2P_centralSlices_store{nn} = S2P_Slices_cell;

                S2P_nucVolume{nn} = ...
                    ones(size(S2P_Volume_array)) ...
                    .*Volume_array(nn);
                S2P_nucClustVolume{nn} = ...
                    cell(size(S2P_Volume_array));
                for object_nn = 1:numel(S2P_Volume_array)
                    S2P_nucClustVolume{nn}{object_nn} = ...
                        S2P_Volume_array;
                end


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

                    S2P_quant_props = regionprops3(...
                        S2P_comps,quant_subImage,'VoxelValues');
                    Quant_ClusterMedian = cellfun(...
                        @(vals)median(vals),S2P_quant_props.VoxelValues);
                    S2P_intensity{qq,nn} = ...
                        Quant_ClusterMedian./Quant_nucleusMedian;
                    S2P_nucIntensity{qq,nn} = ...
                        ones(size(S2P_intensity{qq,nn})) ...
                        .*nuc_intCell{ff}{qq}(nn);

                end

            else

                S2P_volume{nn} = [];
                S2P_solidity{nn} = [];
                S2P_elongation{nn} = [];
                S2P_centralSlices_store{nn} = {};
                S2P_cent_store{nn} = [];
                S2P_nucVolume{nn} = [];
                S2P_nucClustVolume{nn} = {};

                for qq = 1:numQuantChannels
                    S2P_intensity{qq,nn} = [];
                    S2P_nucIntensity{qq,nn} = [];
                end

            end


            if OP_comps.NumObjects>0

                % DBSCAN clustering of labeled regions
                if OP_DBSCAN_epsilon > 0
                    OP_props = regionprops3(OP_comps,OP_subImage,...
                        'Centroid');
                    centroid_coords = ...
                        OP_props.Centroid.*[pixelSize,pixelSize,zStepSize];
                    dbscan_inds = ...
                        dbscan(centroid_coords,OP_DBSCAN_epsilon,1);

                    unique_inds = unique(dbscan_inds);
                    num_inds = numel(unique_inds);
                    updated_comps = OP_comps;
                    updated_comps.NumObjects = num_inds;
                    updated_comps.PixelIdxList = cell(1,num_inds);
                    for ii = 1:num_inds
                        updated_comps.PixelIdxList{ii} = ...
                            sort(vertcat(OP_comps.PixelIdxList{...
                            dbscan_inds==unique_inds(ii)} ...
                            ));
                    end
                    OP_comps = updated_comps;
                end

                LL = labelmatrix(OP_comps);

                subplot(2,2,1)
                centerPlaneInd = round(boxArray(6).*0.5);
                % 			imagesc(squeeze(max(OP_subImage,[],3)))
                imagesc(squeeze(OP_subImage(:,:,centerPlaneInd)))
                axis tight equal
                set(gca,'Colormap',gray)

                subplot(2,2,2)
                % 			imagesc(squeeze(max(OP_subImage,[],3)))
                imagesc(squeeze(OP_subImage(:,:,centerPlaneInd)))
                axis tight equal
                set(gca,'Colormap',gray)

                subplot(2,2,3)
                % 			imagesc(squeeze(max(OP_mask,[],3)))
                imagesc(squeeze(LL(:,:,centerPlaneInd)))
                axis tight equal
                set(gca,'Colormap',lines)

                subplot(2,2,4)
                % 			imagesc(squeeze(max(OP_mask,[],3)))
                imagesc(squeeze(OP_mask(:,:,centerPlaneInd)))
                axis tight equal
                set(gca,'Colormap',gray)

                % Uncomment this waitforbuttonpress command to see the
                % segmentation results for the two types of foci
                % waitforbuttonpress

                OP_props = regionprops3(OP_comps,OP_subImage,...
                    'Volume','Solidity',...
                    'Centroid','Image','BoundingBox');

                OP_Volume_array = ...
                    [OP_props.Volume].*pixelSize.^2.*zStepSize;
                OP_Solidity_array = [OP_props.Solidity];
                OP_Centroid_array = ...
                    OP_props.Centroid.*[pixelSize,pixelSize,zStepSize];

                perNuc_countCell{ff}{3}(nn) = numel(OP_Volume_array);

                perNuc_volCell{ff}{3}{nn} = OP_Volume_array;

                nuc_medianVolCell{ff}{3}(nn) = median(OP_Volume_array);

                OP_Elongation_array = ...
                    zeros(size(OP_Solidity_array));
                OP_Slices_cell = ...
                    cell(size(OP_Solidity_array));

                for object_nn = 1:numel(OP_Volume_array)

                    boundingBox = OP_props.BoundingBox(object_nn,:);
                    thisImage = squeeze(OP_subImage(...
                        boundingBox(2)+0.5:boundingBox(2)+boundingBox(5)-0.5,...
                        boundingBox(1)+0.5:boundingBox(1)+boundingBox(4)-0.5,...
                        ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));
                    thisImage = thisImage-min(thisImage(:));
                    thisImage = thisImage./max(thisImage(:));

                    thisMask = OP_props.Image(object_nn);
                    if iscell(thisMask)
                        % There is some weird inconsistent behavior here,
                        % sometimes .Image(object_nn) results in a cell
                        % output, sometimes in a matrix output. It is not
                        % clear if there is a system to it.
                        thisMask = thisMask{1};
                    else
                        thisMask = OP_props.Image;
                    end
                    centerInd = ceil(size(thisMask,3)./2);
                    thisMask = squeeze(thisMask(:,:,centerInd));
                    thisImage((bwperim(thisMask))) = 0;

                    centroid_1 = round(OP_props.Centroid(object_nn,2));
                    centroid_2 = round(OP_props.Centroid(object_nn,1));
                    center_z = round(OP_props.Centroid(object_nn,3));
                    img_limits = [...
                        centroid_1-ImgSquareExtension,...
                        centroid_1+ImgSquareExtension,...
                        centroid_2-ImgSquareExtension,...
                        centroid_2+ImgSquareExtension];
                    img_limits = [...
                        max(1,img_limits(1)),...
                        min(subImgSize(1),img_limits(2)),...
                        max(1,img_limits(3)),...
                        min(subImgSize(2),img_limits(4))];

                    for color = 1:numStoreChannels
                        OP_Slices_cell{object_nn}{color} = ...
                            squeeze(store_subImages{color}(...
                            img_limits(1):img_limits(2),...
                            img_limits(3):img_limits(4),...
                            center_z));
                    end
                    OP_Slices_cell{object_nn}{numStoreChannels+1} = ...
                        squeeze(S5P_mask(...
                        img_limits(1):img_limits(2),...
                        img_limits(3):img_limits(4),...
                        center_z));
                    OP_Slices_cell{object_nn}{numStoreChannels+2} = ...
                        squeeze(S2P_mask(...
                        img_limits(1):img_limits(2),...
                        img_limits(3):img_limits(4),...
                        center_z));

                    if sum(sum(uint8(thisMask)))>0
                        thisProps = regionprops(uint8(thisMask),...
                            'MajorAxisLength','MinorAxisLength');
                        OP_Elongation_array(object_nn) = ...
                            thisProps.MajorAxisLength./thisProps.MinorAxisLength;
                    else
                        OP_Elongation_array(object_nn) = NaN;
                    end

                end
                OP_volume{nn} = OP_Volume_array;
                OP_solidity{nn} = OP_Solidity_array;
                OP_elongation{nn} = OP_Elongation_array;
                OP_cent_store{nn} = OP_Centroid_array;
                OP_centralSlices_store{nn} = OP_Slices_cell;

                OP_nucVolume{nn} = ...
                    ones(size(OP_Volume_array)) ...
                    .*Volume_array(nn);
                OP_nucClustVolume{nn} = ...
                    cell(size(OP_Volume_array));
                for object_nn = 1:numel(OP_Volume_array)
                    OP_nucClustVolume{nn}{object_nn} = ...
                        OP_Volume_array;
                end


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

                    OP_quant_props = regionprops3(...
                        OP_comps,quant_subImage,'VoxelValues');
                    Quant_ClusterMedian = cellfun(...
                        @(vals)median(vals),OP_quant_props.VoxelValues);
                    OP_intensity{qq,nn} = ...
                        Quant_ClusterMedian./Quant_nucleusMedian;
                    OP_nucIntensity{qq,nn} = ...
                        ones(size(OP_intensity{qq,nn})) ...
                        .*nuc_intCell{ff}{qq}(nn);

                end





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

                    OP_quant_props = regionprops3(...
                        OP_comps,quant_subImage,'VoxelValues');
                    Quant_ClusterMedian = cellfun(...
                        @(vals)median(vals),OP_quant_props.VoxelValues);
                    OP_intensity{qq,nn} = ...
                        Quant_ClusterMedian./Quant_nucleusMedian;
                    OP_nucIntensity{qq,nn} = ...
                        ones(size(OP_intensity{qq,nn})) ...
                        .*nuc_intCell{ff}{qq}(nn);

                end

            else

                OP_volume{nn} = [];
                OP_solidity{nn} = [];
                OP_elongation{nn} = [];
                OP_centralSlices_store{nn} = {};
                OP_cent_store{nn} = [];
                OP_nucVolume{nn} = [];
                OP_nucClustVolume{nn} = {};

                for qq = 1:numQuantChannels
                    OP_intensity{qq,nn} = [];
                    OP_nucIntensity{qq,nn} = [];
                end

            end


            % --- From here on pairing of clusters and oligopaint spots

            Cluster_comps = S5P_comps;

            if Cluster_comps.NumObjects>0

                Cluster_numPxls = cellfun(@(elmt)numel(elmt),...
                    Cluster_comps.PixelIdxList);
                minPixels = paired_cluster_minVol./(pixelSize.^2)./zStepSize;
                Cluster_comps.NumObjects = sum(Cluster_numPxls>minPixels);
                Cluster_comps.PixelIdxList = ...
                    Cluster_comps.PixelIdxList(Cluster_numPxls>minPixels);
                Cluster_numPxls = cellfun(@(elmt)numel(elmt),...
                    Cluster_comps.PixelIdxList);

            end

            if Cluster_comps.NumObjects>1
                Cluster_props = regionprops3(Cluster_comps,S5P_subImage,...
                    'Centroid');
                centroids = ...
                    Cluster_props.Centroid.*[pixelSize,pixelSize,zStepSize];
                % Gets the nearest neighbor distance between large Pol II
                % clusters
                distMatrix = squareform(pdist(centroids));
                distMatrix = sort(distMatrix,2);
                LargeS5P_NNDist{nn} = distMatrix(:,2);
            else
                LargeS5P_NNDist{nn} = [];
            end

            if OP_comps.NumObjects>0 && Cluster_comps.NumObjects>0

                OP_props = regionprops3(OP_comps,S5P_subImage,...
                    'Centroid');

                Cluster_props = regionprops3(Cluster_comps,S5P_subImage,...
                    'Volume','MeanIntensity','Solidity',...
                    'Centroid','Image','BoundingBox');

                pwDist = pdist2(...
                    OP_props.Centroid.*[pixelSize,pixelSize,zStepSize],...
                    Cluster_props.Centroid.*[pixelSize,pixelSize,zStepSize]);
                [minDist,minInds] = min(pwDist,[],2);

                Cluster_OP_minDist = minDist;
                Cluster_OP_minInd = minInds;

                OP_centroid_array = ...
                    OP_props.Centroid.*[pixelSize,pixelSize,zStepSize];
                Cluster_centroid_array = ...
                    Cluster_props.Centroid.*[pixelSize,pixelSize,zStepSize];
                Cluster_centroid_array = Cluster_centroid_array(Cluster_OP_minInd,:);

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
                paired_Cluster_OP_dist{nn} = Cluster_OP_minDist;
                paired_OP_Cluster_disp{nn} = ...
                    OP_centroid_array-Cluster_centroid_array;
                paired_Cluster_volume{nn} = Cluster_Volume_array;
                paired_Cluster_solidity{nn} = Cluster_Solidity_array;
                paired_Cluster_elongation{nn} = Cluster_Elongation_array;

                % --- quantification for all target channels
                Cropped_Img_cell = ...
                    cell(size(Cluster_Solidity_array));
                for cl = 1:numel(Cluster_Volume_array)
                    Cropped_Img_cell{cl} = cell(1,numQuantChannels);
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
                    paired_Cluster_intensity{qq,nn} = ...
                        Quant_ClusterMedian./Quant_nucleusMedian;

                    OP_quant_props = regionprops3(...
                        OP_comps,quant_subImage,'VoxelValues');
                    Quant_OPMedian = cellfun(...
                        @(vals)median(vals),OP_quant_props.VoxelValues);
                    paired_OP_intensity{qq,nn} = ...
                        Quant_OPMedian./Quant_nucleusMedian;

                    if centralSliceExtension==0

                        for cl = 1:numel(Cluster_Volume_array)
                            Cropped_Img_cell{cl}{qq} = [];
                        end

                    else

                        for cl = 1:numel(Cluster_Volume_array)

                            centroid_1 = round(Cluster_props.Centroid(cl,2));
                            centroid_2 = round(Cluster_props.Centroid(cl,1));
                            center_z = round(Cluster_props.Centroid(cl,3));

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
                            Cropped_Img_cell{cl}{qq} = thisImage;

                        end
                    end
                end
                paired_Cropped_Img_store{nn} = Cropped_Img_cell;

                % Store properties of parent nuclei
                paired_nucVolume{nn} = ...
                    ones(size(Cluster_Volume_array)) ...
                    .*nuc_volCell{ff}(nn);
                paired_nucClustVolume{nn} = ...
                    cell(size(OP_Volume_array));
                paired_nucClustNNDist{nn} = ...
                    cell(size(OP_Volume_array));
                for object_nn = 1:numel(OP_Volume_array)
                    paired_nucClustVolume{nn}{object_nn} = ...
                        Cluster_Volume_array;
                    paired_nucClustNNDist{nn}{object_nn} = ...
                        LargeS5P_NNDist{nn};
                end
                % --- Store nuclear intensity for all target channels
                for qq = 1:numQuantChannels

                    paired_nucIntensity{qq,nn} = ...
                        ones(size(paired_OP_intensity{qq,nn})) ...
                        .*nuc_intCell{ff}{qq}(nn);

                end





            else
                paired_Cluster_OP_dist{nn} = [];
                paired_Cluster_volume{nn} = [];
                paired_Cluster_solidity{nn} = [];
                paired_Cluster_elongation{nn} = [];
                paired_Cropped_Img_store{nn} = {};
                for qq = 1:numQuantChannels
                    paired_Cluster_intensity{qq,nn} = [];
                    paired_OP_intensity{qq,nn} = [];
                end
                paired_OP_Cluster_disp{nn} = [];

                paired_nucVolume{nn} = [];
                paired_nucClustVolume{nn} = {};
                for qq = 1:numQuantChannels
                    paired_nucIntensity{qq,nn} = [];
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
                squeeze(max(S5P_subImage,[],3)))
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
                squeeze(max(S5P_mask,[],3)))
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


            if OP_comps.NumObjects>0 && Cluster_comps.NumObjects>0

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

    end % ends the loop over nuclei in this file

    % pool from nuclei and store for this file

    S5P_volCell{ff} = vertcat(S5P_volume{:});
    S5P_solCell{ff} = vertcat(S5P_solidity{:});
    S5P_eloCell{ff} = vertcat(S5P_elongation{:});
    S5P_imgCell{ff} = vertcat(S5P_centralSlices_store{:});
    S5P_centCell{ff} = vertcat(S5P_cent_store{:});
    S5P_intCell{ff} = cell(1,numQuantChannels);
    S5P_nucIntCell{ff} = cell(1,numQuantChannels);
    for qq = 1:numQuantChannels
        S5P_intCell{ff}{qq} = vertcat(S5P_intensity{qq,:});
        S5P_nucIntCell{ff}{qq} = vertcat(S5P_nucIntensity{qq,:});
    end
    S5P_nucVolCell{ff} = vertcat(S5P_nucVolume{:});
    S5P_nucClustVolCell{ff} = vertcat(S5P_nucClustVolume{:});

    S2P_volCell{ff} = vertcat(S2P_volume{:});
    S2P_solCell{ff} = vertcat(S2P_solidity{:});
    S2P_eloCell{ff} = vertcat(S2P_elongation{:});
    S2P_imgCell{ff} = vertcat(S2P_centralSlices_store{:});
    S2P_centCell{ff} = vertcat(S2P_cent_store{:});
    S2P_intCell{ff} = cell(1,numQuantChannels);
    S2P_nucIntCell{ff} = cell(1,numQuantChannels);
    for qq = 1:numQuantChannels
        S2P_intCell{ff}{qq} = vertcat(S2P_intensity{qq,:});
        S2P_nucIntCell{ff}{qq} = vertcat(S2P_nucIntensity{qq,:});
    end
    S2P_nucVolCell{ff} = vertcat(S2P_nucVolume{:});
    S2P_nucClustVolCell{ff} = vertcat(S2P_nucClustVolume{:});

    OP_volCell{ff} = vertcat(OP_volume{:});
    OP_solCell{ff} = vertcat(OP_solidity{:});
    OP_eloCell{ff} = vertcat(OP_elongation{:});
    OP_imgCell{ff} = vertcat(OP_centralSlices_store{:});
    OP_centCell{ff} = vertcat(OP_cent_store{:});
    OP_intCell{ff} = cell(1,numQuantChannels);
    OP_nucIntCell{ff} = cell(1,numQuantChannels);
    for qq = 1:numQuantChannels
        OP_intCell{ff}{qq} = vertcat(OP_intensity{qq,:});
        OP_nucIntCell{ff}{qq} = vertcat(OP_nucIntensity{qq,:});
    end
    OP_nucVolCell{ff} = vertcat(OP_nucVolume{:});
    OP_nucClustVolCell{ff} = vertcat(OP_nucClustVolume{:});




    LargeS5P_NNDistCell{ff} = LargeS5P_NNDist;

    paired_Cluster_distCell{ff} = vertcat(paired_Cluster_OP_dist{:});
    paired_Cluster_volCell{ff} = vertcat(paired_Cluster_volume{:});
    paired_Cluster_solCell{ff} = vertcat(paired_Cluster_solidity{:});
    paired_Cluster_eloCell{ff} = vertcat(paired_Cluster_elongation{:});
    paired_Cropped_ImgCell{ff} = vertcat(paired_Cropped_Img_store{:});
    paired_Cluster_intCell{ff} = cell(1,numQuantChannels);
    for qq = 1:numQuantChannels
        paired_Cluster_intCell{ff}{qq} = vertcat(paired_Cluster_intensity{qq,:});
    end
    paired_OP_intCell{ff} = cell(1,numQuantChannels);
    paired_nucIntCell{ff} = cell(1,numQuantChannels);
    for qq = 1:numQuantChannels
        paired_OP_intCell{ff}{qq} = vertcat(paired_OP_intensity{qq,:});
        paired_nucIntCell{ff}{qq} = ...
            vertcat(paired_nucIntensity{qq,:});
    end
    paired_OP_displCell{ff} = vertcat(paired_OP_Cluster_disp{:});
    paired_nucVolCell{ff} = vertcat(paired_nucVolume{:});
    paired_nucClustVolumeCell{ff} = vertcat(paired_nucClustVolume{:});
    paired_nucNNDistCell{ff} = vertcat(paired_nucClustNNDist{:});



    figure(2)
    clf

    subplot(1,3,1)
    % Count of OP detections per nucleus

    OP_countPerNucleus = ...
        cellfun(@(vals)numel(vals),paired_Cluster_OP_dist);
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
    inContact_inds = paired_Cluster_distCell{ff}<=dist_threshold;
    plot(paired_OP_intCell{ff}{S5P_SegChannel}(inContact_inds),...
        paired_OP_intCell{ff}{OP_SegChannel}(inContact_inds),'ko',...
        'MarkerEdgeColor','none','MarkerFaceColor',[1,0,0])
    hold on

    plot(paired_OP_intCell{ff}{S5P_SegChannel},...
        paired_OP_intCell{ff}{OP_SegChannel},...
        'ko')

    hold on

    xlabel('Pol II Ser5P Intensity (a.u.)')
    ylabel('OP Intensity (a.u.)')
    title(sprintf('%d OP foci from %d nuclei',...
        numel(paired_Cluster_distCell{ff}),numNuclei_vec(ff)))
    %set(gca,'XLim',[0,3],'YLim',[0,3])
    hold on
    plot([0,3],[0,3],'k-','LineWidth',1.5)
    hold off

    subplot(1,3,3)
    % Cross talk assessment within Pol II Ser5P clusters
    plot(paired_Cluster_intCell{ff}{S5P_SegChannel}(inContact_inds),...
        paired_Cluster_intCell{ff}{OP_SegChannel}(inContact_inds),'ko',...
        'MarkerEdgeColor','none','MarkerFaceColor',[1,0,0])
    hold on

    plot(paired_Cluster_intCell{ff}{S5P_SegChannel},...
        paired_Cluster_intCell{ff}{OP_SegChannel},...
        'ko')

    hold on

    xlabel('Pol II Ser5P Intensity (a.u.)')
    ylabel('OP Intensity (a.u.)')
    title(sprintf('Nearest neighbor Pol II clusters',...
        numel(paired_Cluster_distCell{ff}),numNuclei_vec(ff)))
    %set(gca,'XLim',[0,3],'YLim',[0,3])
    hold on
    plot([0,3],[0,3],'k-','LineWidth',1.5)
    hold off


    if numNuclei_vec(ff)>0 && poolsize == 0
        %waitforbuttonpress
    end

    if fractionAboveFour > 0
        validFileFlag(ff) = false;
    else
        validFileFlag(ff) = true;
    end

end

hasNuclei_flag = numNuclei_vec>0;

%% To prepare for sorting into conditions, introduce cluster quality control
% The logic is to check the intensity distributions of the large Pol II
% clusters against the overall data set. For all files, where this
% distribution is looking markedly different, that file will be marked for
% not being included in the next section, where the results are sorted
% based on experimental condition.

Vol_threshold = 0.08; % in micrometers

withinLimits_flag = hasNuclei_flag;

all_S5P_vols = vertcat(S5P_volCell{withinLimits_flag});
all_S5P_ints = vertcat(cellfun(@(xx)xx{2},S5P_intCell(withinLimits_flag),...
    'UniformOutput',false));
all_S2P_ints = vertcat(cellfun(@(xx)xx{1},S5P_intCell(withinLimits_flag),...
    'UniformOutput',false));

all_S5P_ints = vertcat(all_S5P_ints{:});
all_S2P_ints = vertcat(all_S2P_ints{:});

all_S5P_ints = all_S5P_ints(all_S5P_vols>=Vol_threshold);
all_S2P_ints = all_S2P_ints(all_S5P_vols>=Vol_threshold);

S5P_limits = prctile(all_S5P_ints,[5,95]);
S2P_limits = prctile(all_S2P_ints,[5,95]);

numValidFiles = sum(withinLimits_flag);
validInds = find(withinLimits_flag);

for vv = 1:numValidFiles

    ff = validInds(vv);

    file_Vols = S5P_volCell{ff};
    file_S5P_ints = S5P_intCell{ff}{2};
    file_S2P_ints = S5P_intCell{ff}{1};

    file_S5P_ints = file_S5P_ints(file_Vols>=Vol_threshold);
    file_S2P_ints = file_S2P_ints(file_Vols>=Vol_threshold);

    max_fraction_outside = max([
        mean(file_S5P_ints<S5P_limits(1)),...
        mean(file_S5P_ints>S5P_limits(2)),...
        mean(file_S2P_ints<S2P_limits(1)),...
        mean(file_S2P_ints>S2P_limits(2))]);

    if max_fraction_outside < 0.2
        withinLimits_flag(ff) = true;
    end

end

%% Sort into conditions

% ---
% Assemble the indices of files to include for all of the unique
% experimental conditions present in the data set

% Get unique condition names
uniqueCondNames = unique(condNames);
numUniqueConds = numel(uniqueCondNames);

%Prepare containers for file indices
fileIndsCell = cell(1,numUniqueConds);
numFiles_perCond = zeros(1,numUniqueConds);
numRejectedFiles_perCond = zeros(1,numUniqueConds);

for cc = 1:numUniqueConds

    indsBasedOnNames = cellfun(...
        @(elmt)strcmp(elmt,uniqueCondNames{cc}),condNames);
    indsCombinedValidity = indsBasedOnNames ...
        & withinLimits_flag & validFileFlag;
    indsRejectedFiles = ~indsCombinedValidity(indsBasedOnNames);
    
	fileIndsCell{cc} = indsBasedOnNames ...
        & withinLimits_flag & validFileFlag;

	numFiles_perCond(cc) = sum(fileIndsCell{cc});
    numRejectedFiles_perCond(cc) = sum(indsRejectedFiles);

    if numRejectedFiles_perCond(cc)>0
        disp('There were files rejected during quality control,')
        disp('condition:')
        disp(uniqueCondNames{cc})
        disp('File names:')
        for ff = find(indsRejectedFiles)
            disp(listing(ff).name)
        end
    end

end
% ---

% ---
% Actually use the indices for file selective file inclusion for the
% different conditions

% Define container variables for all conditions
sortedCondNames = cell(1,numUniqueConds);
sortedNumFiles = zeros(1,numUniqueConds);

sortedVoxelSize_xy = cell(1,numUniqueConds);
sortedVoxelSize_z = cell(1,numUniqueConds);

sortedNumNuclei = zeros(1,numUniqueConds);
sortedNucVol = cell(1,numUniqueConds);
sortedNucInt = cell(1,numQuantChannels);
sortedCytoInt = cell(1,numQuantChannels);
sortedNucStd = cell(1,numQuantChannels);

sortedS5PNumCell = cell(1,numUniqueConds);
sortedS5PVolCell = cell(1,numUniqueConds);
sortedS5PVolPerNucCell = cell(1,numUniqueConds);
sortedS5PSolCell = cell(1,numUniqueConds);
sortedS5PEloCell = cell(1,numUniqueConds);
sortedS5PCentralSliceCell = cell(1,numUniqueConds);
sortedS5PCentroidsCell = cell(1,numUniqueConds);
sortedS5PIntCell = cell(1,numQuantChannels);
sortedS5PNucIntCell = cell(1,numQuantChannels);
sortedS5PNucVolCell = cell(1,numUniqueConds);  
sortedS5PNucClustVolCell = cell(1,numUniqueConds);

sortedS2PNumCell = cell(1,numUniqueConds);
sortedS2PVolCell = cell(1,numUniqueConds);
sortedS2PVolPerNucCell = cell(1,numUniqueConds);
sortedS2PSolCell = cell(1,numUniqueConds);
sortedS2PEloCell = cell(1,numUniqueConds);
sortedS2PCentralSliceCell = cell(1,numUniqueConds);
sortedS2PCentroidsCell = cell(1,numUniqueConds);
sortedS2PIntCell = cell(1,numQuantChannels);
sortedS2PNucIntCell = cell(1,numQuantChannels);
sortedS2PNucVolCell = cell(1,numUniqueConds);

sortedOPNumCell = cell(1,numUniqueConds);
sortedOPVolCell = cell(1,numUniqueConds);
sortedOPVolPerNucCell = cell(1,numUniqueConds);
sortedOPSolCell = cell(1,numUniqueConds);
sortedOPEloCell = cell(1,numUniqueConds);
sortedOPCentralSliceCell = cell(1,numUniqueConds);
sortedOPCentroidsCell = cell(1,numUniqueConds);
sortedOPIntCell = cell(1,numQuantChannels);
sortedOPNucIntCell = cell(1,numQuantChannels);
sortedOPNucVolCell = cell(1,numUniqueConds);  

sorted_S5P_NNDist = cell(1,numUniqueConds);

sorted_paired_ClustVolCell = cell(1,numUniqueConds);
sorted_paired_VolPerNucCell = cell(1,numUniqueConds);
sorted_paired_S5PNNDistPerNucCell = cell(1,numUniqueConds);
sorted_paired_ClustSolCell = cell(1,numUniqueConds);
sorted_paired_ClustEloCell = cell(1,numUniqueConds);
sorted_paired_ClustIntCell = cell(1,numQuantChannels);
sorted_paired_OPIntCell = cell(1,numQuantChannels);
sorted_paired_OPDisplCell = cell(1,numUniqueConds);
sorted_paired_distCell = cell(1,numUniqueConds);
sorted_paired_ImgCell = cell(1,numUniqueConds);
sorted_paired_NucIntCell = cell(1,numQuantChannels);
sorted_paired_NucVolCell = cell(1,numUniqueConds);  

for qq = 1:numQuantChannels

    sortedNucInt{qq} = cell(1,numUniqueConds);
	sortedCytoInt{qq} = cell(1,numUniqueConds);
    sortedNucStd{qq} = cell(1,numUniqueConds);

    sortedS5PIntCell{qq} = cell(1,numUniqueConds);
    sortedS5PNucIntCell{qq} = cell(1,numUniqueConds);

    sortedS2PIntCell{qq} = cell(1,numUniqueConds);
    sortedS2PNucIntCell{qq} = cell(1,numUniqueConds);

    sortedOPIntCell{qq} = cell(1,numUniqueConds);
    sortedOPNucIntCell{qq} = cell(1,numUniqueConds);

    sorted_paired_ClustIntCell{qq} = cell(1,numUniqueConds);
    sorted_paired_OPIntCell{qq} = cell(1,numUniqueConds);
    sorted_paired_NucIntCell{qq} = cell(1,numUniqueConds);

end




for cc = 1:numUniqueConds

    % These indices will be used to pull all the data from all valid files
    % that are part of this unique experimental condition
    thisCondInds = fileIndsCell{cc};

    % Name of condition and number of included files (with quality control)
    sortedCondNames{cc} = uniqueCondNames{cc};
    sortedNumFiles(cc) = sum(thisCondInds);

    % --- general file properties, here mostly physical distance between
    % voxels, so transfor to physical units becomes possible

    sortedVoxelSize_xy{cc} = voxelSize_array_XY(thisCondInds);
    sortedVoxelSize_z{cc} = voxelSize_array_Z(thisCondInds);

    % --- properties of whole nuclei

    sortedNumNuclei(cc) = sum(numNuclei_vec(thisCondInds));

    sortedNucVol{cc} = vertcat(nuc_volCell{thisCondInds});

    for qq = 1:numQuantChannels

        sortedNucInt{qq}{cc} = vertcat(arrayfun(...
            @(ind)nuc_intCell{ind}{qq},....
            find(thisCondInds),...
            'UniformOutput',false));
        sortedNucInt{qq}{cc} = vertcat(sortedNucInt{qq}{cc}{:})';

        sortedCytoInt{qq}{cc} = vertcat(arrayfun(...
            @(ind)cyto_intCell{ind}{qq},....
            find(thisCondInds),...
            'UniformOutput',false));
        sortedCytoInt{qq}{cc} = vertcat(sortedCytoInt{qq}{cc}{:})';

        sortedNucStd{qq}{cc} = vertcat(arrayfun(...
            @(ind)nuc_stdCell{ind}{qq},....
            find(thisCondInds),...
            'UniformOutput',false));
        sortedNucStd{qq}{cc} = vertcat(sortedNucStd{qq}{cc}{:})';

    end

    % ---
    % Properties of detected objects irrespective of being associated with
    % specific oligopaint channel based objects

    % - properties of Ser5P clusters
	S5P_vols = vertcat(S5P_volCell{thisCondInds});
    S5P_sols = vertcat(S5P_solCell{thisCondInds});
	S5P_elos = vertcat(S5P_eloCell{thisCondInds});
	S5P_slices = vertcat(S5P_imgCell{thisCondInds});
	S5P_centroids = vertcat(S5P_centCell{thisCondInds});
	S5P_ints = S5P_intCell(thisCondInds);
    % - properties connected to parent nuclei
    S5P_nucClustVols = S5P_nucClustVolCell(thisCondInds);
    S5P_numPerNucleus = vertcat(arrayfun(...
        @(val)perNuc_countCell{val}{1},find(thisCondInds),...
		'UniformOutput',false));
    S5P_numPerNucleus = vertcat(S5P_numPerNucleus{:});
    % - properties of parent nuclei
    S5P_nucInts = S5P_nucIntCell(thisCondInds);
    S5P_nucVols = vertcat(S5P_nucVolCell{thisCondInds});

	sortedS5PNumCell{cc} = S5P_numPerNucleus;
    sortedS5PVolCell{cc} = S5P_vols;
    sortedS5PSolCell{cc} = S5P_sols;
	sortedS5PEloCell{cc} = S5P_elos;
	sortedS5PCentralSliceCell{cc} = S5P_slices;
	sortedS5PCentroidsCell{cc} = S5P_centroids;
    sortedS5PNucVolCell{cc} = S5P_nucVols;
    
    % - properties of Ser2P foci
	S2P_vols = vertcat(S2P_volCell{thisCondInds});
    S2P_sols = vertcat(S2P_solCell{thisCondInds});
	S2P_elos = vertcat(S2P_eloCell{thisCondInds});
	S2P_slices = vertcat(S2P_imgCell{thisCondInds});
	S2P_centroids = vertcat(S2P_centCell{thisCondInds});
	S2P_ints = S2P_intCell(thisCondInds);
    % - properties connected to parent nuclei
    S2P_numPerNucleus = vertcat(arrayfun(...
        @(val)perNuc_countCell{val}{2},find(thisCondInds),...
		'UniformOutput',false));
    S2P_numPerNucleus = vertcat(S2P_numPerNucleus{:});
    % - properties of parent nuclei
    S2P_nucInts = S2P_nucIntCell(thisCondInds);
    S2P_nucVols = vertcat(S2P_nucVolCell{thisCondInds});
    
    sortedS2PNumCell{cc} = S2P_numPerNucleus;
    sortedS2PVolCell{cc} = S2P_vols;
    sortedS2PSolCell{cc} = S2P_sols;
	sortedS2PEloCell{cc} = S2P_elos;
	sortedS2PCentralSliceCell{cc} = S2P_slices;
	sortedS2PCentroidsCell{cc} = S2P_centroids;
    sortedS2PNucVolCell{cc} = S2P_nucVols;

    % - properties of Ser2P foci
	OP_vols = vertcat(OP_volCell{thisCondInds});
    OP_sols = vertcat(OP_solCell{thisCondInds});
	OP_elos = vertcat(OP_eloCell{thisCondInds});
	OP_slices = vertcat(OP_imgCell{thisCondInds});
	OP_centroids = vertcat(OP_centCell{thisCondInds});
	OP_ints = paired_OP_intCell(thisCondInds);
    % - properties connected to parent nuclei
    OP_numPerNucleus = vertcat(arrayfun(...
        @(val)perNuc_countCell{val}{3},find(thisCondInds),...
		'UniformOutput',false));
    OP_numPerNucleus = vertcat(OP_numPerNucleus{:});
    % - properties of parent nuclei
    OP_nucInts = OP_nucIntCell(thisCondInds);
    OP_nucVols = vertcat(OP_nucVolCell{thisCondInds});
    
    sortedOPNumCell{cc} = OP_numPerNucleus;
    sortedOPVolCell{cc} = OP_vols;
    sortedOPSolCell{cc} = OP_sols;
	sortedOPEloCell{cc} = OP_elos;
	sortedOPCentralSliceCell{cc} = OP_slices;
	sortedOPCentroidsCell{cc} = OP_centroids;
    sortedOPNucVolCell{cc} = OP_nucVols;

    for qq = 1:numQuantChannels

        sortedS5PIntCell{qq}{cc} = vertcat(arrayfun(...
            @(ind)S5P_intCell{ind}{qq},....
            find(thisCondInds),...
            'UniformOutput',false));
        sortedS5PIntCell{qq}{cc} = vertcat(sortedS5PIntCell{qq}{cc}{:});

        sortedS5PNucIntCell{qq}{cc} = vertcat(arrayfun(...
            @(ind)S5P_nucIntCell{ind}{qq},....
            find(thisCondInds),...
            'UniformOutput',false));
        sortedS5PNucIntCell{qq}{cc} = vertcat(sortedS5PNucIntCell{qq}{cc}{:});

        sortedS2PIntCell{qq}{cc} = vertcat(arrayfun(...
            @(ind)S2P_intCell{ind}{qq},....
            find(thisCondInds),...
            'UniformOutput',false));
        sortedS2PIntCell{qq}{cc} = vertcat(sortedS2PIntCell{qq}{cc}{:});

        sortedS2PNucIntCell{qq}{cc} = vertcat(arrayfun(...
            @(ind)S2P_nucIntCell{ind}{qq},....
            find(thisCondInds),...
            'UniformOutput',false));
        sortedS2PNucIntCell{qq}{cc} = vertcat(sortedS2PNucIntCell{qq}{cc}{:});

        sortedOPIntCell{qq}{cc} = vertcat(arrayfun(...
            @(ind)paired_OP_intCell{ind}{qq},....
            find(thisCondInds),...
            'UniformOutput',false));
        sortedOPIntCell{qq}{cc} = vertcat(sortedOPIntCell{qq}{cc}{:});

        sortedOPNucIntCell{qq}{cc} = vertcat(arrayfun(...
            @(ind)OP_nucIntCell{ind}{qq},....
            find(thisCondInds),...
            'UniformOutput',false));
        sortedOPNucIntCell{qq}{cc} = vertcat(sortedOPNucIntCell{qq}{cc}{:});

    end

    % ---
    % Properties of oligopaint-based observation instances of oligopaint
    % labeled genes connected to the nearest large Pol II Ser5P cluster


    sorted_S5P_NNDist{cc} = LargeS5P_NNDistCell(thisCondInds);
    sorted_S5P_NNDist{cc} = [sorted_S5P_NNDist{cc}{:}];
    
    sorted_paired_distCell{cc} = ...
        vertcat(paired_Cluster_distCell{thisCondInds});
    sorted_paired_ClustVolCell{cc} = ...
        vertcat(paired_Cluster_volCell{thisCondInds});
    sorted_paired_ClustVolCell{cc} = ...
        vertcat(paired_Cluster_volCell{thisCondInds});
    sorted_paired_ClustSolCell{cc} = ...
        vertcat(paired_Cluster_solCell{thisCondInds});
    sorted_paired_ClustEloCell{cc} = ...
        vertcat(paired_Cluster_eloCell{thisCondInds});
    sorted_paired_ImgCell{cc} = ...
        vertcat(paired_Cropped_ImgCell{thisCondInds});
    sorted_paired_OPDisplCell{cc} = ...
        vertcat(paired_OP_displCell{thisCondInds});
    sorted_paired_NucIntCell{cc} = ...
         paired_nucIntCell(thisCondInds);
    sorted_paired_NucVolCell{cc} = ...
         vertcat(paired_nucVolCell{thisCondInds});

    sorted_paired_VolPerNucCell{cc} = ...
        paired_nucClustVolumeCell(thisCondInds);
    sorted_paired_S5PNNDistPerNucCell{cc} = ...
        paired_nucNNDistCell(thisCondInds);


    for qq = 1:numQuantChannels

        sorted_paired_ClustIntCell{qq}{cc} = vertcat(arrayfun(...
            @(ind)paired_Cluster_intCell{ind}{qq},....
            find(thisCondInds),...
            'UniformOutput',false));
        sorted_paired_ClustIntCell{qq}{cc} = vertcat(sorted_paired_ClustIntCell{qq}{cc}{:});

        sorted_paired_OPIntCell{qq}{cc} = vertcat(arrayfun(...
            @(ind)paired_OP_intCell{ind}{qq},....
            find(thisCondInds),...
            'UniformOutput',false));
        sorted_paired_OPIntCell{qq}{cc} = vertcat(sorted_paired_OPIntCell{qq}{cc}{:});

        sorted_paired_NucIntCell{qq}{cc} = vertcat(arrayfun(...
            @(ind)paired_nucIntCell{ind}{qq},....
            find(thisCondInds),...
            'UniformOutput',false));
        sorted_paired_NucIntCell{qq}{cc} = vertcat(sorted_paired_NucIntCell{qq}{cc}{:});

    end

end

% End of sorting individual file data into the unique conditions
% ---


%% Saving of results

if visualizationFlag
    save('Visualization_ConditionSortedResults_Test')
else
    save('ConditionSortedResults_Test')
end

%% Overview plots for all analyzed conditions

figure(4)
clf

prct_dist = cellfun(@(xx)prctile(xx,10),sorted_paired_distCell); % in micrometers
prct_Ser2P = cellfun(@(xx)prctile(xx,90),sorted_paired_ClustIntCell{1});

contact_freq = cellfun(@(xx)mean(xx<=0.50),sorted_paired_distCell); % in micrometers);
avg_dist = cellfun(@mean,sorted_paired_distCell); % in micrometers
avg_Ser5P = cellfun(@mean,sorted_paired_ClustIntCell{2});
avg_Ser2P = cellfun(@mean,sorted_paired_ClustIntCell{1});

scatter(avg_Ser5P,avg_Ser2P,100,prct_dist,'filled')
xlabel('Mean Pol II Ser5P')
ylabel('Mean Pol II Ser2P')

colormap(parula)
colorbar
set(gca,'Box','on')
title('Mean gene-cluster distance [\mum]','FontWeight','normal')

