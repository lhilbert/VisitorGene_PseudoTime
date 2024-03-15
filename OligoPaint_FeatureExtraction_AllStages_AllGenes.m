
clear all

%% Analysis parameters

% Specify the directory that contains the extraced files from the last
% step, where you extracted from the raw files obtained from the microscope

sourceDirectory = './ExtractedStacks/**/'; %For all stages
% sourceDirectory = './ExtractedStacks_Stages/Sphere/**/'; %Specific stage
% sourceDirectory = './ExtractedStacks_Stages/Sphere/Cond_11/'; %Test run

% Channels for segmentation (Order defined by the extraction script, not by
% the original data
NucSegChannel = 2; %Channel used to detect nuclei (S5P)
ClusterSegChannel = 2; %Channel used to detect Pol II S5P clusters
OPSegChannel = 3; % Channel used to detect gene OligoPaint

% Target channels for intensity quantification, applied for all objects
quantChannels = [1,2,OPSegChannel]; %Define specific order
quantBlurSigma = [0,0,0.07]; %Why quantblur sigma=0.07 for 

%Blur for nuclei segmentation
segBlurSigma_small = 1.0; % in microns (original:1)
segBlurSigma_large = 10; % in microns (original:10)

%Blur for S5P segmentation
clusterSegBlurSigma_large = 3.0; %backgrounbd subtraction 3.0
clusterSeg_numStdDev = 2.0;%robust background thresholing: two standard 2.0
% deviations above intensity mean) also used: 5 and 2

%Blur for OP segmentation
OPsegBlurSigma_small = 0.1; % in microns 0.1
OPsegBlurSigma_large = 5.0; % in microns 5
OPseg_numStdDev = 6;%6; % number of standard deviations in robust threshold 6

%Minimal properties for Nuclei, S5P and OP segmentation
Nuc_min_vol = 40; % Minimal nucleus volume [cubic microns] original:40
Nuc_min_sol = 0.7; % Minimal nucleus solidity to ensure round nuclei 0.7
Cluster_minVol = 0.03; %Minimal cluster volume [cubic microns?] 0.03
OP_minVol = 0.05; %Minimal OP colume [cubic microns] 0.05

%Keep track of OP vol and OP intensity signal

%What is it for?
centralSliceExtension = 0; % pixels from centroid

% Cluster connection range: Not fully implemented yet
cluster_DBSCAN_epsilon = 0.65;%0.65; % in microns, choose 0 for no clustering

%PCA parameters
dist_threshold = 0.25;
Vol_threshold = 0.2;
angle_shift = -0.1;

%End of analysis parameter section. Do not change anything else in
%this section, all necessary parameters are listed above.

% Extract data from previous extracted stacks

%Get all extracted microscopy files and their total number
listing = rdir([sourceDirectory,'*Image*.mat']); 
numFiles = numel(listing);

% Condition index and name retrieval
condInds = []; %Condition indices
condNames = {}; %Condition names
%Loop all files
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

% Analyze the image stacks one by one

% Variables to store properties of nuclei
numNuclei_vec = zeros(1,numFiles);
Cluster_distCell = cell(1,numFiles);
Cluster_volCell = cell(1,numFiles);
Cluster_solCell = cell(1,numFiles);
Cluster_eloCell = cell(1,numFiles);
Cluster_maskCell = cell(1,numFiles);
Cluster_intCell = cell(1,numFiles);
% Nucleus_intCell = cell(1,numFiles);
OP_intCell = cell(1,numFiles);
OP_BG_intCell = cell(1,numFiles);
OP_displCell = cell(1,numFiles);
pixelSize_array = zeros(1,numFiles);
zStep_array = zeros(1,numFiles);

% Number of channels for quantification
numQuantChannels = numel(quantChannels);

%Loop trough all files one by one. For parallel loop use parfor instead of 
%for. But within parfor no plotting is possible.

figure(6)
clf

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
	
    %Perform feature enhancement - Noise removal by small sigma Gaussian 
    %and background subtraction by large sigma Gaussian -> highlights the
    %features between the two sizes

    segImg = ...
		+ imgaussfilt(segImg,segBlurSigma_small./pixelSize) ...
		- imgaussfilt(segImg,segBlurSigma_large./pixelSize);

	%Histogram of intensity values
    [bin_counts,bin_centers] = hist(segImg(:),1000);
	
    %Apply otsu threshold to the histogram values
    [nuc_seg_thresh,~] = otsuLimit(bin_centers,bin_counts,[0,Inf]);
	
    %Define nucleus segmentation mask with calculates otsu threshold
    NucSegMask = segImg > 1.0.*nuc_seg_thresh;
	
%     %Plot different steps of image processing in combined subplot
%     
%     %Original image
% 	subplot(2,2,1)
% 	imagesc(squeeze(imgStack{NucSegChannel}(:,:,ceil(imgSize(3)./2)))) 
% 	axis tight equal
%     title('Original',FontSize=20); 
%     
%     %Processed image
% 	subplot(2,2,2)
% 	imagesc(squeeze(segImg(:,:,ceil(imgSize(3)./2))))
% 	axis tight equal
%     title('Processed',FontSize=20);
%     
%     %Segmentation mask
% 	subplot(2,2,3)
% 	imagesc(squeeze(NucSegMask(:,:,ceil(imgSize(3)./2))))
% 	axis tight equal
%     title('Segmented',FontSize=20)
% 
%     %Original image with boundary of segmentation mask
%     subplot(2,2,4)
% 	imagesc(squeeze(imgStack{NucSegChannel}(:,:,ceil(imgSize(3)./2)))) 
% 	axis tight equal;
%     hold on;
%     visboundaries(squeeze(NucSegMask(:,:,ceil(imgSize(3)./2))), ...
%         'Color', 'white')
%     title('Original + Segmented',FontSize=20);
%     hold off;
% 
%  	waitforbuttonpress
	
	% --- Connected component segmentation of nuclei
	comps = bwconncomp(NucSegMask,18);
	numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList); %count number
    % of elements/pixels for every component
% 	fprintf('No. objects unfiltered: %d \n', max(size(numPxls)))
    
    %Volume filter

    minPixels = Nuc_min_vol./(pixelSize.^2)./zStepSize;
	comps.NumObjects = sum(numPxls>=minPixels);
% 	fprintf('No. objects vol. filter: %d \n', comps.NumObjects)
    comps.PixelIdxList = comps.PixelIdxList(numPxls>=minPixels);
	numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList);
	
	%Solidity filter
    props = regionprops3(comps,imgStack{NucSegChannel},...
		'Solidity');
	
	Solidity_array = [props.Solidity]; %Faster to filter first for solidity
    %and then again for the rest?
	
	comps.NumObjects = sum(Solidity_array>=Nuc_min_sol);
	comps.PixelIdxList = comps.PixelIdxList(Solidity_array>=Nuc_min_sol);
	numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList);
%     fprintf('No. objects sol. filter: %d \n', comps.NumObjects)

    %waitforbuttonpress
    
%     %Plotof all segmented z-planes
%     subplot(1,1,1)
%     montage(reshape(imgStack{NucSegChannel},[size(NucSegMask,1),size(NucSegMask,2),1,size(NucSegMask,3)]),'displayRange' , []);
%     title(['Original image all z-planes file', num2str(ff)], FontSize=20);
%     waitforbuttonpress
% 
%     subplot(1,1,1)
%     montage(reshape(NucSegMask,[size(NucSegMask,1),size(NucSegMask,2),1,size(NucSegMask,3)]),'displayRange' , []);
%     title(['Segmentation mask all z-planes file', num2str(ff)], FontSize=20);
%     waitforbuttonpress
    
%     % Show all z-planes of file with segmented nuclei
%     nS   = sqrt(imgSize(3));
%     nCol = ceil(nS);
%     nRow = nCol - (nCol * nCol - imgSize(3) > nCol - 1);
%     t = tiledlayout(nRow,nCol,'TileSpacing','none');
%     for i = 1:imgSize(3)
%         %subplot(nRow, nCol, i);
%         nexttile
%         imagesc(squeeze(imgStack{NucSegChannel}(:,:,i)));
%         axis equal;
%         axis off;
%         hold on;
%         hold on;
%         visboundaries(squeeze(NucSegMask(:,:,i)),'Color','white');
%     end
%     sgtitle(['Label control (original + boundaries) file',num2str(ff)],FontSize=20);
% 
%     waitforbuttonpress

    %Analysis of segmented 3D connected components 
    props = regionprops3(comps,imgStack{NucSegChannel},...
		'Volume','VoxelValues','Solidity','VoxelIdxList',...
		'BoundingBox');
	
    %Component volume, median intensity and solidity
	Volume_array = [props.Volume].*pixelSize.^2.*zStepSize;
	Intensity_array = cellfun(@(vals)median(vals),props.VoxelValues);
	Solidity_array = [props.Solidity];
	
	%Final number of detected nuclei
	numNuclei = numel(Intensity_array); %Get actual number
    %fprintf('Number of segmented nuclei: %d \n',numNuclei) %Display
	numNuclei_vec(ff) = numNuclei; %Save number within array
	
	% --- For each nucleus, get clusters and oligopaint images
	Cluster_img = imgStack{ClusterSegChannel};
	OP_img = imgStack{OPSegChannel};
	
    %Empty cells for quantifications
	Cluster_OP_dist = cell(1,numNuclei);
	Cluster_volume = cell(1,numNuclei);
	Cluster_solidity = cell(1,numNuclei);
	Cluster_elongation = cell(1,numNuclei);
	Cluster_mask_store = cell(1,numNuclei);
% 	Nucleus_intensity = cell(numQuantChannels,numNuclei); %Added
	Cluster_intensity = cell(numQuantChannels,numNuclei);
	OP_intensity = cell(numQuantChannels,numNuclei);
    OP_BG_intensity = cell(numQuantChannels,numNuclei);
	
    %Loop through every nuclei
	for nn = 1:numNuclei
		%Get bounding box of segmented nuclei
		boxArray = props.BoundingBox(nn,:); %array with box borders
		
        %Extract subimages for cluster and op detection
        %Change to clusterSegChannel
        Cluster_subImage = imgStack{NucSegChannel}(...
			boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
			boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
			boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
		OP_subImage = OP_img(...
			boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
			boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
			boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
		
        %Processing of cluster image (only background - large filter)
		Cluster_subImage = Cluster_subImage ...
			- imgaussfilt(Cluster_subImage,...
			clusterSegBlurSigma_large./pixelSize);		
		

        %Processing of OP image (small and large filter)
		OP_subImage = imgaussfilt(OP_subImage,...
			OPsegBlurSigma_small./pixelSize);
		OP_subImage = OP_subImage - imgaussfilt(OP_subImage,...
			OPsegBlurSigma_large./pixelSize);
		
        %Create nuclear segmentation masks
		NucMask = false(imgSize);
		NucMask(props.VoxelIdxList{nn}) = true;
		NucMask_subImage = NucMask(...
			boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
			boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
			boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
		
        %Get S5P cluster intensities, mean, std and mask above threshold
		Cluster_intensities = Cluster_subImage(NucMask_subImage);
		Cluster_mean = mean(Cluster_intensities);
		Cluster_std = std(Cluster_intensities);
		Cluster_mask = (Cluster_subImage.*NucMask_subImage)...
			>(Cluster_mean+clusterSeg_numStdDev.*Cluster_std);
		
        %Get OP cluster intensities, mean, std and mask above threshold
		OP_intensities = OP_subImage(NucMask_subImage);
		OP_mean = mean(OP_intensities);
		OP_std = std(OP_intensities);
		OP_mask = (OP_subImage.*NucMask_subImage)...
			>(OP_mean+OPseg_numStdDev.*OP_std);
		
		subImgSize = size(Cluster_subImage);
		if numel(subImgSize)==2
			subImgSize(3)=1;
        else

%             %Plot subimages and masks (Max projection)
%             %S5P subimage
% 			subplot(2,2,1)
% 			imagesc(squeeze(max(Cluster_subImage,[],3)))
% 			axis tight equal
%             title('S5P subimage',FontSize=20);
% 			
%             %OP subimage
% 			subplot(2,2,2)
% 			imagesc(squeeze(max(OP_subImage,[],3)))
% 			axis tight equal
%             title('OP subimage',FontSize=20);
% 			
% 			%S5P mask
% 			subplot(2,2,3)
% 			imagesc(squeeze(max(Cluster_mask,[],3)))
% 			axis tight equal
%             title('S5P segmentation mask',FontSize=20);
% 			
%             %OP mask
% 			subplot(2,2,4)
% 			imagesc(squeeze(max(OP_mask,[],3)))
% 			axis tight equal
%             title('OP segmentation mask',FontSize=20);
%             
%             waitforbuttonpress

%             %Plot subimages and masks (mid section)
%             centerPlaneInd = round(boxArray(6).*0.5);
%     
%             %S5P subimage
% 			subplot(2,2,1)
% 			imagesc(squeeze(Cluster_subImage(:,:,centerPlaneInd)))
% 			axis tight equal
%             title('S5P subimage',FontSize=20);
% 			
%             %OP subimage
% 			subplot(2,2,2)
% 			imagesc(squeeze(OP_subImage(:,:,centerPlaneInd)))
% 			axis tight equal
%             title('OP subimage',FontSize=20);
% 			
% 			%S5P mask
% 			subplot(2,2,3)
% 			imagesc(squeeze(Cluster_mask(:,:,centerPlaneInd)))
% 			axis tight equal
%             title('S5P segmentation mask',FontSize=20);
% 			
%             %OP mask
% 			subplot(2,2,4)
% 			imagesc(squeeze(OP_mask(:,:,centerPlaneInd)))
% 			axis tight equal
%             title('OP segmentation mask',FontSize=20);
% 
%             waitforbuttonpress

		end
				
 		%waitforbuttonpress
	
        % Connected components of OP and Cluster masks

        %Oligopaint
		OP_comps = bwconncomp(OP_mask,18);
		OP_numPxls = cellfun(@(elmt)numel(elmt),OP_comps.PixelIdxList);
		minPixels = OP_minVol./(pixelSize.^2)./zStepSize;
		OP_comps.NumObjects = sum(OP_numPxls>minPixels);
% 		fprintf('Number of segmented Oligopaint foci: %d\n',OP_comps.NumObjects)
        OP_comps.PixelIdxList = OP_comps.PixelIdxList(OP_numPxls>minPixels);
		OP_numPxls = cellfun(@(elmt)numel(elmt),OP_comps.PixelIdxList);
		
        %S5P cluster
		Cluster_comps = bwconncomp(Cluster_mask,18);
		Cluster_numPxls = cellfun(@(elmt)numel(elmt),Cluster_comps.PixelIdxList);
		minPixels = Cluster_minVol./(pixelSize.^2)./zStepSize;
		Cluster_comps.NumObjects = sum(Cluster_numPxls>minPixels);
%         fprintf('Number of segmented clusters: %d\n',Cluster_comps.NumObjects)
		Cluster_comps.PixelIdxList = Cluster_comps.PixelIdxList(Cluster_numPxls>minPixels);
		Cluster_numPxls = cellfun(@(elmt)numel(elmt),Cluster_comps.PixelIdxList);
		
		%Check if S5P cluster and OP are detected and if yes, measure their
        %properties
		if OP_comps.NumObjects>0 && Cluster_comps.NumObjects>0
			
            %Apply hierarchical clustering with dbscan
            if cluster_DBSCAN_epsilon > 0
                Cluster_props = regionprops3(Cluster_comps,Cluster_subImage,...
                    'Centroid');
                centroid_coords = ...
                    Cluster_props.Centroid.*[pixelSize,pixelSize,zStepSize];
                dbscan_inds = ...
                    dbscan(centroid_coords,cluster_DBSCAN_epsilon,1);

                unique_inds = unique(dbscan_inds);
                num_inds = numel(unique_inds);
                updated_comps = Cluster_comps;
                updated_comps.NumObjects = num_inds;
                updated_comps.PixelIdxList = cell(1,num_inds);
                for ii = 1:num_inds
                    updated_comps.PixelIdxList{ii} = ...
                        sort(vertcat(Cluster_comps.PixelIdxList{...
                        dbscan_inds==unique_inds(ii)} ...
                        ));
                end
                Cluster_comps = updated_comps;
            end
            
%             %Display OP with CLuster
%             figure(6)
%             LL1 = labelmatrix(OP_comps); %Label OP
%             LL2 = labelmatrix(Cluster_comps); %Label Cluster
%             centerPlaneInd = round(boxArray(6).*0.5);
%     
%             subplot(3,3,1)
% 	        imagesc(squeeze(max(OP_mask,[],3)))
% 	        axis tight equal
%             title('Max OP mask',FontSize=10);
%             
% 	        subplot(3,3,2)
% 	        imagesc(squeeze(max(LL1,[],3)))
% 	        axis tight equal
%             title('Max OP comps',FontSize=10);
%             
%             subplot(3,3,3)
%             imagesc(squeeze(LL1(:,:,centerPlaneInd)))
%             axis tight equal
%             title('Mid OP comps',FontSize=10);
%             %set(gca,'Colormap',lines)
% 
%             subplot(3,3,4)
% 	        imagesc(squeeze(max(Cluster_mask,[],3)))
% 	        axis tight equal
%             title('Max Cluster mask',FontSize=10);
%             
% 	        subplot(3,3,5)
% 	        imagesc(squeeze(max(LL2,[],3)))
% 	        axis tight equal
%             title('Max Cluster comps',FontSize=10);
%             
%             subplot(3,3,6)
%             imagesc(squeeze(LL2(:,:,centerPlaneInd)))
%             axis tight equal
%             title('Mid Cluster comps',FontSize=10);
% 
%             subplot(3,3,7)
% 	        imagesc(squeeze(max(OP_mask,[],3)))
% 	        axis tight equal
%             hold on;
%             visboundaries(squeeze(max(Cluster_mask,[],3)), ...
%                 'Color', 'white')
%             title('Max OP mask with Cluster',FontSize=10);
%             hold off;
%             
% 	        subplot(3,3,8)
% 	        imagesc(squeeze(max(LL1,[],3)))
% 	        axis tight equal
%             hold on;
%             visboundaries(squeeze(max(LL2,[],3)), ...
%                 'Color', 'white')
%             title('Max OP comps with Cluster',FontSize=10);
%             hold off;
% 
%             subplot(3,3,9)
%             imagesc(squeeze(LL1(:,:,centerPlaneInd)))
%             axis tight equal
%             hold on;
%             visboundaries(squeeze(LL2(:,:,centerPlaneInd)), ...
%                 'Color', 'white')
%             title('Mid OP comps with Cluster',FontSize=10);
%             hold off;

%             %Display hierarchical clustering
%             LL = labelmatrix(Cluster_comps);
%             
%             centerPlaneInd = round(boxArray(6).*0.5);
% 
%             subplot(1,3,1)
% 		    imagesc(squeeze(max(Cluster_mask,[],3)))
% 		    axis tight equal
%             title('Max cluster mask',FontSize=20);
%             
% 		    subplot(1,3,2)
% 		    imagesc(squeeze(max(LL,[],3)))
% 		    axis tight equal
%             title('Max cluster comps',FontSize=20);
%             
%             subplot(1,3,3)
%             imagesc(squeeze(LL(:,:,centerPlaneInd)))
%             axis tight equal
%             title('Mid cluster comps',FontSize=20);
%             %set(gca,'Colormap',lines)
%             
%             waitforbuttonpress
           
            
            
            
            
            %OP properties
            OP_props = regionprops3(OP_comps,Cluster_subImage,...
				'Centroid','VoxelIdxList');
			
            %Cluster properties
			Cluster_props = regionprops3(Cluster_comps,Cluster_subImage,...
				'Volume','MeanIntensity','Solidity',...
				'Centroid','Image','BoundingBox');
			
            %Measure euclidean distance between center of Cluster and Gene
			pwDist = pdist2(...
				OP_props.Centroid.*[pixelSize,pixelSize,zStepSize],...
				Cluster_props.Centroid.*[pixelSize,pixelSize,zStepSize]);
			%Get minimal distances and according indices
            [minDist,minInds] = min(pwDist,[],2);
           

            %Define minimal distances and indices
			Cluster_OP_minDist = minDist;
			Cluster_OP_minInd = minInds;
			
            %Filter for cluster properties of those which are closest to
            %genes
			Cluster_props = Cluster_props(Cluster_OP_minInd,:);
			Cluster_comps.NumObjects = numel(Cluster_OP_minDist);
			Cluster_comps.PixelIdxList = ...
				Cluster_comps.PixelIdxList(Cluster_OP_minInd);
			Cluster_Volume_array = ...
				[Cluster_props.Volume].*pixelSize.^2.*zStepSize;
			Cluster_Solidity_array = [Cluster_props.Solidity];
			
			%Get cluster central plane and elongation in this plane
            %Empty array for elongation
            Cluster_Elongation_array = ...
				zeros(size(Cluster_Solidity_array));
			
            %Cluster_Volume_array

            %Loop through all clusters, measure major and minor axis and
            %calculate finally the elongation (major/minor)
            for cl = 1:numel(Cluster_Volume_array)
				thisMask = Cluster_props.Image{cl};
% 				centerInd = ceil(size(thisMask,3)./2);
% 				thisMask = squeeze(thisMask(:,:,centerInd)); %Center Plane
				thisMask = squeeze(max(thisMask,[],3)); %Center Plane  
                thisProps = regionprops(uint8(thisMask),...
					'MajorAxisLength','MinorAxisLength');
				Cluster_Elongation_array(cl) = ...
					thisProps.MajorAxisLength./thisProps.MinorAxisLength;	
            end

            %Save all measured properties for each nuclei within the 
            %respective arrays for each nucleus
			Cluster_OP_dist{nn} = Cluster_OP_minDist;
            Cluster_OP_dist{nn} = Cluster_OP_minDist;
%             fprintf('Minimal distance: %d \n',Cluster_OP_minDist)
%             waitforbuttonpress

			Cluster_volume{nn} = Cluster_Volume_array;
			Cluster_solidity{nn} = Cluster_Solidity_array;
			Cluster_elongation{nn} = Cluster_Elongation_array;

			%Quantification for all target channels
			Cluster_Mask_cell = ...
				cell(size(Cluster_Solidity_array));
			
            for cl = 1:numel(Cluster_Volume_array)
				Cluster_Mask_cell{cl} = cell(1,numQuantChannels);
			end
			
            %Loop through all quant channels
            for qq = 1:numQuantChannels
			    %Get index of channel to quantify
                channelInd = quantChannels(qq);
				%Get corresponding subimage
                quant_subImage = imgStack{channelInd}(...
					boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
					boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
					boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
				
                %If sigma for channel given/larger 0, perform blurr
				if quantBlurSigma(qq)>0
					quant_subImage = imgaussfilt(quant_subImage,...
						quantBlurSigma(qq)./pixelSize);
				end
				
                %Calculate nucleus median
				Quant_nucleusMedian = ...
					median(quant_subImage(NucMask_subImage));
				
                %Get cluster intensities
				Cluster_quant_props = regionprops3(...
					Cluster_comps,quant_subImage,'VoxelValues');
				
                %Get cluster median intensities
                Quant_ClusterMedian = cellfun(...
					@(vals)median(vals),Cluster_quant_props.VoxelValues);

                %Normalize cluster intenities by nucleus median
                Cluster_intensity{qq,nn} = ...
					Quant_ClusterMedian./Quant_nucleusMedian;

                %Normalize nucleus intenities by nucleus median
%                 Nucleus_intensity{qq,nn} = ...
% 					quant_subImage(NucMask_subImage)./Quant_nucleusMedian;
                    
                %quant_subImage(NucMask_subImage)./Quant_nucleusMedian
				
                %Get OP intensities
				OP_quant_props = regionprops3(...
					OP_comps,quant_subImage,'VoxelValues');
				
                %Get OP median intensities
                Quant_OPMedian = cellfun(...
					@(vals)median(vals),OP_quant_props.VoxelValues);
				
                %Normalize OP intenities by nucleus median
                OP_intensity{qq,nn} = ...
					Quant_OPMedian./Quant_nucleusMedian;
                
                %Get nuclear mask intensities
                Nuc_int = quant_subImage(NucMask_subImage);
                
                %Get OP Background (random choice) median intensities
                OP_BG_intensity{qq,nn} = arrayfun(...
                    @(kk)median(Nuc_int(...
                    randi(numel(Nuc_int),[1,kk]))),...
                    cellfun(@(vals)numel(vals),...
                    OP_quant_props.VoxelValues))...
                    ./Quant_nucleusMedian;
                
				
				%What is it good for?
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
			%Save the cluster masks
            Cluster_mask_store{nn} = Cluster_Mask_cell;
		
        %If no OP or S5P cluster detected, leave empty
		else
			Cluster_OP_dist{nn} = [];
			Cluster_volume{nn} = [];
			Cluster_solidity{nn} = [];
			Cluster_elongation{nn} = [];
			Cluster_mask_store{nn} = {};
			for qq = 1:numQuantChannels
%                 Nucleus_intensity{qq,nn} = [];
				Cluster_intensity{qq,nn} = [];
				OP_intensity{qq,nn} = [];
                OP_BG_intensity{qq,nn} = [];
			end
			
		end
		
	end
	
	Cluster_distCell{ff} = vertcat(Cluster_OP_dist{:});
	Cluster_volCell{ff} = vertcat(Cluster_volume{:});
	Cluster_solCell{ff} = vertcat(Cluster_solidity{:});
	Cluster_eloCell{ff} = vertcat(Cluster_elongation{:});
	Cluster_maskCell{ff} = vertcat(Cluster_mask_store{:});

% 	Nucleus_intCell{ff} = cell(1,numQuantChannels);
% 	for qq = 1:numQuantChannels
% 		Nucleus_intCell{ff}{qq} = vertcat(Nucleus_intensity{qq,:});
%     end
    
%     test = vertcat(Nucleus_intensity{qq,:});

	Cluster_intCell{ff} = cell(1,numQuantChannels);
	for qq = 1:numQuantChannels
		Cluster_intCell{ff}{qq} = vertcat(Cluster_intensity{qq,:});
    end

	OP_intCell{ff} = cell(1,numQuantChannels);
	for qq = 1:numQuantChannels
		OP_intCell{ff}{qq} = vertcat(OP_intensity{qq,:});
    end

    OP_BG_intCell{ff} = cell(1,numQuantChannels);
	for qq = 1:numQuantChannels
		OP_BG_intCell{ff}{qq} = vertcat(OP_BG_intensity{qq,:});
	end
	
end

save('variables_threshold_cluster2std_withOPBG_DB0.65.mat','-v7.3')

%% Load saved workspace
load('variables_threshold_cluster2std_withOPBG_DB0.65.mat') %_DB0.65

%% Sort into conditions


%Empty cells/arrays for sorted variables
sortedCondNames = cell(1,numConds);
sortedNumNuclei = zeros(1,numConds);
sortedNumFiles = zeros(1,numConds);
sortedDistCell = cell(1,numConds);
sortedVolCell = cell(1,numConds);
sortedSolCell = cell(1,numConds);
sortedEloCell = cell(1,numConds);
sortedMaskCell = cell(1,numConds);
sortedNucIntCell = cell(1,numQuantChannels);
sortedIntCell = cell(1,numQuantChannels);
sortedOPIntCell = cell(1,numQuantChannels);
for qq = 1:numQuantChannels
	%sortedNucIntCell{qq} = cell(1,numConds);
    sortedIntCell{qq} = cell(1,numConds);
	sortedOPIntCell{qq} = cell(1,numConds);
    sortedOPBGIntCell{qq} = cell(1,numConds);
end
sortedPixelSize = zeros(1,numConds);
sortedZStep = zeros(1,numConds);


%Actual sorting
for cc = 1:numConds
	
    %Condition names
	sortedCondNames{cc} = ...
		condNames(condInds==uniqueCondInds(cc));
	sortedCondNames{cc} = sortedCondNames{cc}{1};
	
    %Number of nuclei per condition
	sortedNumNuclei(cc) = ...
		sum(numNuclei_vec(condInds==uniqueCondInds(cc)));
	sortedNumFiles(cc) = sum(condInds==uniqueCondInds(cc));
	
    %Vertical concatenation of cluster and op properties of each condition
	Cluster_dists = vertcat(Cluster_distCell{condInds==uniqueCondInds(cc)});
	Cluster_vols = vertcat(Cluster_volCell{condInds==uniqueCondInds(cc)});
	Cluster_sols = vertcat(Cluster_solCell{condInds==uniqueCondInds(cc)});
	Cluster_elos = vertcat(Cluster_eloCell{condInds==uniqueCondInds(cc)});
	Cluster_masks = vertcat(Cluster_maskCell{condInds==uniqueCondInds(cc)});
% 	Nucleus_ints = Nucleus_intCell(condInds==uniqueCondInds(cc));
    Cluster_ints = Cluster_intCell(condInds==uniqueCondInds(cc));
	OP_ints = OP_intCell(condInds==uniqueCondInds(cc));
    OP_BG_ints = OP_BG_intCell(condInds==uniqueCondInds(cc));
    
    %Add concatenated properties to cell
	sortedDistCell{cc} = Cluster_dists;
	sortedVolCell{cc} = Cluster_vols;
	sortedSolCell{cc} = Cluster_sols;
	sortedEloCell{cc} = Cluster_elos;
	sortedMaskCell{cc} = Cluster_masks;

	%Pixel sizes and z-steps
    pixelSizes = pixelSize_array(condInds==uniqueCondInds(cc));
	zStepSizes = zStep_array(condInds==uniqueCondInds(cc));	
	sortedPixelSize(cc) = pixelSizes(1);
	sortedZStep(cc) = zStepSizes(1);
    
    %Loop through all quantification channels
	for qq = 1:numQuantChannels
		
        %Sort OP intensities
        sortedOPIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)OP_intCell{ind}{qq},....
			find(condInds==uniqueCondInds(cc)),...
			'UniformOutput',false));
		sortedOPIntCell{qq}{cc} = vertcat(sortedOPIntCell{qq}{cc}{:});
        
        %Sort OP BG intensities
        sortedOPBGIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)OP_BG_intCell{ind}{qq},....
			find(condInds==uniqueCondInds(cc)),...
			'UniformOutput',false));
		sortedOPBGIntCell{qq}{cc} = vertcat(sortedOPBGIntCell{qq}{cc}{:});
		
        %Sort cluster intensities
        sortedIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)Cluster_intCell{ind}{qq},....
			find(condInds==uniqueCondInds(cc)),...
			'UniformOutput',false));
		sortedIntCell{qq}{cc} = vertcat(sortedIntCell{qq}{cc}{:});

%         %Sort nucleus intensities
%         sortedNucIntCell{qq}{cc} = vertcat(arrayfun(...
% 			@(ind)Nucleus_intCell{ind}{qq},....
% 			find(condInds==uniqueCondInds(cc)),...
% 			'UniformOutput',false));
% 		sortedNucIntCell{qq}{cc} = vertcat(sortedNucIntCell{qq}{cc}{:});

	end
	
end


%% Principal Component based sorting of observations

% figure(1)
% clf
% 
% figure(2)
% clf

%Empty arrays/cells
sorted_central_slices = cell(1,numConds);
radius_median = zeros(1,numConds);
radius_CI = zeros(2,numConds);
in_range_perc = zeros(1,numConds);
radii_cell = cell(1,numConds);
crossCorr_vals = zeros(1,numConds);
crossCorr_CI = zeros(2,numConds);



%Loop through all conditions
for cc = 1:numConds
	
% 	figure(1)
	
	inclInds = ...
		sortedDistCell{cc}<=Inf & sortedVolCell{cc}>=Vol_threshold;
		
	dist_vals = [sortedDistCell{cc}(inclInds)];
    OP_S5P_vals = [sortedOPIntCell{1}{cc}(inclInds)];
	OP_S2P_vals = [sortedOPIntCell{2}{cc}(inclInds)];
    %Add for third channel 
	Cluster_S5P_vals = [sortedIntCell{1}{cc}(inclInds)];
	Cluster_S2P_vals = [sortedIntCell{2}{cc}(inclInds)];
    Vol_vals = [sortedVolCell{cc}(inclInds)];
	Elo_vals = [sortedEloCell{cc}(inclInds)];
	Sol_vals = [sortedSolCell{cc}(inclInds)];
	Central_slices = [sortedMaskCell{cc}(inclInds)];

    %OP background
    OP_BG_S5P_vals = [sortedOPBGIntCell{1}{cc}(inclInds)]; %(inclInds)
	OP_BG_S2P_vals = [sortedOPBGIntCell{2}{cc}(inclInds)];

    %Additional Nucleus intensities
%     Nuc_S5P_vals = [sortedNucIntCell{1}{cc}(randi(numel(sortedNucIntCell{1}{cc}),[10000,1]))];
% 	Nuc_S2P_vals = [sortedNucIntCell{2}{cc}(randi(numel(sortedNucIntCell{1}{cc}),[10000,1]))];

	frac_close = sum(dist_vals<=dist_threshold)./numel(dist_vals);
	in_range_perc(cc) = frac_close.*100;

	[ff,xx] = ksdensity(dist_vals,'BandWidth',0.07,'support','unbounded');
	
  % Plot Gene-cluster distance vs. probability. Optimize automatic plotting
  % routine. Automatically give number of rows and colums of subplot.

  %Also change auomatic adjustment of text placement and y-axis limits.

% 	%subplot(4,numConds,0.*numConds+cc)
% 	subplot(5,numConds/5,0.*numConds+cc)
%     plot(xx,ff,'k-','LineWidth',1)
% 	hold on
% 	if max(ff)>1
%         y_lim = get(gca, 'YLim');
%         plot([1,1].*dist_threshold,[0,ceil(y_lim(2) * 5)/5],'k-','LineWidth',1) %Plot threshold
%         text(2.3,ceil(y_lim(2) * 5)/5*0.95,sprintf('f(d< %d nm)=%.2f',dist_threshold*1000,in_range_perc(cc)))
%     else
%         plot([1,1].*dist_threshold,[0,1],'k-','LineWidth',1) %Plot threshold
%         text(2.3,0.9,sprintf('f(d< %d nm)=%.2f',dist_threshold*1000,in_range_perc(cc)))
%     end 
%     xlabel('Gene-cluster distance')
%     ylabel('Probability')
% 
%     title(strrep(sprintf('%s (n=%d)',...
% 		sortedCondNames{cc},numel(sortedDistCell{cc})), '_', ' '),'Interpreter','none')
% 	set(gca,'XLim',[-0.05,4.5])

    
%     %Plot Ser5p/Ser2P intensity with distance to cluster
% 	figure(2)
%     subplot(5,8,0.*numConds+cc)
%     %subplot(5,numConds/5,0.*numConds+cc)
% %     scatter3(OP_S5P_vals,OP_S2P_vals,dist_vals,10,dist_vals,'filled')
%     
% %     hold on
%     %S5P intensity cluster mask
%     [ff_O5,xx_O5] = ksdensity(OP_S5P_vals,'BandWidth',0.07,'support','unbounded');
%     %plot(xx,ff./5,'k-','LineWidth',1)
%     plot(xx_O5,ff_O5,'b-','LineWidth',1,'DisplayName','OP mask')
%     hold on
% %     pts = linspace(0,5,100);
% %     [ff_O5,xx_O5] = ksdensity(OP_S5P_vals,pts,'BandWidth',0.07,'support','unbounded');
% %     plot(xx_O5,ff_O5,'r--','LineWidth',1,'DisplayName','OP mask')
% %     hold on
% 
%     %S2P intensity cluster mask
%     [ff,xx] = ksdensity(OP_S2P_vals,'BandWidth',0.07,'support','unbounded');
%     %plot(ff./5,xx,'k-','LineWidth',1)
%     plot(ff,xx,'b-','LineWidth',1)
%     hold on
%     %S5P intensity OP BG nucleus 
%     %[ff_B5,xx_B5] = ksdensity(OP_BG_S5P_vals,pts,'BandWidth',0.07,'support','unbounded');
%     [ff_B5,xx_B5] = ksdensity(OP_BG_S5P_vals,'BandWidth',0.07,'support','unbounded');
%     %plot(xx,ff./5,'r-','LineWidth',1)
%     plot(xx_B5,ff_B5,'k-','LineWidth',1,'DisplayName','OP bg.')
%     hold on
%     %S2P intensity OP BG nucleus
%     [ff,xx] = ksdensity(OP_BG_S2P_vals,'BandWidth',0.07,'support','unbounded');
%     %plot(ff./5,xx,'r-','LineWidth',1)
%     plot(ff,xx,'k-','LineWidth',1)
% 
%     hold on 
%     
%     scatter(OP_S5P_vals,OP_S2P_vals,10,dist_vals,'filled','DisplayName',...
%         'Dist.')
%     colorbar; %comment, when scatterhist
%     caxis([0 2]); %comment, when scatterhist

    %Kolmogorov Smirnov Test
%     [h,p,k] = kstest2(OP_S2P_vals,OP_BG_S2P_vals)

    %Convolution for Ser5P
%     ff_c = conv(ff_B5,ff_O5)

%     %S5P intensity nucleus mask
%     [ff,xx] = ksdensity(Nuc_S5P_vals,'BandWidth',0.07,'support','unbounded');
%     %plot(xx,ff./5,'r-','LineWidth',1)
%     plot(xx,ff,'r-','LineWidth',1)
%     hold on
%     %S2P intensity nucleus mask
%     [ff,xx] = ksdensity(Nuc_S2P_vals,'BandWidth',0.07,'support','unbounded');
%     %plot(ff./5,xx,'r-','LineWidth',1)
%     plot(ff,xx,'r-','LineWidth',1)


%     %     scatterhist(OP_S5P_vals,OP_S2P_vals,'Kernel','on','Location','SouthWest','Direction','out');
%     xlabel('OP Ser5P (a.u.)');
%     ylabel('OP Ser2P (a.u.)');
%     hold all
%     title(strrep(sprintf('%s (n=%d)',...
% 		sortedCondNames{cc},numel(sortedDistCell{cc})), '_', ' '),'Interpreter','none')
% 	set(gca,'XLim',[0,5])
%     set(gca,'YLim',[0,2.5])
    
%     %Cluster solidity
%     figure(3)
%     subplot(5,8,0.*numConds+cc)
%     [ff,xx] = ksdensity(Sol_vals,'BandWidth',0.07,'support','unbounded');
%     plot(xx,ff,'b-','LineWidth',1)
%     xlabel('Cluster solidity')
%     ylabel('Probability')
%     title(strrep(sprintf('%s (n=%d)',...
% 		sortedCondNames{cc},numel(sortedDistCell{cc})), '_', ' '),'Interpreter','none')
% 	%set(gca,'XLim',[-0.05,4.5])
% 
% 
%     %Cluster size
%     figure(4)
%     subplot(5,8,0.*numConds+cc)
%     [ff,xx] = ksdensity(Vol_vals,'BandWidth',0.07,'support','unbounded');
%     plot(xx,ff,'b-','LineWidth',1)
%     xlabel('Cluster volume')
%     ylabel('Probability')
%     title(strrep(sprintf('%s (n=%d)',...
% 		sortedCondNames{cc},numel(sortedDistCell{cc})), '_', ' '),'Interpreter','none')

%     %Cluster size
%     figure(5)
%     titles = {'drll2';'foxd5';'gadd45ga';'iscub';'klf2b';'ripply1';'vamp2'};
%     if mod(cc,8) == 1 | mod(cc,8) == 2 | mod(cc,8) == 3 | mod(cc,8) == 4 | ...
%             mod(cc,8) == 5 | mod(cc,8) == 6 | mod(cc,8) == 7 
%         subplot(4,8,mod(cc,8))
%         title(titles(mod(cc,8))) 
%         scatter(ceil(cc/8),mean(dist_vals),'k','filled')
%         xlim([0 6])
%         xticks([1 2 3 4 5])
%         xticklabels({'Ob','Sp','Do','30%','50%'})
%         ylabel('Mean distance [µm]')
%         ylim([0.8 1.8])
%         hold on
%         
%         %S5P
%         subplot(4,8,8 + mod(cc,8))
%         scatter(ceil(cc/8),mean(OP_S5P_vals),'k','filled')
%         xlim([0 6])
%         xticks([1 2 3 4 5])
%         xticklabels({'Ob','Sp','Do','30%','50%'})
%         ylabel('Mean S5P int. (a.u.)')
%         ylim([0.9 2.0])
%         hold on
% 
%         %S2P
%         subplot(4,8,16 + mod(cc,8))
%         scatter(ceil(cc/8),mean(OP_S2P_vals),'k','filled')
%         xlim([0 6])
%         xticks([1 2 3 4 5])
%         xticklabels({'Ob','Sp','Do','30%','50%'})
%         ylabel('Mean S2P int. (a.u.)')
%         ylim([0.9 2.0])
%         hold on
% 
%         %Contact Prob
%         subplot(4,8,24 + mod(cc,8))
%         scatter(ceil(cc/8),in_range_perc(cc),'k','filled')
%         xlim([0 6])
%         xticks([1 2 3 4 5])
%         xticklabels({'Ob','Sp','Do','30%','50%'})
%         ylabel("Contact prob. f(d<" + dist_threshold*1000 + "nm)")
%         ylim([0 30])
%         hold on
% 
%     elseif mod(cc,8) == 0
%         %Distance
%         subplot(4,8,8)
%         title('zgc:64022') 
%         scatter(ceil(cc/8),mean(dist_vals),'k','filled')
%         xlim([0 6])
%         xticks([1 2 3 4 5])
%         xticklabels({'Ob','Sp','Do','30%','50%'})
%         ylabel('Mean distance [µm]')
%         ylim([0.8 1.8])
%         hold on
% 
%         %S5P
%         subplot(4,8,8 + 8)
%         scatter(ceil(cc/8),mean(OP_S5P_vals),'k','filled')
%         xlim([0 6])
%         xticks([1 2 3 4 5])
%         xticklabels({'Ob','Sp','Do','30%','50%'})
%         ylabel('Mean S5P int. (a.u.)')
%         ylim([0.8 1.8])
%         hold on
% 
%         %S2P
%         subplot(4,8,16 + 8)
%         scatter(ceil(cc/8),mean(OP_S2P_vals),'k','filled')
%         xlim([0 6])
%         xticks([1 2 3 4 5])
%         xticklabels({'Ob','Sp','Do','30%','50%'})
%         ylabel('Mean S2P int. (a.u.)')
%         ylim([0.8 1.8])
%         hold on
% 
%         %Contact Prob
%         subplot(4,8,24 + 8)
%         scatter(ceil(cc/8),in_range_perc(cc),'k','filled')
%         xlim([0 6])
%         xticks([1 2 3 4 5])
%         xticklabels({'Ob','Sp','Do','30%','50%'})
%         ylabel("Contact prob. f(d<" + dist_threshold*1000 + "nm)")
%         ylim([0 30])
%         hold on

    end

%     legend([],'FontSize',5)
% 
%     sprintf('Nuc Median S5P %d',median(Nuc_S5P_vals))
%     sprintf('Nuc Median S2P %d',median(Nuc_S2P_vals))
    
		
	
% 	% PCA, input: Rows of X are observations, columns the variables
% 	
%     %Create observation matrix with all measured properties
%     observationMatrix = [dist_vals,...
% 		OP_S5P_vals,OP_S2P_vals, ...
% 		Cluster_S5P_vals,Cluster_S2P_vals,...
% 		Vol_vals,Elo_vals,Sol_vals]; 
% 	
%     %Perform PCA with three components
% 	[PCA_coeffs,PCA_scores,~,~,PCA_percExplained] = ...
% 		pca(observationMatrix,...
% 		'NumComponents',3);
% 	
%     %Get maximal values
% 	[maxVal,maxInd] = max(abs(PCA_coeffs),[],2);	
% 	[~,maxIndCluster] = max(PCA_coeffs(6,:));
% 	
% 	maxValOP = max(abs(PCA_coeffs(1:3,:)),[],1);
% 	[~,maxIndOP] = max(maxValOP);
% 	maxIndOP = maxIndOP + (maxIndOP==maxIndCluster);
% 	
%     %Set order for PCA (according to max vol cluster)
% 	PCA_order = [maxIndCluster,maxIndOP,...
% 		setdiff([1,2,3],[maxIndCluster,maxIndOP])];
% 	if PCA_order(1) == PCA_order(2)
% 		PCA_order = PCA_order(2:end);
% 	end
% 	PCA_coeffs_plot = PCA_coeffs(:,PCA_order);
% 	
%     figure(3);
% 	subplot(4,numConds,1.*numConds+cc)
% 	imagesc(PCA_coeffs_plot',[-1,+1])
% % 	
% 
% 	PCA_labels = arrayfun(@(nn) ...
% 		sprintf('PC %d (%1.1f%%)',...
% 		PCA_order(nn),PCA_percExplained(PCA_order(nn))),1:3,...
% 		'UniformOutput',false);
% 
% % 	set(gca,'YTick',1:3,'YTickLabel',PCA_labels)
% % 	
% % 	title(sprintf('f(d<%d nm)=%1.1f%%',...
% % 		dist_threshold.*1000,frac_close.*100),...
% % 		'FontWeight','normal')
% 
% 	
%     %Define colormap matrix
% 	colormap_matrix = ones(65,3);
% 	colormap_matrix(1:33,2) = linspace(0,1,33);
% 	colormap_matrix(1:33,3) = linspace(0,1,33);
% 	colormap_matrix(33:65,1) = linspace(1,0,33);
% 	colormap_matrix(33:65,2) = linspace(1,0,33);
% 	colormap(colormap_matrix);
% 	
%     %Order for first two components
% 	PCA_coeffs = PCA_coeffs(:,PCA_order([1,2]));
% 	PCA_scores = PCA_scores(:,PCA_order([1,2]));
% 	PCA_percExplained = PCA_percExplained(PCA_order([1,2]));
% 	PCA_labels = PCA_labels([1,2]);
% 	
% 	topCount = 10;
% 	[~,S5P_sortInds] = sort(OP_S5P_vals,'descend');
% 	S5P_topInds = S5P_sortInds(1:topCount);
% 	[~,Vol_sortInds] = sort(Vol_vals,'descend');
% 	Vol_topInds = Vol_sortInds(1:topCount);
% 	
% 	%Vector for S5P intensity
%     S5P_vec = [...
% 		mean(PCA_scores(S5P_topInds,1)),...
% 		mean(PCA_scores(S5P_topInds,2))];
% 	
%     %Vector for cluster volume
%     Vol_vec = [...
% 		mean(PCA_scores(Vol_topInds,1)),...
% 		mean(PCA_scores(Vol_topInds,2))];
% 	
% % 	subplot(4,numConds,2.*numConds+cc)
% % 	plot(PCA_scores(:,1),PCA_scores(:,2),'k.',...
% % 		'MarkerEdgeColor',[0.6,0.6,0.6])
% % 	hold on
% % 	plot([0,S5P_vec(1)],[0,S5P_vec(2)],'m-','LineWidth',1)
% % 	plot([0,Vol_vec(1)],[0,Vol_vec(2)],'b-','LineWidth',1)
% % 	xlabel(PCA_labels{1})
% % 	ylabel(PCA_labels{2})
% % 	axis equal
% % % 	set(gca,'XLim',[-6,3],'YLim',[-2,8])
% 
% 	if S5P_vec(2)>0
% 		Vol_ortho_vec = Vol_vec*[0,-1;+1,0];
% 	else
% 		Vol_ortho_vec = Vol_vec*[0,+1;-1,0];
% 	end
% 		
% % 	plot([0,Vol_ortho_vec(1)],[0,Vol_ortho_vec(2)],'k-','LineWidth',1)
% 	
% 
%     %Perform transformation, such that volume points to north
% 	
%     %Define unit vectors
%     unit_vec_2 = Vol_vec./norm(Vol_vec);
% 	unit_vec_1 = Vol_ortho_vec./norm(Vol_ortho_vec);
%     
%     %Define transformation matrix
% 	Trafo_matrix = [unit_vec_1',unit_vec_2'];
% 	
%     %Perform transformation
% 	trafo_scores = PCA_scores*Trafo_matrix;
% 	trafo_Vol_vec = Vol_vec*Trafo_matrix;
% 	trafo_S5P_vec = S5P_vec*Trafo_matrix;	
% 	
% 	
% % 	subplot(4,numConds,3.*numConds+cc)
% % 	plot(trafo_scores(:,1),trafo_scores(:,2),'k.',...
% % 		'MarkerEdgeColor',[0.6,0.6,0.6])
% % 	hold on
% % 	plot([0,trafo_S5P_vec(1)],[0,trafo_S5P_vec(2)],'m-','LineWidth',1)
% % 	plot([0,trafo_Vol_vec(1)],[0,trafo_Vol_vec(2)],'b-','LineWidth',1)
% % 	xlabel('Transformed 1')
% % 	ylabel('Transformed 2')
% % 	axis equal
% % 	set(gca,'XLim',[-6,3],'YLim',[-2,8])
% 	
% 
%     %Get angels from x,y coordinate
% 	angles = -atan2(trafo_scores(:,2),trafo_scores(:,1))./2./pi; % in radians
% 	%apply modulo operator
%     angles = mod(angles,1);
% 	%Perform angle shift
%     angles = mod(angles+angle_shift,1);
% 	%Get number of angles
%     numPoints = numel(angles);
% 	%Define pseudo-time coordinate
%     coord_s = ((1:numPoints)-1)./numPoints;
% 	
% 	[angles,sortInds] = sort(angles);
% 	
%     %Calculate moving mean of Windowsize 20
% 	windowSize = 20;
% 	curveSmooth = @(xx) movmean(...
% 		padarray(xx,windowSize,'circular','both'),...
% 		windowSize,'Endpoints','discard');
% % 	plot(curveSmooth(trafo_scores(sortInds,1)),...
% % 		curveSmooth(trafo_scores(sortInds,2)),'k-','LineWidth',1)
% 	
%     %Get radius of PCA
% 	radii = sqrt(sum(trafo_scores.^2,2));
% 	%Median radius
%     radius_median(cc) = median(radii);
% 	%Radius confidence interval
%     radius_CI(:,cc) = bootci(1000,@median,radii);
% 	radii_cell{cc} = radii;
% 	
%     %Sort all variables
% 	dist_vals = dist_vals(sortInds);
% 	OP_S5P_vals = OP_S5P_vals(sortInds);
% 	OP_S2P_vals = OP_S2P_vals(sortInds);
% 	Cluster_S5P_vals = Cluster_S5P_vals(sortInds);
% 	Cluster_S2P_vals = Cluster_S2P_vals(sortInds);
% 	Vol_vals = Vol_vals(sortInds);
% 	Elo_vals = Elo_vals(sortInds);
% 	Sol_vals = Sol_vals(sortInds);
% 	sorted_central_slices{cc} = Central_slices(sortInds);
% 	
%     %Calculate cross correlation
% 	crossCorr_vals(cc) = corr(...
% 		OP_S5P_vals,circshift(Elo_vals,-ceil(numPoints.*0.2)));
% 	
%     %Bootstrap confidence interval or cross correlation
%     crossCorr_CI(:,cc) = bootci(1000,@corr,...
% 		OP_S5P_vals,circshift(Elo_vals,-ceil(numPoints.*0.2)));
% 	
%     crossCorr_vals(cc) = corr(...
% 		OP_S5P_vals,circshift(Cluster_S2P_vals,-ceil(numPoints.*0.25)));
% 	crossCorr_CI(:,cc) = bootci(1000,@corr,...
% 		OP_S5P_vals,circshift(Cluster_S2P_vals,-ceil(numPoints.*0.25)));
% 	
% 	title(sprintf('\\rho=%2.2f [%2.2f,%2.2f]',crossCorr_vals(cc),...
% 		crossCorr_CI(1,cc),crossCorr_CI(2,cc)),...
% 		'FontWeight','normal')
% 
% 	
% 	figure(2)
% 	
% 	%Bin discretization for plotting of properties vs. pseudo time s
% 	
%     %Define bin boarders
% 	windowWidth = 0.2;
% 	numWindows = 100;
% 	windowCenters = linspace(0,1,numWindows);
% 	leftEdges = windowCenters-windowWidth./2;
% 	rightEdges = windowCenters+windowWidth./2;
% 
% 	%Empty arrays/cells for each plot
% 	mean_dist = zeros(1,numWindows);
% 	mean_OP_S5P = zeros(1,numWindows);
% 	mean_OP_S2P = zeros(1,numWindows);
% 	mean_Cluster_S5P = zeros(1,numWindows);
% 	mean_Cluster_S2P = zeros(1,numWindows);
% 	mean_Vol = zeros(1,numWindows);
% 	mean_Elo = zeros(1,numWindows);
% 	mean_Sol = zeros(1,numWindows);
%     
%     %Loop through all windows
% 	for nn = 1:numWindows
% 		
%         %Get window edges
% 		thisLeftEdge = leftEdges(nn);
% 		thisRightEdge = rightEdges(nn);
% 		
%         %Get window indices
% 		windowInds = find(coord_s>=thisLeftEdge & coord_s<thisRightEdge);
% 		if thisRightEdge>1
% 			windowInds = [windowInds,find(coord_s<(thisRightEdge-1))];
% 		end
% 		if thisLeftEdge<0
% 			windowInds = [windowInds,find(coord_s>(thisLeftEdge+1))];
%         end
%         
%         %Get mean properties of each window
% 		mean_dist(nn) = median(dist_vals(windowInds));
% 		mean_OP_S5P(nn) = median(OP_S5P_vals(windowInds));
% 		mean_OP_S2P(nn) = median(OP_S2P_vals(windowInds));
% 		mean_Cluster_S5P(nn) = median(Cluster_S5P_vals(windowInds));
% 		mean_Cluster_S2P(nn) = median(Cluster_S2P_vals(windowInds));
% 		mean_Vol(nn) = median(Vol_vals(windowInds));
% 		mean_Elo(nn) = median(Elo_vals(windowInds));
% 		mean_Sol(nn) = median(Sol_vals(windowInds));
% 	end
	
%     %Plot all properties within subplot
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
	
% 	% Cross-correlation analysis of S5P, S2P and Solidity vs. Elongation
	
%     %Define details for cross-correlation
% 	shiftRange = numPoints;
% 	shiftVals = (1:5:shiftRange)-1-ceil(numPoints./2);
% 	numShifts = numel(shiftVals);
% 	padRange = numPoints;
% 	n_boot = 10;
% 
% 	%Pad the arrays 
% 	S5P_vals = padarray(Cluster_S5P_vals,padRange,'circular','both');
% 	S2P_vals = padarray(Cluster_S2P_vals,padRange,'circular','both');
% 	Sol_vals = padarray(Sol_vals,padRange,'circular','both');
% 	Elo_vals = padarray(Elo_vals,padRange,'circular','both');
% 	
%     %Empty arrays for correlation results
% 	corr_S5P = zeros(1,numShifts);
% 	corr_S5P_CI = zeros(2,numShifts);
% 	corr_S2P = zeros(1,numShifts);
% 	corr_S2P_CI = zeros(2,numShifts);
% 	corr_Sol = zeros(1,numShifts);
% 	corr_Sol_CI = zeros(2,numShifts);
% 	shiftDist = zeros(1,numShifts);
% 
%     %Loop through all shifts and calculate correlation with CI
% 	for nn = 1:numShifts
% 		
%         %Get shift index
% 		shiftInds = numPoints+(1:numPoints)+shiftVals(nn);
% 		
%         %Actual correlation calculation
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
% 
%     %Plot moving mean on top of cross correlation plot (
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
	
% end







% figure(3)
% clf
% 
% errorbar(in_range_perc,crossCorr_vals,...
% 	+crossCorr_CI(2,:)-crossCorr_vals,...
% 	-crossCorr_CI(1,:)+crossCorr_vals,'ko');
% xlabel(sprintf('f(d<%d nm) [percent]',1000.*dist_threshold))
% ylabel('1/4 shift correlation')
% set(gca,'YLim',[-1,1])


%% Overview plot of changes with developmental stages

figure(1)
clf

figure(2)
clf

figure(3)
clf

figure(4)
clf

figure(10)
clf

% plotSets = {...
% 	{'Oblong_foxd5','Sphere_foxd5','Dome_foxd5','Epi30_foxd5','Epi50_foxd5'},...
% 	{'Oblong_ripply1','Sphere_ripply1','Dome_ripply1','Epi30_ripply1','Epi50_ripply1'},...
% 	{'Oblong_vamp2','Sphere_vamp2','Dome_vamp2','Epi30_vamp2','Epi50_vamp2'},...
% 	{'Oblong_drll2','Sphere_drll2','Dome_drll2','Epi30_drll2','Epi50_drll2'},...
% 	{'Oblong_klf2b','Sphere_klf2b','Dome_klf2b','Epi30_klf2b','Epi50_klf2b'},...
% 	{'Oblong_zgc:64022','Sphere_zgc:64022','Dome_zgc:64022','Epi30_zgc:64022','Epi50_zgc:64022'},...
% 	{'Oblong_gadd45ga','Sphere_gadd45ga','Dome_gadd45ga','Epi30_gadd45ga','Epi50_gadd45ga'},...
% 	{'Oblong_iscub','Sphere_iscub','Dome_iscub','Epi30_iscub','Epi50_iscub'},...
% 	};
% setNames = {'foxd5','ripply1','vamp2','drll.2',...
%     'klf2b','zgc:64022','gadd45ga','iscub'};

plotSets = {...
	{'Oblong_zgc:64022','Sphere_zgc:64022','Dome_zgc:64022','Epi30_zgc:64022','Epi50_zgc:64022'},...
	};
setNames = {'zgc:64022'};



stageNames = {'Obl','Sph','Dome','Epi30','Epi50'};
setColors = {[0,0,1],[1,0,0],[0.5,0.5,0.5]};
numSets = numel(plotSets);

sorted_central_slices = cell(1,numConds);

%For interpolation of data
Contact_mean_all = [];
S5P_mean_all = [];
S2P_mean_all = [];
%End for interpolation



for nn_set = 1:numSets
	
	condNamesInSet = plotSets{nn_set};
	
	numCondsInSet = numel(condNamesInSet);
	
	in_range_perc = zeros(1,numCondsInSet);
	in_range_perc_CI = zeros(2,numCondsInSet);
	dist_median = zeros(1,numCondsInSet);
    dist_mean = zeros(1,numCondsInSet);
	dist_CI = zeros(2,numCondsInSet);
	S5P_median = zeros(1,numCondsInSet);
    S5P_mean = zeros(1,numCondsInSet);
	S5P_CI = zeros(2,numCondsInSet);
    S5P_CI_mean = zeros(2,numCondsInSet);
	S2P_median = zeros(1,numCondsInSet);
    S2P_mean = zeros(1,numCondsInSet);
	S2P_CI = zeros(2,numCondsInSet);
    S2P_CI_mean = zeros(2,numCondsInSet);
    S2P_top = zeros(1,numCondsInSet);
    S2P_top_CI = zeros(2,numCondsInSet);
    S2P_perc = zeros(1,numCondsInSet);
    S2P_perc_CI = zeros(2,numCondsInSet);
    S5P_perc = zeros(1,numCondsInSet);
    S5P_perc_CI = zeros(2,numCondsInSet);

    %setCrossCorr = zeros(1,numCondsInSet);
    %setCrossCorr_CI = zeros(2,numCondsInSet);

	dist_value_array = [];
	S5P_value_array = [];
	S2P_value_array = [];
	grouping_array = [];

	
	for cc = 1:numCondsInSet
		
		thisCondInd = find(...
			cellfun(@(elmt)strcmp(elmt,plotSets{nn_set}{cc}),sortedCondNames));
		
        %setCrossCorr(cc) = crossCorr_vals(thisCondInd);
        %setCrossCorr_CI(:,cc) = crossCorr_CI(:,thisCondInd);

        top_prctile = 5;
% 		dist_threshold = 0.3;
% 		Vol_threshold = 0.02;0.4;0.02;
		inclInds = ...
			sortedDistCell{thisCondInd}<=Inf ...
			& sortedVolCell{thisCondInd}>=Vol_threshold;
		
		dist_vals = [sortedDistCell{thisCondInd}(inclInds)];
		OP_S5P_vals = [sortedOPIntCell{1}{thisCondInd}(inclInds)];
		OP_S2P_vals = [sortedOPIntCell{2}{thisCondInd}(inclInds)];
        OP_OP_vals = [sortedOPIntCell{3}{thisCondInd}(inclInds)];
		Cluster_S5P_vals = [sortedIntCell{1}{thisCondInd}(inclInds)];
		Cluster_S2P_vals = [sortedIntCell{2}{thisCondInd}(inclInds)];
        Cluster_OP_vals = [sortedIntCell{3}{thisCondInd}(inclInds)];
		Vol_vals = [sortedVolCell{thisCondInd}(inclInds)];
		Elo_vals = [sortedEloCell{thisCondInd}(inclInds)];
		Sol_vals = [sortedSolCell{thisCondInd}(inclInds)];
		
        
        %figure(10)
        subplot(numCondsInSet,numSets,(cc-1)*numSets + nn_set)
%         scatter(OP_S5P_vals,OP_OP_vals,'filled','MarkerFaceColor',[.5 .5 .5],...
%             'MarkerEdgeColor',[.5 .5 .5])
        scatter(OP_S5P_vals(dist_vals>0.25),OP_OP_vals(dist_vals>0.25), ...
            'filled','MarkerFaceColor',[.5 .5 .5],...
            'MarkerEdgeColor',[.5 .5 .5])
        hold on
%         scatter(Cluster_S5P_vals,Cluster_OP_vals)
        scatter(Cluster_S5P_vals(dist_vals>0.25),Cluster_OP_vals(dist_vals>0.25))
        xlabel('S5P int. (a.u.)')
        ylabel('OP int. (a.u.)')
        str = sprintf('Gene %s Stage %s', setNames{nn_set}, stageNames{cc});
        title(str)
        hold on;
        r = corrcoef(Cluster_S5P_vals,Cluster_OP_vals);
        Fit = polyfit(Cluster_S5P_vals,Cluster_OP_vals,1);
        f = polyval(Fit,Cluster_S5P_vals);
        plot(Cluster_S5P_vals,f)
        legend('OP mask', 'Cluster mask', 'Fit')
        hold on;
        xlim([0.5,8])
        ylim([0.9,1.7])
        sgtitle('No Contact')
        %waitforbuttonpress
        hold off;
        %clf
        
        % 	Central_slices = [sortedCentralSliceCell{thisCondInd}(inclInds)];
		
%         figure(2)
%         
%         subplot(1,numSets,nn_set)
%         
% %         plot3(dist_vals,OP_S5P_vals,OP_S2P_vals,'k.')
% %         hold on
% %         set(gca,'XLim',[0,2],'YLim',[0,7],'ZLim',[0,3])
% %         ylabel('Pol II Ser5P')
% %         xlabel('Distance [\mum]')
% %         zlabel('Pol II Ser2P')
%         
%         plot(dist_vals,OP_S5P_vals,'k.')
%         hold on
%         
%         set(gca,'XLim',[0,5],'YLim',[0,6])
%         
%         
%         
%         ylabel('Pol II Ser5P')
%         xlabel('Distance [\mum]')
%         
%         title(setNames{nn_set})

        
		n_boot = 1000;
		dist_median(cc) = median(dist_vals);
        dist_mean(cc) = mean(dist_vals);
		dist_CI(:,cc) = bootci(n_boot,@median,dist_vals);
		S5P_median(cc) = median(OP_S5P_vals);
        S5P_mean(cc) = mean(OP_S5P_vals);
		S5P_CI(:,cc) = bootci(n_boot,@median,OP_S5P_vals);
        S5P_CI_mean(:,cc) = bootci(n_boot,@mean,OP_S5P_vals);
		S2P_median(cc) = median(OP_S2P_vals);
        S2P_mean(cc) = mean(OP_S2P_vals);
		S2P_CI(:,cc) = bootci(n_boot,@median,OP_S2P_vals);
        S2P_CI_mean(:,cc) = bootci(n_boot,@mean,OP_S2P_vals);
		
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
    
    %For Interpolation
    Contact_mean_all = [Contact_mean_all, in_range_perc];
    S5P_mean_all = [S5P_mean_all,S5P_mean];
	S2P_mean_all = [S2P_mean_all,S2P_mean];
    %End for interpolation


    %%% Main figure start here%%%
	figure(1)
	
    subplot(3,numSets,0.*numSets+nn_set)

% 	errorbar(1:numCondsInSet,dist_median,...
% 		dist_CI(1,:)-dist_median,...
% 		dist_median-dist_CI(2,:),...
% 		'k-o')

%     errorbar(1:numCondsInSet,dist_mean,...
% 		dist_CI(1,:)-dist_mean,...
% 		dist_mean-dist_CI(2,:),...
% 		'k-o')

%     hold on
% %     scatter(1:numCondsInSet,dist_mean,'red','filled')
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
		'YLim',[0.0,40])

    title(setNames{nn_set})

    
	subplot(3,numSets,numSets+nn_set)

% 	errorbar(1:numCondsInSet,S5P_median,...
% 		S5P_CI(1,:)-S5P_median,...
% 		S5P_median-S5P_CI(2,:),...
% 		'k-o')

    errorbar(1:numCondsInSet,S5P_mean,...
		S5P_CI_mean(1,:)-S5P_mean,...
		S5P_mean-S5P_CI_mean(2,:),...
		'k-o')

    hold on 
%     scatter(1:numCondsInSet,S5P_mean,'red','filled')
%     errorbar(1:numCondsInSet,S5P_perc,...
% 		S5P_perc_CI(1,:)-S5P_perc,...
% 		S5P_perc-S5P_perc_CI(2,:),...
% 		'k-o')
	if nn_set == 1
        ylabel('Mean S5P int. (a.u.)')
%         ylabel('Median S5P int. (a.u.)')
%         ylabel('S5P %')
    else
        set(gca,'YTickLabel',[])
    end
    
	set(gca,'XTick',1:numCondsInSet,'XTickLabel',stageNames,...
		'XLim',([1,numCondsInSet])+[-0.5,+0.5],...
		'YLim',[0.9,2.0])
	

	subplot(3,numSets,2.*numSets+nn_set) %use 4 in case of 4 rows
% 	errorbar(1:numCondsInSet,S2P_median,...
% 		S2P_CI(1,:)-S2P_median,...
% 		S2P_median-S2P_CI(2,:),...
% 		'k-o')
% % 
    errorbar(1:numCondsInSet,S2P_mean,...
		S2P_CI_mean(1,:)-S2P_mean,...
		S2P_mean-S2P_CI_mean(2,:),...
		'k-o')

    hold on 
%     scatter(1:numCondsInSet,S2P_mean,'red','filled')
%     errorbar(1:numCondsInSet,S2P_perc,...
% 		S2P_perc_CI(1,:)-S2P_perc,...
% 		S2P_perc-S2P_perc_CI(2,:),...
% 		'k-o')
    if nn_set == 1
        ylabel('Mean S2P int. (a.u.)')
%         ylabel('Median S2P int. (a.u.)')
%         ylabel('S2P %')
    else
        set(gca,'YTickLabel',[])
    end
	
	set(gca,'XTick',1:numCondsInSet,'XTickLabel',stageNames,...
		'XLim',([1,numCondsInSet])+[-0.5,+0.5],...
		'YLim',[0.9,1.6])

    
%     subplot(4,numSets,3.*numSets+nn_set)
% % 	errorbar(1:numCondsInSet,S2P_median,...
% % 		S2P_CI(1,:)-S2P_median,...
% % 		S2P_median-S2P_CI(2,:),...
% % 		'k-o')
% %     errorbar(1:numCondsInSet,setCrossCorr,...
% % 		setCrossCorr_CI(1,:)-setCrossCorr,...
% % 		setCrossCorr-setCrossCorr_CI(2,:),...
% % 		'k-o')
%     plot(S5P_mean,S2P_mean)
%     hold on 
%     scatter(S5P_mean(1),S2P_mean(1),150,'s','k','LineWidth',1.0)
%     hold on 
%     scatter(S5P_mean,S2P_mean,36,in_range_perc,'filled')
%     colorbar; %comment
%     caxis([0 30]);
%     
%     if nn_set == 1
%         ylabel('Mean S2P int. (a.u.)')
% %         ylabel('Clock radius')
%     else
%         set(gca,'YTickLabel',[])
%     end
% 	
%     xlabel('Mean S5P int. (a.u.)')
% 
% 	set(gca,'XLim',[0.8,2],...
% 		'YLim',[0.8,1.4])
	figure(2)
    h(nn_set) = plot(S5P_mean,S2P_mean,'DisplayName',setNames{nn_set})
    hold on 
    scatter(S5P_mean(1),S2P_mean(1),150,'s','k','LineWidth',1.0)
    hold on 
    scatter(S5P_mean,S2P_mean,36,in_range_perc,'filled')
    c = colorbar; %comment
    c.Label.String='Contact %'
    caxis([0 30]);

    ylabel('Mean S2P int. (a.u.)')
    xlabel('Mean S5P int. (a.u.)')

	set(gca,'XLim',[0.9,2.0],...
		'YLim',[0.9,1.3])

    %%%Main figures end here%%%
    

    figure(3)
    h(nn_set) = plot(S5P_mean,S2P_mean,'DisplayName',setNames{nn_set})
    hold on 
    scatter(S5P_mean(1),S2P_mean(1),150,'s','k','LineWidth',1.0)
    hold on 
    scatter(S5P_mean,S2P_mean,36,dist_mean,'filled')
    c = colorbar; %comment
    c.Label.String='Mean gene-cluster dist. [\mum] '
%     caxis([0 30]);

    ylabel('Mean S2P int. (a.u.)')
    xlabel('Mean S5P int. (a.u.)')

	set(gca,'XLim',[0.8,2.0],...
		'YLim',[0.8,1.4])
    legend(h)
    
    %Filter and plot zgc data only
    if nn_set == 6
        figure(4)
        colororder({'k','k'})
        
        yyaxis left

        errorbar(1:numCondsInSet,S5P_mean,...
		S5P_CI_mean(1,:)-S5P_mean,...
		S5P_mean-S5P_CI_mean(2,:),...
		'r-o', 'LineWidth',1.0)

        hold on

        errorbar(1:numCondsInSet,S2P_mean,...
		S2P_CI_mean(1,:)-S2P_mean,...
		S2P_mean-S2P_CI_mean(2,:),...
		'-o', 'color', [.5 .5 .5], ...
        'LineWidth',1.0)

        set(gca,'XTick',1:numCondsInSet,'XTickLabel',stageNames,...
		'XLim',([1,numCondsInSet])+[-0.5,+0.5],...
		'YLim',[0.9,2.0])

        ylabel('Mean Pol II S5P and S2P int. (a.u.)')

        yyaxis right

        errorbar(1:numCondsInSet,in_range_perc,...
		    in_range_perc_CI(1,:)-in_range_perc,...
		    in_range_perc-in_range_perc_CI(2,:),...
		    'k-o', 'LineWidth',1.0)
        ylabel('Gene-cluster contact %')
	    set(gca,'XTick',1:numCondsInSet,'XTickLabel',stageNames,...
		    'XLim',([1,numCondsInSet])+[-0.5,+0.5],...
		    'YLim',[0.0,40])
       
        hold on 

        title(setNames{nn_set})
        legend('Pol II S5P', 'Pol II S2P', 'Contact %')
    end

end
figure(2)
legend('zgc:64022','','')

%Save figures
figure(2)

% savefig('Combined_Surface_Contact_MeanS5P_MeanS2P.fig')
% saveas(gcf,'Combined_Surface_Contact_MeanS5P_MeanS2P.png')

figure(4)
% savefig('ZGC_Quantification_Contact_MeanS5P_MeanS2P.fig')
% saveas(gcf,'ZGC_Quantification_Contact_MeanS5P_MeanS2P.png')

% %% -- Interpolation figure
% 
% figure(3)
% clf
% 
% grid_N = 150; %150
% grid_vals = zeros(grid_N,grid_N);
% averaging_N = numel(S2P_mean_all); %8
% averaging_KK = 0.4;
% 
% %Min and max values obtained from figure 2
% S2P_min = 0.9; %0.95;
% %S2P_min = floor(min(S2P_mean_all))
% S2P_max = 1.3; %1.3;
% S5P_min = 0.9;
% S5P_max = 1.9;
% S2P_vec = linspace(S2P_min,S2P_max,grid_N);
% S5P_vec = linspace(S5P_min,S5P_max,grid_N);
% 
% for mm = 1:grid_N
%     for nn = 1:grid_N
% 
%         S2P_val = S2P_vec(mm);
%         S5P_val = S5P_vec(nn);
% 
%         grid_dist_vec = ...
%             + (S2P_mean_all-S2P_val).^2./var(S2P_mean_all) ...
%             + (S5P_mean_all-S5P_val).^2./var(S5P_mean_all);
% 
%         %[~,sort_inds] = sort(grid_dist_vec,'ascend');
%         %grid_vals(mm,nn) = mean(Contact_mean_all(sort_inds(1:averaging_N)));
% 
%         weights = averaging_KK./(averaging_KK+grid_dist_vec);
%         grid_vals(mm,nn) = sum(weights.*Contact_mean_all)./sum(weights);
% 
%     end
% end
% 
% cla
% imagesc(S5P_vec,S2P_vec,grid_vals)
% set(gca,'YDir','normal')
% colorbar
% hold on
% 
% 
% % --- clustering
% 
% S2P_mean_zScore = (S2P_mean_all-mean(S2P_mean_all))...
%     ./std(S2P_mean_all);
% S5P_mean_zScore = (S5P_mean_all-mean(S5P_mean_all))...
%     ./std(S5P_mean_all);
% 
% %observMatrix = [S2P_mean_all',S5P_mean_all',Contact_mean_all'];
% observMatrix = [S2P_mean_zScore',S5P_mean_zScore'];
% ClustEval = evalclusters(observMatrix,"kmeans","gap","KList",1:10);
% 
% numClust = ClustEval.OptimalK;
% 
% [clustIdx,clusterCentr,sumdist] = ...
%     kmeans(observMatrix,numClust);
% 
% %markerStyles = {'ko','ks','kd'};
% markerStyles = {'o','o','o'};
% 
% for cc = 1:numClust
% 
%     S2P_clust = S2P_mean_all(cc==clustIdx);
%     S5P_clust = S5P_mean_all(cc==clustIdx);
%     
%     hullInds = convhull(S5P_clust,S2P_clust);
% 
%     %     plot(S5P_clust(hullInds),S2P_clust(hullInds),'k-',...
% %         'LineWidth',1)
% 
%     plot(S5P_clust(hullInds),S2P_clust(hullInds),'k-',...
%          'LineWidth',1)
% 
%     plot(S5P_clust,S2P_clust,...
%         markerStyles{cc},'MarkerSize',8,'MarkerEdgeColor',[0,0,0],...
%         'MarkerFaceColor',[1,1,1], 'LineWidth',1.0)
% 
% %     plot(S5P_clust,S2P_clust,...
% %         markerStyles{cc},'MarkerSize',8,'MarkerEdgeColor',[1,1,1],...
% %         'MarkerFaceColor',[1,1,1])
% %     plot(S5P_clust,S2P_clust,...
% %         markerStyles{cc},'MarkerSize',8,'MarkerEdgeColor',[0,0,0],...
% %         'LineWidth',1.0)
% 
% 
% 
% 
% end
% 
% %title('Contact%')
% xlabel('Mean S5P Int. (a.u.)')
% ylabel('Mean S2P Int. (a.u.)')
% c = colorbar;
% c.Label.String='Contact %'
% %caxis([0,16])
% %colormap(flipud(parula))
% colormap(parula)
% set(gca,'Box','on')
% set(gca,'XLim',[0.9,1.9],...
% 		'YLim',[0.9,1.3], ...
%         'Box','on')
%     
% 
% figure(3)
% % savefig('Surface_Interpolation_Contact_MeanS5P_MeanS2P.fig')
% 
