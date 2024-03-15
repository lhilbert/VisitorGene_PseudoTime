clear all
	
%% Plot the example images

SourceFileCell = {...
	'./ExtractedStacks_Sphere/Cond_5/Image_6.mat',...
	};

OP_blurRange = 0.03; % micrometers
Ser2P_blurRange = 0.01; % micrometers
Ser5P_blurRange = 0.01; % micrometers

imgRanges = {...
	[30,30,66.5,66.5]+[-6,+6,-6,+6]};
zCoordinate = [22];
rotate_flag = [false];

plotTitles = {'foxd5'};

scaleBar = 4.0; % in microns

numPlots = numel(imgRanges);

figure(1)
clf

for pp = 1:numPlots
	
	pp
	
	thisFilePath = SourceFileCell{pp};
	loadStruct = load(thisFilePath,...
		'imgStack','imgSize','pixelSize','zStepSize');
	imgStack = loadStruct.imgStack;
	imgSize = loadStruct.imgSize;
	pixelSize = loadStruct.pixelSize;
	zStepSize = loadStruct.zStepSize;

	thisImg = imgStack;
	thisPixelSize = pixelSize;
	
	% --- Removal of DNA background
	OP_img = double(thisImg{1}(:,:,zCoordinate(pp)));
	Ser2P_img = double(thisImg{2}(:,:,zCoordinate(pp)));
	Ser5P_img = double(thisImg{3}(:,:,zCoordinate(pp)));
	
	if rotate_flag(pp)
		OP_img = OP_img';
		Ser2P_img = Ser2P_img';
		Ser5P_img = Ser5P_img';
	end

	OP_img = imgaussfilt(OP_img,...
		OP_blurRange./thisPixelSize);
	Ser5P_img = imgaussfilt(Ser5P_img,...
		Ser5P_blurRange./thisPixelSize);
	Ser2P_img = imgaussfilt(Ser2P_img,...
		Ser2P_blurRange./thisPixelSize);
	
	thisImgRange = round(imgRanges{pp}./thisPixelSize)+1;
	OP_img = OP_img(...
		thisImgRange(1):thisImgRange(2),...
		thisImgRange(3):thisImgRange(4));
	Ser2P_img = Ser2P_img(...
		thisImgRange(1):thisImgRange(2),...
		thisImgRange(3):thisImgRange(4));
	Ser5P_img = Ser5P_img(...
		thisImgRange(1):thisImgRange(2),...
		thisImgRange(3):thisImgRange(4));
	thisSize = size(Ser5P_img);
	
	OP_lims = prctile(OP_img(:),[10,99.99]);
	Ser5P_img = Ser5P_img./median(Ser5P_img(:));
	Ser5P_lims = [0.0,4.5];
	Ser2P_lims = prctile(Ser2P_img(:),[10,99.9]);
	
	
	subplot(3,numPlots.*2,numPlots.*0+1+(pp-1).*2)
	imagesc([0,thisSize(2)].*thisPixelSize,...
		[0,thisSize(1)].*thisPixelSize,...
		OP_img,OP_lims)
	axis equal tight
	set(gca,'XTick',[],'YTick',[])
	set(gca,'YDir','normal')
	hold on
	plot([1,1+scaleBar],thisSize(1).*thisPixelSize-[1,1],...
		'w-','LineWidth',3) % Plot scale bar
	textObject = ...
		text(1+0.5.*scaleBar,thisSize(2).*thisPixelSize-2.0,...
		sprintf('%d \\mum',scaleBar));
	set(textObject,'Color',[1,1,1],'HorizontalAlignment','center',...
		'FontSize',12)
	ylabel('Oligopaint')
	
	
	subplot(3,numPlots.*2,numPlots.*2+1+(pp-1).*2)
	imagesc([0,thisSize(2)].*thisPixelSize,...
		[0,thisSize(1)].*thisPixelSize,...
		Ser5P_img,Ser5P_lims)
	axis equal tight
	set(gca,'XTick',[],'YTick',[])
	set(gca,'YDir','normal')
	ylabel('Pol II S5P')
	colormap(gray)

	
	subplot(3,numPlots.*2,numPlots.*4+1+(pp-1).*2)
	imagesc([0,thisSize(2)].*thisPixelSize,...
		[0,thisSize(1)].*thisPixelSize,...
		Ser2P_img,Ser2P_lims)
	axis equal tight
	set(gca,'XTick',[],'YTick',[])
	set(gca,'YDir','normal')
	ylabel('Pol II S2P')
	colormap((gray))
	
	
	
	% --- OP / Ser5P
	% Prepare the color channels for RGB plotting
	this_magenta_plot = ...
		(Ser5P_img-Ser5P_lims(1))./diff(Ser5P_lims);
	this_green_plot = ...
		(OP_img-OP_lims(1))./diff(OP_lims);
	thisSize = size(OP_img);
	
	subplot(3,numPlots.*2,numPlots.*0+2+(pp-1).*2)
	
	redChannel = this_magenta_plot;
	blueChannel = this_magenta_plot;
	greenChannel = this_green_plot;
	rgb_img = zeros(thisSize(1),thisSize(2),3);
	rgb_img(:,:,1) = redChannel;
	rgb_img(:,:,2) = greenChannel;
	rgb_img(:,:,3) = blueChannel;
	
	image([0,thisSize(2)].*thisPixelSize,...
		[0,thisSize(1)].*thisPixelSize,...
		rgb_img)
	set(gca,'XTick',[],'YTick',[])
	set(gca,'YDir','normal')
	ylabel('OP / Pol II Ser5P')
	axis equal tight
	
	% --- OP / Ser2P
	% Prepare the color channels for RGB plotting
	this_magenta_plot = ...
		(Ser2P_img-Ser2P_lims(1))./diff(Ser2P_lims);
	this_green_plot = ...
		(OP_img-OP_lims(1))./diff(OP_lims);
	thisSize = size(OP_img);
	
	subplot(3,numPlots.*2,numPlots.*2+2+(pp-1).*2)
	
	redChannel = this_magenta_plot;
	blueChannel = this_magenta_plot;
	greenChannel = this_green_plot;
	rgb_img = zeros(thisSize(1),thisSize(2),3);
	rgb_img(:,:,1) = redChannel;
	rgb_img(:,:,2) = greenChannel;
	rgb_img(:,:,3) = blueChannel;
	
	image([0,thisSize(2)].*thisPixelSize,...
		[0,thisSize(1)].*thisPixelSize,...
		rgb_img)
	set(gca,'XTick',[],'YTick',[])
	set(gca,'YDir','normal')
	ylabel('OP / Pol II Ser2P')
	axis equal tight
	
	% --- OP / Ser2P
	% Prepare the color channels for RGB plotting
	this_magenta_plot = ...
		(Ser5P_img-Ser5P_lims(1))./diff(Ser5P_lims);
	this_green_plot = ...
		(Ser2P_img-Ser2P_lims(1))./diff(Ser2P_lims);
	thisSize = size(OP_img);
	
	subplot(3,numPlots.*2,numPlots.*4+2+(pp-1).*2)
	
	redChannel = this_magenta_plot;
	blueChannel = this_magenta_plot;
	greenChannel = this_green_plot;
	rgb_img = zeros(thisSize(1),thisSize(2),3);
	rgb_img(:,:,1) = redChannel;
	rgb_img(:,:,2) = greenChannel;
	rgb_img(:,:,3) = blueChannel;
	
	image([0,thisSize(2)].*thisPixelSize,...
		[0,thisSize(1)].*thisPixelSize,...
		rgb_img)
	set(gca,'XTick',[],'YTick',[])
	set(gca,'YDir','normal')
	ylabel('Pol iI Ser5P / Pol II Ser2P')
	axis equal tight
	
end