function MRI_Abdominal_Segmentation(file, dirsave, saveName, slices, spineExtract, debug)
% MRI fat segmentation
% Performs bias-correction and segmentation of abdominal MRI images. 
% The compartments extracted are the visceral and the sub-cutaneous fat, 
% with an option for splitting the sub-cutaneous into a deep and superficial
% part.
% ------------------------------
% INPUT:
% file: 
%       Full path to data file. 
%       Data file MUST be named 'data' and MUST be formatted as struct:
%       -   data.data: Scan volume in double
%           data.info: DICOM header
%           data.SliceLocation: [x, y, z] in scanner coordiante system
% dirsave:
%       Full path to where results are saved
% saveName:
%       Name results are saved under
% slices:
%       Slices of the volume to segment. Leave empty [], if all slices are
%       to be used
% SpineExtract:
%       0: no spine extraction
%       1: spine extraction is performed
% debug:
%       0: no debug images outputtet
%       1: debug images outputtet
%       WARNING: LOTS OF IMAGES OPENEND!
% ------------------------------
% OUTPUT:
%       Images: Images of all slices with the segmentations shown
%           fig_intra: Visceral/intra-abdominal fat
%           fig_layer: All regions shown
%           fig_sub: Subcutaneous fat
%           fig: Images with no segmentation
%       Abdomen.mat: Struct with following fields
%           data: the bias-corrected image data
%           sub: Subcutaneous mask
%           subS: Superficial suncutaneous mask
%           subD: Deep subcutaneous mask
%           intra: Intra abdominal/visceral fat
%           total: Mask of entire slice
%           VoxelSize: The voxel dimensions
%       Abdomen.xls: Excel-file with areas for each slice
%
% ------------------------------
% Please cite the following work when using:
%
% Christensen A.N. et al. (2017) 
% Automatic Segmentation of Abdominal Fat in MRI-Scans, 
% Using Graph-Cuts and Image Derived Energies. 
% In: Sharma P., Bianchi F. (eds) Image Analysis. SCIA 2017. 
% Lecture Notes in Computer Science, vol 10270. Springer, Cham
% DOI: 10.1007/978-3-319-59129-2_10
%
%
% Copyright Anders Nymark Christensen, anym@dtu.dk, DTU Compute 2008
%
%
% The software in 'GraphCut' and 'BiasCorr' are under SEPARATE licensces! 
%
%
% The "Software" is limited to this function and the software contained 
% in the following folders:
% 'grid_cut'
% 'cost_functions'
% 'auxiliary_functions'
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.




%% Check inputs
if ~isstruct(load(file))
    error('First input must be the full path to a data file STRUCT')
end
    
if nargin < 5
    spineExtract = 0;
end
if nargin < 6
    debug = 0;
end

%% Add paths
addpath 'grid_cut' 
addpath 'GraphCut' 
addpath 'cost_functions'
addpath 'auxiliary_functions' 
addpath(genpath('BiasCorr'))

%% Find System path style
if ispc
    slash = '\';
elseif isunix
    slash = '/';
else
    error('System unknown. Try hardcoding slash and septoken')
end

%% Load data and create folder for results
load(file)

if ~exist([dirsave,slash,saveName],'dir')
    mkdir([dirsave,slash,saveName])
end

%% stiffness in radial and height direction
drdz = [3 10; 4 10]; % [inSliceInner betweenSliceInner; inSliceOuter betweenSliceOuter]
dldu = [2 150]; % Distance between inner and outer [min max]

% hard coding resolution for polar coordinates
angle_resolution = 540;
radius_resolution = 300;

%% Bias Correction I - masking

% Apply slice slection
abdomen = double(data.data);
if nargin < 3
    slices = 1:size(abdomen,3);
elseif isempty(slices)
    slices = 1:size(abdomen,3);
end
abdomen = abdomen(:,:,slices);

% Make rough mask of data
[~, C] = kmeans(abdomen(:),5);
BW = abdomen > min(C)*3;

% Additional opening and erosion if object touches image border
if isempty(bwareaopen(imclearborder(BW),2500))
    BW = imopen(BW,strel('disk',7));
    BW = imerode(BW,strel('disk',7));
end

% For each slice, select largest area as the mask
% Dilate and fill holes - find centroid of mask
for i=1:size(BW,3)
    [L, num] = bwlabel(BW(:,:,i));
    if num > 1
        areas = zeros(num,1);
        for labels = 1:num
            areas(labels) = bwarea(L == labels);
        end
        [~, label] = max(areas);
        BW(:,:,i) = (L == label);
    end
    BW(:,:,i) = imdilate(BW(:,:,i),strel('disk',4));
    BW(:,:,i) = imfill(BW(:,:,i),'holes');
    C = regionprops(BW,'centroid');
end
% Adjust centroid in the posterior direction
cent = zeros(length(C),2);
for i = 1:length(C)
    cent(i,:) = C(i).Centroid(1:2);
end
cent(2) = 0.95*cent(2);
cent = fliplr(mean(cent,1));
abdomenTemp = abdomen.*BW;


%% Bias-Correction II - correction

% Bias-correction using the rough mask
[correctedData, ~, ~] = iic(abdomenTemp, ...
            'stepSize', [data.info.PixelSpacing(1), data.info.PixelSpacing(2), data.info.SliceThickness], ...
            'desiredStepSize', 4, ...
            'numberOfGaussians', 8, ...
            'smoothingDistance', 25, ...
            'smoothingRegularization', 5e-1, ...            
            'showPlots', false, ...
            'verbosity', 'none',...
            'mask', double(abdomenTemp > 0) );% 'backgroundThresholdType', 'otzu'

% Refine mask using Chan-Vese
for index = 1:size(BW,3)
    BW(:,:,index) = activecontour(correctedData(:,:,index),BW(:,:,index),100,'chan-vese','SmoothFactor',0.5,'ContractionBias',0.2);
end

% Output mask as images
if debug
    for index = 1:size(BW,3)
        if mod(index,9) == 1
            figure
        end
        subplot(3,3,mod(index,9))
        imagesc(BW)
        title(['Slice: ',num2str(index)])
    end
end

% Bias-Correction using refined mask
abdomen = abdomen.*BW;
[correctedData, ~, ~] = iic(abdomen, ...
            'stepSize', [data.info.PixelSpacing(1), data.info.PixelSpacing(2), data.info.SliceThickness], ...
            'desiredStepSize', 4, ...
            'numberOfGaussians', 8, ...
            'smoothingDistance', 25, ...
            'smoothingRegularization', 5e-1, ...            
            'showPlots', false, ...
            'verbosity', 'none',...
            'mask', double(abdomen > 0) );% 'backgroundThresholdType', 'otzu'

% Change orientation of data to match graph-cut routines        
V = permute(correctedData,[2 1 3]);  
dimV = size(V);

% volume coordinates: x,y,height; terrain coordinates: angle,height,radius
unfolded = unfold_volume(V,angle_resolution,radius_resolution,cent);

% Smoothing data
medMask = [round(10/data.info.PixelSpacing(1)) round(10/data.info.PixelSpacing(2))];

unfoldedM = zeros(size(unfolded));
unfoldedS = zeros(size(unfolded));
unfoldedG = zeros(size(unfolded));

for i = 1:size(unfolded,2)
    unfoldedM(:,i,:) = medfilt2(squeeze(unfolded(:,i,:)),medMask,'symmetric');
    unfoldedS(:,i,:) = imfilter(medfilt2(squeeze(unfolded(:,i,:)),medMask), fspecial('Gaussian', round(1.5*medMask), 3/data.info.PixelSpacing(1)),'symmetric');
    unfoldedG(:,i,:) = imfilter(squeeze(unfolded(:,i,:)), fspecial('Gaussian', 5, 2),'symmetric');
    
    if debug
        figure
        title(['Slice: ',num2str(i)])
        subplot(2,2,1)
        imagesc(squeeze(unfolded(:,i,:)))
        subplot(2,2,2)
        imagesc(squeeze(unfoldedM(:,i,:)))
        subplot(2,2,3)
        imagesc(squeeze(unfoldedS(:,i,:)))
        subplot(2,2,4)
        imagesc(squeeze(unfoldedG(:,i,:)))
    end
end


%% Find edge measures

% outer surface edge energies
eo = detectEdge(unfolded, 1);
eo_S = detectEdge(unfoldedS, 1);
eo_M = detectEdge(unfoldedM, 1);
eo = (eo-min(eo(:)))/(max(eo(:))-min(eo(:)));
eo_S = (eo_S-min(eo_S(:)))/(max(eo_S(:))-min(eo_S(:)));
eo_M = (eo_M-min(eo_M(:)))/(max(eo_M(:))-min(eo_M(:)));


% Inner surface edge energies
ei = detectEdge(unfolded, -1 );
ei_G = detectEdge(unfoldedG, -1 );
ei_M = detectEdge(unfoldedM, -1 );

ei = (ei-min(ei(:)))/(max(ei(:))-min(ei(:)));
ei_G = (ei_G-min(ei_G(:)))/(max(ei_G(:))-min(ei_G(:)));
for i = 1:size(ei,2)
   ei_M(:,i,:) = medfilt2(squeeze(ei(:,i,:)),[round(70/data.info.PixelSpacing(2)) round(4/data.info.PixelSpacing(1))],'symmetric'); %ei_M(:,i,:) = medfilt2(squeeze(ei(:,i,:)),[55 3],'symmetric');
end
ei_M = (ei_M-min(ei_M(:)))/(max(ei_M(:))-min(ei_M(:)));

%% Find region measures

% Outer surface region energy
outerSurf = detectEdge(unfolded, -1);
outerSurf(outerSurf < 0) = 0;
outerSurf = on_surface_cumulative_cost(outerSurf,1,0);

% Inner surface region energy
edgeIm = detectEdge(unfolded, 1);
edgeIm(edgeIm < 0) = 0;

surfCost = on_surface_cumulative_cost(edgeIm,1,0);
for i = 1:size(surfCost,2)
   surfCost(:,i,:) = medfilt2(squeeze(surfCost(:,i,:)),medMask,'symmetric');
   surfCost(:,i,:) = (surfCost(:,i,:)-min(min(surfCost(:,i,:))))/(max(max(surfCost(:,i,:)))-min(min(surfCost(:,i,:))));
end
surfCost2 = on_surface_cumulative_cost(unfolded,1,0);

surfCost = (surfCost-min(surfCost(:)))/(max(surfCost(:))-min(surfCost(:)));
surfCost2 = (surfCost2-min(surfCost2(:)))/(max(surfCost2(:))-min(surfCost2(:)));
surfCost3 = surfCost.*surfCost2;
for i = 1:size(surfCost3,2)
   surfCost3(:,i,:) = medfilt2(squeeze(surfCost3(:,i,:)),[round(70/data.info.PixelSpacing(2)) round(4/data.info.PixelSpacing(1))],'symmetric');
end
surfCost3 = (surfCost3-min(surfCost3(:)))/(max(surfCost3(:))-min(surfCost3(:)));

%% radial energy dependence
latFactor = 1.1;
LateralIncrease = ones(size(surfCost,1),1);
LateralIncrease(5*round(size(surfCost,1)/6):end) = latFactor; % 5/6 to end
LateralIncrease(1:round(size(surfCost,1)/10)) = latFactor; % begin to 1/10
LateralIncrease(round(4*size(surfCost,1)/10):round(2*size(surfCost,1)/3)) = latFactor; % 4/10 to 2/3
LateralIncrease = repmat(LateralIncrease,1,size(surfCost,3));

for i=1:size(surfCost3,2)
    surfCost3(:,i,:) = squeeze(surfCost3(:,i,:)) .* LateralIncrease;
    surfCost(:,i,:) = squeeze(surfCost(:,i,:)) .* LateralIncrease;
    if i > round(size(surfCost3,2)/2) + round(size(surfCost3,2)/8)
        surfCost3(:,i,:) = surfCost3(:,i,:)*3;
        surfCost(:,i,:) = surfCost(:,i,:)*3;
    elseif i > round(size(surfCost3,2)/2)
        surfCost3(:,i,:) = surfCost3(:,i,:)*1.5;
        surfCost(:,i,:) = surfCost(:,i,:)*1.5;
    end
end

%% Final energy
ei2 = 0.5*ei + 0.1*ei_G + 0.1*ei_M + 0.05*surfCost + 0*surfCost2 + 0.2*surfCost3;
eo2 = 0.0*eo + 0.25*eo_S + 0.0*eo_M + 0.75*outerSurf;

% Find cut
testCost = cat(4,ei2,eo2);
so2 = grid_cut(testCost,[],drdz,[1 0],dldu);

% Smooth inner cut
for i = 1:size(so2,2)
   so2(:,i,1) = medfilt1(so2(:,i,1),5);
    dummy = medfilt1(so2(:,i,1),35); 
    so2(round(size(so2,1)*5/8):round(size(so2,1)*7/8),i,1) = dummy(round(size(so2,1)*5/8):round(size(so2,1)*7/8)); 
end

% Fold back to image coordiantes
So2 = fold_back(so2,dimV,radius_resolution,angle_resolution,cent);


if debug
    figure, surf_tubular({So2})
    figure, show_tubular(V,{So2},'yes')
    
    figure, show_terrain(unfolded,{so2},'yes')
    figure, show_terrain(ei2,{so2},'yes')
    figure, show_terrain(ei,{so2},'yes')
    figure, show_terrain(ei_G,{so2},'yes')
    figure, show_terrain(ei_M,{so2},'yes')
    figure, show_terrain(surfCost,{so2},'yes')
    figure, show_terrain(surfCost2,{so2},'yes')
    figure, show_terrain(surfCost3,{so2},'yes')
end

%% SCARPAS FASCIA
h = fspecial('log',[9 9],1);
em = imfilter(unfolded,h);
em = log(1./(em-min(em(:))));
em = em - min(em(:));

em_cut = ones(size(em));
em_cut_raw = max(unfolded(:))* ones(size(em));

surfCostScarpa = zeros(size(surfCost));
outer_cut = 2;
newso2 = so2;

for point = 1:size(so2,1)
    % loop over slices
    for slice = 1:size(so2,2)
        thickness = so2(:,slice,2)-so2(:,slice,1);
        if point > 7*size(so2,1)/32 && point < 9*size(so2,1)/32 % posterior part
            inner_cut = 0;
            surfScale = 0.1;
        elseif point > size(so2,1)/2 % Anterior part
            inner_cut = ceil(min(thickness)/5);
            surfScale = 2;
        else
            inner_cut = ceil(min(thickness)/10); % the rest
            surfScale = 1;
        end
        
        em_cut(point,slice,newso2(point,slice,1)+inner_cut:newso2(point,slice,2)-outer_cut) = em(point,slice,newso2(point,slice,1)+inner_cut:newso2(point,slice,2)-outer_cut);
        em_cut_raw(point,slice,newso2(point,slice,1)+inner_cut:newso2(point,slice,2)-outer_cut) = unfolded(point,slice,newso2(point,slice,1)+inner_cut:newso2(point,slice,2)-outer_cut);
        
        surfCostScarpa(point,slice,:) = surfScale*surfCost(point,slice,:);
    end
end

em_cut = (em_cut-min(em_cut(:)))/(max(em_cut(:))-min(em_cut(:)));
em_cut_raw = (em_cut_raw-min(em_cut_raw(:)))/(max(em_cut_raw(:))-min(em_cut_raw(:)));

testCost = 0.1 * surfCostScarpa + 0.8 * em_cut + 0.1 *em_cut_raw;
io2 = grid_cut(testCost,[],[2 3],[1 0],[1 3]);
dummy = io2;
for i = 1:size(io2,2)
   io2(:,i) = medfilt1(io2(:,i),15); 
   io2(1,i) = dummy(1,i);
end

Io2 = fold_back(io2,dimV,radius_resolution,angle_resolution,cent);

if debug
    figure, show_tubular(V,{Io2},'yes')
    figure, show_tubular(V,{So2,Io2},'yes')
    
    figure, show_terrain(testCost,{io2},'yes')
    figure, show_terrain(surfCostScarpa,{io2},'yes')
    figure, show_terrain(em_cut,{io2},'yes')
    figure, show_terrain(em_cut_raw,{io2,so2},'yes')
end

%% Spine extraction
if spineExtract
    
    % Adjust center more posterior
    newCent = cent;
    newCent(1) = newCent(1)*1.2;
    unfolded2 = unfold_volume(V,angle_resolution,radius_resolution,newCent);
    unfolded2M = zeros(size(unfolded2));
    for i = 1:size(unfolded2,2)
        unfolded2M(:,i,:) = medfilt2(squeeze(unfolded2(:,i,:)),medMask,'symmetric');
    end
    
    % Define energies
    spineCut = zeros(size(unfolded2));
    spineCut(1:50,:,:) = unfolded2(1:50,:,:);
    spineCut(230:end,:,:) = unfolded2(230:end,:,:);
    
    spineCutM = zeros(size(unfolded2M));
    spineCutM(1:50,:,:) = unfolded2M(1:50,:,:);
    spineCutM(230:end,:,:) = unfolded2M(230:end,:,:);
    
    surfCostSpine = on_surface_cumulative_cost(spineCut,0,0);
    surfCostSpine2 = on_surface_cumulative_cost(spineCut,1,0);

    test = detectEdge(spineCut, -1);
    testM = detectEdge(spineCutM, -1);

    spineCut = (spineCut-min(spineCut(:)))/(max(spineCut(:))-min(spineCut(:)));
    test = (test-min(test(:)))/(max(test(:))-min(test(:)));
    testM = (testM-min(testM(:)))/(max(testM(:))-min(testM(:)));
    surfCostSpine = (surfCostSpine-min(surfCostSpine(:)))/(max(surfCostSpine(:))-min(surfCostSpine(:)));
    surfCostSpine2 = (surfCostSpine2-min(surfCostSpine2(:)))/(max(surfCostSpine2(:))-min(surfCostSpine2(:)));
    
    % Set hard geometrical/anatomical constraints
    surfCostSpine(:,:,1:round(size(surfCostSpine,3)/30)) = inf; 
    surfCostSpine2(:,:,round(size(surfCostSpine2,3)/4):end) = inf;  
    surfCostSpine2(round(2*size(surfCostSpine2,1)/3):round(4*size(surfCostSpine2,1)/5),:,1:round(size(surfCostSpine2,3)/9)) = inf; 
    surfCostSpine2(round(4*size(surfCostSpine2,1)/7):round(14*size(surfCostSpine2,1)/15),:,1:round(size(surfCostSpine2,3)/15)) = inf; 
    
    
    % Energy, cut and fold back
    cutCost = 0.15*spineCut  + 0.15*test + 0.25*testM + 0.25*surfCostSpine + 0.20*surfCostSpine2;%cat(4,test+2*surfCostSpine2,test+surfCostSpine);
    sp2 = grid_cut(cutCost,[],[3 3],[1 0],[5 300]);
    Sp2 = fold_back(sp2,dimV,radius_resolution,angle_resolution,newCent);
    
    % Show images if debug
    if debug
        figure, show_tubular(V,{Sp2},'yes')
        figure, show_terrain(spineCut,{sp2},'yes')
        
        figure, show_terrain(cutCost,{sp2},'yes')
        
        figure, show_terrain(spineCut,{sp2},'yes')
        figure, show_terrain(test,{sp2},'yes')
        figure, show_terrain(testM,{sp2},'yes')
        figure, show_terrain(surfCostSpine,{sp2},'yes')
        figure, show_terrain(surfCostSpine2,{sp2},'yes')
    end
    
end

%% VISCERAL FAT SEGMENTATION

voxelSize = data.info.PixelSpacing(1)*data.info.PixelSpacing(2)*data.info.SliceThickness;

if exist('Sp2','var')
    spine = zeros(size(V));
end
total = zeros(size(V));
sub = zeros(size(V));
subS = zeros(size(V));
subD = zeros(size(V));
intra = zeros(size(V));

for index = 1:size(So2,2)
    % mask of subcutaneous fat
    BWinner = poly2mask(So2(:,index,1,1), So2(:,index,2,1)+1, size(V,1), size(V,2));
    BWouter = poly2mask(So2(:,index,1,2), So2(:,index,2,2), size(V,1), size(V,2));
    BWscarpa = poly2mask(Io2(:,index,1), Io2(:,index,2), size(V,1), size(V,2));
    
    % Check for spine compartment
    if exist('Sp2','var')
        BWspine = poly2mask(Sp2(:,index,1), Sp2(:,index,2), size(V,1), size(V,2));
        BWvat = BWinner & ~BWspine;
        spine(:,:,index) = BWspine;
    else
        BWvat = BWinner;
    end
    
    % Find areas
    dummy = V(:,:,index);
    total(:,:,index) = BWouter;
    sub(:,:,index) = BWouter-BWinner;
    subS(:,:,index) = BWouter-BWscarpa;
    subD(:,:,index) = BWscarpa-BWinner;
    
    % Find fat using k-means
    grayLevel = mean(dummy(BWvat));
    [~, C] = kmeans(dummy(BWvat),5,'MaxIter',1000,'start',[grayLevel-100 grayLevel-50 grayLevel grayLevel+50 grayLevel+100]');
    C = sort(C);

    intra(:,:,index) = dummy.*BWvat > C(3); % 4
end

%% OUTPUT
Abdomen.data.data = permute(V,[2 1 3]); 
Abdomen.data.total = permute(total,[2 1 3]); 
Abdomen.data.sub = permute(sub,[2 1 3]);
Abdomen.data.subS = permute(subS,[2 1 3]); 
Abdomen.data.subD = permute(subD,[2 1 3]); 
Abdomen.data.intra = permute(intra,[2 1 3]); 
Abdomen.data.voxelSize = voxelSize;
if spineExtract == 1
    Abdomen.data.spine = permute(spine,[2 1 3]);
end
  
V = permute(V,[2 1 3]);   
figure
for layer = 1:size(So2,2)
    B_total = bwboundaries(Abdomen.data.total(:,:,layer));
    B_scarpa = bwboundaries(Abdomen.data.subD(:,:,layer));
    B_sub = bwboundaries(Abdomen.data.sub(:,:,layer));
    B_intra = bwboundaries(Abdomen.data.intra(:,:,layer));
    
    
    
    clf;
    imagesc(V(:,:,layer))
    colormap gray
    axis image
    axis off
    hold on
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    axis off
    saveas(gcf,[dirsave,slash,saveName,slash,saveName,'fig',num2str(layer),'.png'])
    
    for i=1:length(B_sub)
        boundary2 = B_sub{i};
        plot(boundary2(:,2), boundary2(:,1), 'r', 'LineWidth', 1)
    end
    axis off
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    saveas(gcf,[dirsave,slash,saveName,slash,saveName,'fig_sub',num2str(layer),'.png'])
    
    if spineExtract == 1
        B_spine = bwboundaries(Abdomen.data.spine(:,:,layer));
        clf;
        imagesc(V(:,:,layer))
        colormap gray
        axis image
        axis off
        hold on
        for i=1:length(B_spine)
            boundary2 = B_spine{i};
            plot(boundary2(:,2), boundary2(:,1), 'r', 'LineWidth', 1)
        end
        axis off
        set(gca, 'LooseInset', get(gca, 'TightInset'));
        saveas(gcf,[dirsave,slash,saveName,slash,saveName,'fig_spine',num2str(layer),'.png'])
    end
    
    clf;
    imagesc(V(:,:,layer))
    colormap gray
    axis image
    axis off
    hold on
    for i=1:length(B_intra)
        boundary3 = B_intra{i};
        plot(boundary3(:,2), boundary3(:,1), 'y', 'LineWidth', 1)
    end
    axis off
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    saveas(gcf,[dirsave,slash,saveName,slash,saveName,'fig_intra',num2str(layer),'.png'])
    
    for i=1:length(B_scarpa)
        boundary4 = B_scarpa{i};
        plot(boundary4(:,2), boundary4(:,1), 'g', 'LineWidth', 1)
    end
    for i=1:length(B_intra)
        boundary3 = B_intra{i};
        plot(boundary3(:,2), boundary3(:,1), 'y', 'LineWidth', 1)
    end
    for i=1:length(B_sub)
        boundary2 = B_sub{i};
        plot(boundary2(:,2), boundary2(:,1), 'r', 'LineWidth', 1)
    end
    for i=1:length(B_total)
        boundary1 = B_total{i};
        plot(boundary1(:,2), boundary1(:,1), 'b', 'LineWidth', 1)
    end

    axis off
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    saveas(gcf,[dirsave,slash,saveName,slash,saveName,'fig_layer',num2str(layer),'.png'])
end
close all

% save mat file
save([dirsave,slash,saveName,slash,saveName,'Abdomen.mat'],'Abdomen')

% generate xls-file
volumeData = {'name' 'slice' 'Total (Voxel)' 'Sub (Voxels)' 'Intra (Voxels)'...
    'Total (mm3)' 'Sub (mm3)' 'Intra (mm3)'};
for index = 1:size(So2,2)
   volumeData{index+1,1} = saveName;
   volumeData{index+1,2} = index;
   volumeData{index+1,3} = nnz(Abdomen.data.total(:,:,index));
   volumeData{index+1,4} = nnz(Abdomen.data.sub(:,:,index));
   volumeData{index+1,5} = nnz(Abdomen.data.intra(:,:,index));
   volumeData{index+1,6} = nnz(Abdomen.data.total(:,:,index))*Abdomen.data.voxelSize;
   volumeData{index+1,7} = nnz(Abdomen.data.sub(:,:,index))*Abdomen.data.voxelSize;
   volumeData{index+1,8} = nnz(Abdomen.data.intra(:,:,index))*Abdomen.data.voxelSize;
end

xlswrite([dirsave,slash,saveName,slash,saveName,'Abdomen.xls'],volumeData);