%% This is structural quality analysis for quality assurance.
% transversal T1-weighted multi-slice spin-echo sequence:
% write_QA_Tran_T1.m
% clear all ; close all; clc ;
magFactor = 100 ;
lineWidth = 1.5 ;
fontSize = 24 ;

%% Step 1: load image data (two acquisitions)
imfile = dir('*.nii') ;
im_all = double(niftiread(imfile(1).name)) ;
im_all = reshape(im_all, [256, 256, 11, 2]) ;

figure ;
montage(mat2gray(im_all(:,:,:,1))) ;
title('11 slices of the first acquisition') ;
colormap default ;
figure ;
montage(mat2gray(im_all(:,:,:,1)-im_all(:,:,:,2))) ;
title('The difference between the first and second acquisitions') ;
colormap default ;
% the center slice of the first acquisition is used for structural quality analysis
im = squeeze(im_all(:,:,6,1)) ;
% calculate pixel size for background detection
fovx = 250 ; fovy = 250 ; % FOV in [mm]
nx = 256 ; ny = 256 ; % matrix size in x and y directions
pxsize = (fovx/nx) * (fovy/ny) ; % pixel size in mm*mm

%% generate phantom and large ROI masks
im_edge = double(edge(im) ) ; % detect circular edge
[im_X, im_Y] = ind2sub(size(im_edge), find (im_edge)) ;
% fit for circular center (Xc, Yc), and radius: R
[Xc, Yc, R] = circfit(im_X, im_Y) ;
pmsk = makeCircleMask([size(im,1) size(im,2)], Xc, Yc, R) ; % phantom mask
% large ROI mask, whose radius is 0.8 * R
largeROI_mask = makeCircleMask([size(pmsk,1) size(pmsk,2)], Xc, Yc, R*0.8) ; % large ROI mask
largeROI = largeROI_mask .* im ;
Sroi = sum(largeROI(:)) / nnz(largeROI_mask) ; % the mean signal intensity inside the large ROI
Xc = round(Xc) ;
Yc = round(Yc) ;
R = round(R) ;

figure ;
t = tiledlayout(1, 3) ;
nexttile
im_edge_fuse = labeloverlay(mat2gray(im), im_edge,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_edge_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Edge mask','FontSize', fontSize) ;
nexttile
im_phantom_mask_fuse = labeloverlay(mat2gray(im), pmsk,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_phantom_mask_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Phantom mask','FontSize', fontSize) ;

nexttile
im_largeROI_mask_fuse = labeloverlay(mat2gray(im), largeROI_mask,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_largeROI_mask_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Large ROI mask','FontSize', fontSize) ;
t.Padding = 'none';
t.TileSpacing = 'none';

%% generate four background masks
% bkg_left and bkg_right: background noise in readout direction
% bkg_btm and bkg_top: background noise in phase-encoded direction
ROI_area = 10 * 100 ; % each backgorund ROI area = 10 cm2
% the background ROI lies 0.1*R from the edge of the phantom as well as
% from the border of the image.
R1 = ceil(R * 0.1) ;
% left background ROI
Wl = size(R1:Xc-R-R1+1, 2) ;
Hl2 = round(ROI_area / Wl / pxsize / 2) ; % half of the height
bkg_left = zeros(size(im)) ;
bkg_left(R1:Xc-R-R1+1,Yc-Hl2+1:Yc+Hl2) = 1 ;
ROIbkg_left = bkg_left .* im ;
Sbkg_left = sum(ROIbkg_left(:)) / nnz(bkg_left) ;
ROIbkg_left_std = std(ROIbkg_left(:) ) ;

% right background ROI
Wr = size(Xc+R+R1-1:nx-R1+1, 2) ;
Hr2 = round(ROI_area / Wr / pxsize / 2) ;
bkg_right = zeros(size(im)) ;
bkg_right(Xc+R+R1-1:nx-R1+1,Yc-Hr2+1:Yc+Hr2) = 1 ;
ROIbkg_right = bkg_right .* im ;
Sbkg_right = sum(ROIbkg_right(:)) / nnz(bkg_right) ;
ROIbkg_right_std = std(ROIbkg_right(:) ) ;

% top background ROI
Wt = size(Yc+R+R1-1:ny-R1+1, 2) ;
Ht2 = round(ROI_area / Wt / pxsize / 2) ;
bkg_top = zeros(size(im)) ;
bkg_top(Xc-Ht2+1:Xc+Ht2,Yc+R+R1-1:ny-R1+1) = 1 ;
ROIbkg_btm = bkg_top .* im ;
Sbkg_btm = sum(ROIbkg_btm(:)) / nnz(bkg_top) ;

% Bottom background ROI
Wb = size(R1:Yc-R-R1+1, 2) ;
Hb2 = round(ROI_area / Wb / pxsize / 2) ;
bkg_btm = zeros(size(im)) ;
bkg_btm(Xc-Hb2+1:Xc+Hb2,R1:Yc-R-R1+1) = 1 ;
ROIbkg_top = bkg_btm .* im ;
Sbkg_top = sum(ROIbkg_top(:)) / nnz(bkg_btm) ;

figure ;
t = tiledlayout(1, 4) ;
nexttile
im_bkg_left_fuse = labeloverlay(mat2gray(im), bkg_left,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_bkg_left_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Left bkg ROI','FontSize', fontSize) ;
nexttile
im_bkg_right_fuse = labeloverlay(mat2gray(im), bkg_right,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_bkg_right_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Right bkg ROI','FontSize', fontSize) ;

nexttile
im_bkg_btm_fuse = labeloverlay(mat2gray(im), bkg_btm,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_bkg_btm_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Bottm bkg ROI','FontSize', fontSize) ;
nexttile
im_bkg_top_fuse = labeloverlay(mat2gray(im), bkg_top,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_bkg_top_fuse),[],'InitialMagnification', magFactor) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Top bkg ROI','FontSize', fontSize) ;

t.Padding = 'none';
t.TileSpacing = 'none';
% exportgraphics(gcf,'se.png',...
%     'ContentType','vector',...
%     'BackgroundColor','none') ;

%% Find the small ROI (area = 100 mm*mm) with maximum/minimum mean intensity within the large ROI
piu_area = 100 ; % the kernel area of PIU calcualtion, [mm*mm]
R2 = ceil(sqrt(piu_area / pxsize / pi)) ; % the radius of the circular kernel
convKernel = makeCircleMask([2*ceil(R2) 2*ceil(R2)], R2, R2, R2) ; % make a kernel with a circular mask for convolution
% convolve large ROI with the kernel, divided by the number of pixels
% inside the kernel to get the mean intensity within the running kernel
smallROI_mean = convn(largeROI, convKernel, 'same') ./ convn(ones(size(largeROI)), convKernel, 'same') ;
% Find a region with lowest mean intensity inside the large ROI
voxelNum = convn(largeROI_mask, convKernel, 'same') ;
isInLargeROI = voxelNum ;
isInLargeROI(isInLargeROI<sum(convKernel(:))) = 0 ; % is the moving box totally located within the large ROI?
isInLargeROI(isInLargeROI==sum(convKernel(:))) = 1 ;
smallROI_mean = smallROI_mean .* isInLargeROI ;
smallROI_mean(smallROI_mean == 0) = NaN ;
% find a region with the highest mean intensity
smallROI_meanMax = max(smallROI_mean(:)) ;
[X_ROImax, Y_ROImax] = ind2sub(size(smallROI_mean), find (smallROI_mean == smallROI_meanMax));
% find a region with the lowest mean intensity
smallROI_meanMin = min(smallROI_mean(:)) ;
[X_ROImin, Y_ROImin] = ind2sub(size(smallROI_mean), find (smallROI_mean == smallROI_meanMin));

largeROI_mask_edge = double(edge(largeROI_mask)) ;
smallROI_min_mask = zeros([size(im,1) size(im,2)]) ;
smallROI_max_mask = zeros([size(im,1) size(im,2)]) ;
smallROI_min_mask(X_ROImin:X_ROImin+2*R2-1, Y_ROImin:Y_ROImin+2*R2-1) = convKernel ;
smallROI_max_mask(X_ROImax:X_ROImax+2*R2-1, Y_ROImax:Y_ROImax+2*R2-1) = convKernel ;
figure ;
t = tiledlayout(1, 2) ;
nexttile
mask_temp = largeROI_mask_edge+smallROI_max_mask ;
mask_temp(mask_temp>0) = 1 ;
im_smallROI_max_fuse = labeloverlay(mat2gray(im), mask_temp,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_smallROI_max_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('small ROI with maximum mean intensity','FontSize', fontSize) ;
nexttile
mask_temp = largeROI_mask_edge + smallROI_min_mask ;
mask_temp(mask_temp>0) = 1 ;
im_smallROI_min_fuse = labeloverlay(mat2gray(im), mask_temp,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_smallROI_min_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('small ROI with minimum mean intensity','FontSize', fontSize) ;
t.Padding = 'none';
t.TileSpacing = 'none';

%% Calculate percent signal ghosting (PSG)
PSG = 100 * abs( (Sbkg_left+Sbkg_right) - (Sbkg_btm + Sbkg_top) ) / (2*Sroi) ; % in percent
disp(['Percent signal ghosting (PSG) = ', num2str(PSG), '%']) ;

%% Calculate percentage Image Uniformity (PIU)
PIU = 100 * (1 - (smallROI_meanMax - smallROI_meanMin) / (smallROI_meanMax + smallROI_meanMax)) ;
disp(['percent image uniformity (PIU) = ', num2str(PIU), '%']) ;

%% Calculate signal-to-noise ratio
% NEMA SNR, 1 slice
SNR_NEMA1 = 2 * 0.655 * Sroi / (ROIbkg_left_std + ROIbkg_right_std) ;
disp(['NEMA SNR1 = ', num2str(SNR_NEMA1), ' using the first acquisition']) ;

% NEMA SNR, 2 slices
im_mean = (im + squeeze(im_all(:,:,6,2))) / 2 ;
im_combine_ROI_mean = largeROI_mask .* im_mean ;
ROI_mean_mean = sum(im_combine_ROI_mean(:)) / nnz(largeROI_mask) ;

im_diff = abs(im(:,:)) - abs(squeeze(im_all(:,:,6,2))) ;
im_ROI_diff = largeROI_mask .* im_diff ;
im_ROI_diff(im_ROI_diff==0) = nan ;
ROI_diff_std = std(im_ROI_diff(:), 'omitnan') ;
SNR_NEMA2 = 1.41 * ROI_mean_mean / ROI_diff_std ;
disp(['NEMA SNR2 = ', num2str(SNR_NEMA2), ' using 2 acquisitions']) ;
figure ;
subplot 221 ;
imshow(im_mean, []) ;
title('mean of 2 acq.') ;
subplot 222 ;
imshow(im_combine_ROI_mean, []) ;
title('ROI mean of 2 acq.') ;
subplot 223 ;
imshow(im_diff, []) ;
title('im diff of 2 acq.') ;
subplot 224 ;
imshow(im_ROI_diff, []) ;
title('ROI diff of 2 acq.') ;
NEMA_SNR_ratio = SNR_NEMA1 / SNR_NEMA2 ;
disp(['NEMA SNR ratio = ', num2str(NEMA_SNR_ratio),]) ;
% spoiler in PE direction at the end of each acquisition, or add a delay between two acqusitions???

%%
allmask = 2*bkg_left+2*bkg_right+2*bkg_top+2*bkg_btm + edge(largeROI_mask) + im_edge ;
figure ;
im_allMask_fuse = labeloverlay(mat2gray(im), allmask,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_allMask_fuse),[],'InitialMagnification', magFactor) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('allMask','FontSize', fontSize) ;

%%
smallROImask = 1*edge(largeROI_mask) + 1*smallROI_min_mask + 2*smallROI_max_mask ;
figure ;
im_allMask_fuse = labeloverlay(mat2gray(im), smallROImask,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_allMask_fuse),[],'InitialMagnification', magFactor) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('allMask','FontSize', fontSize) ;