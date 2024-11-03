%% This is structural quality analysis for quality assurance.
% transversal T1-weighted multi-slice spin-echo sequence:
% write_QA_Tran_T1.m
function [PSG, PIU, SNR_1acq, SNR_2acq, mask_all] = structuralQuality(im_all)
%% Step 1: load image data (two acquisitions)
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
mask_all.im_edge = im_edge ;
mask_all.phantom_mask = pmsk ;
mask_all.largeROI_mask = largeROI_mask ;

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
ROIbkg_top = bkg_top .* im ;
Sbkg_top = sum(ROIbkg_top(:)) / nnz(bkg_top) ;

% Bottom background ROI
Wb = size(R1:Yc-R-R1+1, 2) ;
Hb2 = round(ROI_area / Wb / pxsize / 2) ;
bkg_btm = zeros(size(im)) ;
bkg_btm(Xc-Hb2+1:Xc+Hb2,R1:Yc-R-R1+1) = 1 ;
ROIbkg_btm = bkg_btm .* im ;
Sbkg_btm = sum(ROIbkg_btm(:)) / nnz(bkg_btm) ;

mask_all.bkg_left = bkg_left ;
mask_all.bkg_right = bkg_right ;
mask_all.bkg_top = bkg_top ;
mask_all.bkg_btm = bkg_btm ;

%% Find the small ROI (area = 100 mm*mm) with maximum/minimum mean intensity within the large ROI
piu_area = 100 ; % the kernel area of PIU calcualtion, [mm*mm]
R2 = ceil(sqrt(piu_area / pxsize / pi)) ; % the radius of the circular kernel
convKernel = makeCircleMask([2*ceil(R2) 2*ceil(R2)], R2, R2, R2) ; % make a kernel with a circular mask for convolution
% convolve large ROI with the kernel, divided by the number of pixels
% inside the kernel to get the mean intensity within the running kernel
smallROI_mean = convn(largeROI, convKernel, 'valid') ./ convn(ones(size(largeROI)), convKernel, 'valid') ;
% Find a region with lowest mean intensity inside the large ROI
voxelNum = convn(largeROI_mask, convKernel, 'valid') ;
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
smallROI_min_mask(X_ROImin:X_ROImin+2*R2-1, Y_ROImin:Y_ROImin+2*R2-1) = flip(flip(convKernel,1),2) ;
smallROI_max_mask(X_ROImax:X_ROImax+2*R2-1, Y_ROImax:Y_ROImax+2*R2-1) = flip(flip(convKernel,1),2) ;

mask_all.smallROI_min_mask = smallROI_min_mask ;
mask_all.smallROI_max_mask = smallROI_max_mask ;
mask_all.largeROI_mask_edge = largeROI_mask_edge ;

%% Calculate percent signal ghosting (PSG)
PSG = 100 * abs( (Sbkg_left+Sbkg_right) - (Sbkg_btm + Sbkg_top) ) / (2*Sroi) ; % in percent

%% Calculate percentage Image Uniformity (PIU)
PIU = 100 * (1 - (smallROI_meanMax - smallROI_meanMin) / (smallROI_meanMax + smallROI_meanMax)) ;

%% Calculate signal-to-noise ratio based on NEMA standard
% Method 1: use 2 acquisitions for SNR calculcation
im_diff = im - squeeze(im_all(:,:,6,2)) ;
im_diff = im_diff .* largeROI_mask ;
im_diff(im_diff==0) = nan ;
im_diff_mean = mean(im_diff(:), 'omitnan') ;
sd1 = std(im_diff(:), 'omitnan') ;
SNR_2acq = sqrt(2) * Sroi / sd1 ;

% Method 2: use 1 acquisition for SNR calculcation
% bkg_noise = im .* (bkg_left + bkg_right + bkg_btm + bkg_top) ;
bkg_noise = im .* (bkg_left + bkg_right ) ;
bkg_noise(bkg_noise==0) = nan ;
bkg_noise_std = std(bkg_noise(:), 'omitnan') ;
SNR_1acq = sqrt((4-pi)/2) * Sroi / bkg_noise_std ;
end