%% This is the temporal quality analysis for the quality assurance protocol.
% transversal multi-slice EPI sequence with ramp smapling and navigator:
% write_QA_Tran_EPIrs.m
% clear all ; close all; clc ;
magFactor = 100 ;
lineWidth = 1.5 ;
fontSize = 24 ;

FOVx = 220 ; FOVy = 220 ;
%% load images
imfile = dir('*.nii') ;
im = double(niftiread(imfile(1).name)) ; % [Nx, Ny, Nslice, Nrun]

centerSlc = ceil(size(im, 3) / 2) ;
Nx = size(im, 1) ; Ny = size(im, 2) ; % matrix size
Nrun_discard = 2 ;
Nrun = size(im, 4) - Nrun_discard ;
centerRun = round(Nrun/2) ;
% The first two volumes are discarded to allow NMR and eddy0current
% equilibrium to be achieved.
im = squeeze(im(:, :, centerSlc, (Nrun_discard+1):(Nrun_discard+Nrun)) ) ;

%% mask making
% detect edge of the 2D epi image
im_edge = double(edge(im(:,:,centerRun)) ) ;

[im_X, im_Y] = ind2sub(size(im_edge), find (im_edge)) ;
% fit for the center and radius of the circular phantom
[Xc, Yc, R] = circfit(im_X, im_Y) ;
phantom_mask = makeCircleMask([Nx, Ny], Xc, Yc, R) ;
Xc = round(Xc) ;
Yc = round(Yc) ;
R = round(R) ;
ROImask = zeros(Nx, Ny) ;
% 21*21 voxel ROI centered in the phantom for matrix size = 64*64
ROI_halfwidthX = round(10*Nx/64) ;
ROI_halfwidthY = round(10*Ny/64) ;
ROImask(Xc-ROI_halfwidthX:Xc+ROI_halfwidthX, Yc-ROI_halfwidthY:Yc+ROI_halfwidthY) = 1 ;

figure ;
t = tiledlayout(1, 3) ;
nexttile
im_edge_fuse = labeloverlay(mat2gray(im(:,:,centerRun)), im_edge,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_edge_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([Nx Nx],[0 Ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Edge mask','FontSize', fontSize) ;
nexttile
im_phantom_mask_fuse = labeloverlay(mat2gray(im(:,:,centerRun)), phantom_mask,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_phantom_mask_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([Nx Nx],[0 Ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Phantom mask','FontSize', fontSize) ;
nexttile
im_largeROI_mask_fuse = labeloverlay(mat2gray(im(:,:,centerRun)), ROImask,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_largeROI_mask_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('ROI mask','FontSize', fontSize) ;
t.Padding = 'none' ;
t.TileSpacing = 'none' ;

%% generate two background masks in the readout direction
pxsize = FOVx/Nx * FOVy/Ny ;
ROI_area = 10 * 100 ; % backgorund ROI area = 10 cm2
R1 = ceil(R * 0.1) ;
% left background ROI
% Wl = size(R1:Xc-R-R1+1, 2) ;
Wl = size(R1:Xc-R-R1, 2) ;
Hl2 = round(ROI_area / Wl / pxsize / 2) ;
bkg_left = zeros(Nx,Ny) ;
% bkg_left(R1:Xc-R-R1+1,Yc-Hl2+1:Yc+Hl2) = 1 ;
bkg_left(R1:Xc-R-R1,Yc-Hl2+1:Yc+Hl2) = 1 ;
ROIbkg_left = bkg_left .* im ;
Sbkg_left = sum(ROIbkg_left(:)) / nnz(bkg_left) ;
ROIbkg_left_std = std(ROIbkg_left(:) ) ;

% right background ROI
Wr = size(Xc+R+R1-1:Nx-R1+1, 2) ;
%Wr = size(Xc+R+R1:Nx-R1+1, 2) ;
Hr2 = round(ROI_area / Wr / pxsize / 2) ;
bkg_right = zeros(Nx,Ny) ;
% bkg_right(Xc+R+R1-1:Nx-R1+1,Yc-Hr2+1:Yc+Hr2) = 1 ;
bkg_right(Xc+R+R1:Nx-R1+1,Yc-Hr2+1:Yc+Hr2) = 1 ;
ROIbkg_right = bkg_right .* im ;
Sbkg_right = sum(ROIbkg_right(:)) / nnz(bkg_right) ;
ROIbkg_right_std = std(ROIbkg_right(:) ) ;

figure ;
t = tiledlayout(1, 2) ;
nexttile
im_bkg_left_fuse = labeloverlay(mat2gray(im(:,:,centerRun)), bkg_left,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_bkg_left_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([Nx Nx],[0 Ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Left bkg ROI','FontSize', fontSize) ;
nexttile
im_bkg_right_fuse = labeloverlay(mat2gray(im(:,:,centerRun)), bkg_right,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_bkg_right_fuse),[],'InitialMagnification', magFactor) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Right bkg ROI','FontSize', fontSize) ;
t.Padding = 'none';
t.TileSpacing = 'none';

%% create background ROI matrix: roi_bkg
clear roi_bkg ;
for i = 1:Nrun
    bkg = [bkg_left; bkg_right] ;
    im_bkg = [im(:,:,i) .* bkg_left; im(:,:,i) .* bkg_right] ; 
    im_bkg_nonzero = im_bkg(bkg~=0) ;
    roi_bkg(i,:) = im_bkg_nonzero(:) ;
end

%% create signal ROI matrix: roi_sig
clear roi_sig ;
numvox = (2*ROI_halfwidthX+1) * (2*ROI_halfwidthY+1) ;
roi_sig = zeros(Nrun, numvox) ;
im_ROI = im .* ROImask ;
for i = 1:Nrun
    im_ROI_temp = im((Xc-ROI_halfwidthX):(Xc+ROI_halfwidthX), (Yc-ROI_halfwidthY):(Yc+ROI_halfwidthY), i) ;
    roi_sig(i,:) = im_ROI_temp(:) ;
end

%% produce signal image: the simple average, voxel by voxel, across the 198 volumes.
signalImage = mean(im, 3) ;
figure ;
imshow(rot90(signalImage), [],'InitialMagnification', magFactor*4) ;
colormap default ;
cbar = colorbar ;
set(cbar,'FontSize',18, 'FontWeight', 'bold');
title('signal image') ;
%%   Signal-to-Fluctuation-Noise Ratio (SFNR)
%    = mean voxel signal across time divided by temporal standard deviation
%      of residuals obtained by fitting second-order polynomial to voxel
%      signal data

disp('Calculating SFNR using quadratic detrending.')
sig_detrended = roi_sig ; % initialize
x = (1:Nrun)' ;
for i = 1:numvox
    p = polyfit(x, roi_sig(:,i), 2) ;
    predicted = polyval(p, x) ;
    sig_detrended(:, i) = roi_sig(:,i) - predicted ;
end
noise_SD = std(sig_detrended) ;
figure ;
plot(predicted, 'r') ;
hold on ;
plot(roi_sig(:,end),'g') ;
legend('second-order polynomial trend', 'before detrended') ;
xlabel('#volume') ;
ylabel('Intensity') ;
mean_sig = mean(roi_sig,1) ;
mean_sig_mean = mean(mean_sig) ;
% Avoid Inf's or extreme outliers, in case signal ROI included any zero voxels:
if any(noise_SD<.001)
    ind = noise_SD<.001;
    noise_SD(ind) = nan;
    mean_sig(ind) = nan;
end
% Calculate Local:
SFNR = mean_sig ./ noise_SD;
funct_SFNR = mean(SFNR, 'omitnan') ;
disp(['Signal-to-Fluctuation-Noise Ratio (SFNR) = ', num2str(funct_SFNR)]) ;

% detrending...
im_detrended = zeros(size(im_ROI)) ;
im_detrended_std = zeros(Nx, Ny) ;
for i = 1:Nx
    for j = 1:Ny
        p = polyfit(x, squeeze(im(i,j,:)), 2) ;
        predicted = polyval(p, x) ;
        im_detrended(i,j,:) = squeeze(im(i,j,:)) - predicted ;
        im_detrended_std(i,j) = std(im_detrended(i,j,:)) ;
    end
end
% Avoid Inf's or extreme outliers, in case signal ROI included any zero voxels
if any(im_detrended_std < 0.001)
    ind = im_detrended_std < 0.001 ;
    im_detrended_std(ind) = nan ;
    im_detrended(ind) = nan ;
end
im_SFNR = signalImage ./ im_detrended_std ;
im_SFNR(isinf(im_SFNR)) = nan ;
im_SFNR(isnan(im_SFNR)) = nan ;
im_SFNR(im_SFNR<0) = nan ;
im_SFNR(im_SFNR>9999) = nan ;
im_SFNR_ROI = im_SFNR .* ROImask ;
funct_SFNR1 = sum(im_SFNR_ROI(:), 'omitnan') / sum(ROImask(:)) ; % should be the same as SFNR
% disp(['Signal-to-Fluctuation-Noise Ratio (SFNR1) = ', num2str(funct_SFNR1)]) ;

%% Temporal Signal-to-Noise Ratio (tSNR)
%       = mean voxel signal across time divided by SD across time
%       *References: 
%                  Triantafyllou, C., Hoge, R. D., Krueger, G., Wiggins, 
%                   C. J., Potthast, A., Wiggins, G. C., & Wald, L. L. (2005). 
%                   Comparison of physiological noise at 1.5 T, 3 T and 7 T 
%                   and optimization of fMRI acquisition parameters. 
%                   Neuroimage, 26(1), 243-250.
%                  Murphy, K., Bodurka, J., & Bandettini, P. A. (2007). 
%                   How long to scan? The relationship between fMRI temporal 
%                   signal to noise ratio and necessary scan duration. 
%                   Neuroimage, 34(2), 565-574. 
%                   <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2223273/>
noise_SD = std(roi_sig) ;
mean_sig = mean(roi_sig) ;
% Avoid Inf's or extreme outliers, in case signal ROI included any zero voxels:
if any(noise_SD < 0.01)
    ind = noise_SD < 0.01 ;
    noise_SD(ind) = nan ;
    mean_sig(ind) = nan ;
end

% Calculate for Signal ROI:
tSNR_ROI = mean_sig ./ noise_SD ;
funct_tSNR = mean(tSNR_ROI, 'omitnan') ;
disp(['Temporal Signal-to-Noise Ratio (tSNR) = ',num2str(funct_tSNR)]) ;

% Set NaN's, Inf's, and out of range values (negatives) to 0:
im_tSNR = mean(im, 3, 'omitnan') ./ std(im, 0, 3, 'omitnan') ;
im_tSNR(isinf(im_tSNR)) = 0 ;
im_tSNR(isnan(im_tSNR)) = 0 ;
im_tSNR(im_tSNR<0) = 0 ;
im_tSNR(im_tSNR>9999) = 0 ;

im_tSNR_ROI = im_tSNR .* ROImask ;
funct_tSNR1 = sum(im_tSNR_ROI(:), 'omitnan') / (sum(ROImask(:)) - sum(isnan(im_tSNR_ROI(:)))) ;
% disp(['Temporal Signal-to-Noise Ratio (tSNR1) = ',num2str(funct_tSNR1)]) ;


%%   Temporal Signal-to-Background Noise Ratio (tSBNR)
%       = spatial average across a signal ROI of temporal mean voxel
%       signals / spatial average across temporal standard deviation of 
%       noise/background ROI voxels
%       -whereas tSNR divides mean signal by the temporal SD of the same
%       voxels, this measure compares the mean signal to the average
%       temporal SD of background/noise voxels, which may provide a more
%       accurate measure of temporal fluctuations due specifically to noise

noise_SD = std(roi_bkg) ; % temporal standard deviation of background ROI voxels
mean_sig = mean(roi_sig) ; % temporal mean voxel signals
% Avoid Inf's or extreme outliers, in case signal ROI included any zero voxels:
if any(noise_SD < 0.0001)
    ind = noise_SD < 0.0001 ;
    noise_SD(ind) = nan ;
end
noise_SD = mean(noise_SD, 'omitnan') ;
tSBNR = mean(im, 3) ./ noise_SD ;

% Calculate for Signal ROI
tSBNR_ROI = mean_sig ./ noise_SD ;
funct_tSBNR = mean(tSBNR_ROI, 'omitnan') ;
disp(['Temporal Signal-to-Background Noise Ratio (tSBNR) = ',num2str(funct_tSBNR)]) ;

% Set NaN's, Inf's, and other out of range values to 0;
tSBNR(isinf(tSBNR)) = 0 ;
tSBNR(isnan(tSBNR)) = 0 ;
tSBNR(tSBNR<0) = 0 ;
tSBNR(tSBNR>9999) = 0 ;
im_tSBNR_ROI = tSBNR .* ROImask ;
funct_tSBNR1 = sum(im_tSBNR_ROI(:), 'omitnan') / sum(ROImask(:)) ;
% disp(['Temporal Signal-to-Background Noise Ratio (tSBNR1) = ',num2str(funct_tSBNR1)]) ;

%% Signal-to-Noise Ratio (SNR-functional)
%       = mean voxel signal of signal ROI (averaged across time) divided by
%         the spatial SD of a noise ROI signal (averaged across time),
%         multiplied by a correction factor
%           Reference:
%               http://wagerlab.colorado.edu/wiki/doku.php/help/fmri_quality_control_overview

noise_SD = std(roi_bkg, 0, 2) ;
denominator = mean(noise_SD(:)) ; % sqrt(2/(4-pi))*noise_SD;
im_funct_SNR = signalImage ./ denominator ;
signalROI = signalImage .* ROImask ;
signalROI(ROImask==0) = nan ;
funct_SNR_sd = mean(signalROI(:), 'omitnan') ./ denominator ;
disp(['Signal-to-Noise Ratio (SNR-functional) = ',num2str(funct_SNR_sd)]);

%% plot metric images
figure ;
t = tiledlayout(2, 3) ;

nexttile
imshow(rot90(im_detrended_std),[],'InitialMagnification', magFactor) ;
hold on ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('temporal fluctuation noise image','FontSize', fontSize) ;
cbar = colorbar ;
set(cbar,'FontSize',18, 'FontWeight', 'bold');

nexttile
imshow(rot90(im_SFNR),[],'InitialMagnification', magFactor) ;
hold on ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('SFNR','FontSize', fontSize) ;
cbar = colorbar ;
set(cbar,'FontSize',18, 'FontWeight', 'bold');

nexttile
imshow(rot90(im_tSNR),[],'InitialMagnification', magFactor) ;
hold on ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('tSNR','FontSize', fontSize) ;
cbar = colorbar ;
set(cbar,'FontSize',18, 'FontWeight', 'bold');

nexttile
imshow(rot90(tSBNR),[],'InitialMagnification', magFactor) ;
hold on ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('tSBNR','FontSize', fontSize) ;
cbar = colorbar ;
set(cbar,'FontSize',18, 'FontWeight', 'bold');

nexttile
imshow(rot90(im_funct_SNR),[],'InitialMagnification', magFactor) ;
hold on ;
ax = gca;
set(ax,'XTick',[], 'YTick', []) ;
colormap default ;
title('SNR-functional','FontSize', fontSize) ;
cbar = colorbar ;
set(cbar,'FontSize',18, 'FontWeight', 'bold');

%% Static Spatial noise image, Friedman et al, MRM, 2006
sumODD = sum(im(:,:,1:2:end), 3) ;
sumEVEN = sum(im(:,:,2:2:end), 3) ;
DIFF = sumODD - sumEVEN ;
DIFF_ROI = DIFF .* ROImask ;
DIFF_ROI(ROImask==0) = nan ;
varianceSummaryVallue = var(DIFF_ROI(:), 'omitnan') ;
figure ;
imshow(rot90(DIFF), []) ;
title('DIFF=sumODD-sumEVEN') ;

%% SNR Summary value, Friedman et al, MRM, 2006
signalSummaryValue = mean(signalROI(:), 'omitnan') ;
SNRsummaryvalue = signalSummaryValue / sqrt(varianceSummaryVallue/Nrun ) ;
disp(['SNRsummaryvalue = ', num2str(SNRsummaryvalue)]) ;

%% Percent fluctuation and percent drift
spatialAverageIntensity = mean(roi_sig, 2) ;
spatialAverageIntensity_mean = mean(spatialAverageIntensity) ;
p = polyfit(x, squeeze(spatialAverageIntensity), 2) ;
trend = polyval(p, x) ;
residuals = spatialAverageIntensity - trend ;
residuals_std = std(residuals) ;
percentFluactuation = 100 * residuals_std / spatialAverageIntensity_mean ;
percentDrift = 100 * (max(trend) - min(trend)) / spatialAverageIntensity_mean ;
figure ;
plot(spatialAverageIntensity) ;
hold on ;
plot(trend) ;
legend('spatialAverageIntensity', 'trend') ;
xlabel('#volume') ;
ylabel('Intensity') ;
disp(['percentFluactuation = ', num2str(percentFluactuation)]) ;
disp(['percentDrift = ', num2str(percentDrift)]) ;

%% Fourier analysis of the residuals
TR = 2.99 ; % [sec]
L = length(residuals) ;
TA = L * TR ;
residuals_fft = fftshift(fft(fftshift(residuals))) ; % forward transform, Rader's / Mixed-Radix FFT
residuals_fft_p2 = abs(residuals_fft(L/2+1:end)) ;
w_axis = (1:(L/2))/TA ;
figure ;
plot(w_axis, residuals_fft_p2) ;
title('spectrum of residuals') ;
xlabel('Frequency (Hz)') ;
ylabel('Magnitude') ;

%% Weisskoff analysis: to be continued...
% coefficient of variation (CV, the SD of a time-series divided by the mean of the time-series)
% log(CV) vs. log(N)
% N = 21 ;
% CV = zeros(N, 1) ;
% 
% for i = 1:N
%     ROImask_wk = zeros(nx, ny) ;
%     ROImask_wk(Xc-i:Xc+i, Yc-i:Yc+i) = 1 ;
%     im_wk = im .* ROImask_wk ;
%     mean_wk = mean(im_wk, 3) ;
%     SD_wk = im_detrended_std * ROImask_wk ;
%     temp = 100 * mean_wk ./ SD_wk ;
%     temp(isnan(temp)) = 0 ;
%     CV(i) = sum(temp(:)) / sum(ROImask_wk(:)) ;
% end
% RDC = CV(1) - CV(end) ;
% figure ;
% plot(log(1:N), log(CV), 'k-') ;
% hold on ;
% plot(log(1:N), log(CV), 'r.') ;
% xlim([1 100]) ;
% ylim([1 10]) ;
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')