clear all ; close all; clc ;
imfile = dir([pwd '/product_se_ice/*.nii']) ;
Ndata = size(imfile, 1) ; Nrecon = 4 ;
nx = 256 ; ny = 256 ; nslice = 11; nacq = 2;
im_all_4recon = zeros(nx,ny,nslice,nacq,Ndata,Nrecon) ;
PSG_all = zeros(Ndata, Nrecon) ;
PIU_all = zeros(Ndata, Nrecon) ;
SNR_1acq_all = zeros(Ndata, Nrecon) ;
SNR_2acq_all = zeros(Ndata, Nrecon) ;
load([pwd '/product_se_ice/results.mat']) ;
im_all_4recon(:,:,:,:,:,1) = im_all ;
PSG_all(:,1) = PSG ;
PIU_all(:,1) = PIU ;
SNR_1acq_all(:,1) = SNR_1acq ;
SNR_2acq_all(:,1) = SNR_2acq ;

load([pwd '/pulseq_se_ice/results.mat']) ;
im_all_4recon(:,:,:,:,:,2) = im_all ;
PSG_all(:,2) = PSG ;
PIU_all(:,2) = PIU ;
SNR_1acq_all(:,2) = SNR_1acq ;
SNR_2acq_all(:,2) = SNR_2acq ;

load([pwd '/product_se_gt/results.mat']) ;
im_all_4recon(:,:,:,:,:,3) = im_all ;
PSG_all(:,3) = PSG ;
PIU_all(:,3) = PIU ;
SNR_1acq_all(:,3) = SNR_1acq ;
SNR_2acq_all(:,3) = SNR_2acq ;
mask_all1 = mask_all ;

load([pwd '/pulseq_se_gt/results.mat']) ;
im_all_4recon(:,:,:,:,:,4) = im_all ;
PSG_all(:,4) = PSG ;
PIU_all(:,4) = PIU ;
SNR_1acq_all(:,4) = SNR_1acq ;
SNR_2acq_all(:,4) = SNR_2acq ;

SNR_ratio = SNR_1acq_all ./ SNR_2acq_all ;

%% statistical analysis
% mean and std
PSG_all_mean = 10*mean(PSG_all) ;
PSG_all_std = 10*std(PSG_all) ;
PSG_all_percent = PSG_all_std./PSG_all_mean*100 ;

PSG_all_mean_ice_percent = (PSG_all_mean(2) - PSG_all_mean(1)) / PSG_all_mean(1) * 100 ;
PSG_all_mean_gt_percent = (PSG_all_mean(4) - PSG_all_mean(3)) / PSG_all_mean(3) * 100 ;
PSG_all_mean_product_percent = (PSG_all_mean(3) - PSG_all_mean(1)) / PSG_all_mean(1) * 100 ;
PSG_all_mean_pulseq_percent = (PSG_all_mean(4) - PSG_all_mean(2)) / PSG_all_mean(2) * 100 ;

PSG_all_std_ice_percent = (PSG_all_percent(2) - PSG_all_percent(1)) ;
PSG_all_std_gt_percent = (PSG_all_percent(4) - PSG_all_percent(3)) ;
PSG_all_std_product_percent = (PSG_all_percent(3) - PSG_all_percent(1)) ;
PSG_all_std_pulseq_percent = (PSG_all_percent(4) - PSG_all_percent(2)) ;

PIU_all_mean = mean(PIU_all) ;
PIU_all_std = std(PIU_all) ;
PIU_all_percent = PIU_all_std./PIU_all_mean*100 ;

PIU_all_mean_ice_percent = (PIU_all_mean(2) - PIU_all_mean(1)) / PIU_all_mean(1) * 100 ;
PIU_all_mean_gt_percent = (PIU_all_mean(4) - PIU_all_mean(3)) / PIU_all_mean(3) * 100 ;
PIU_all_mean_product_percent = (PIU_all_mean(3) - PIU_all_mean(1)) / PIU_all_mean(1) * 100 ;
PIU_all_mean_pulseq_percent = (PIU_all_mean(4) - PIU_all_mean(2)) / PIU_all_mean(2) * 100 ;

PIU_all_std_ice_percent = (PIU_all_percent(2) - PIU_all_percent(1)) ;
PIU_all_std_gt_percent = (PIU_all_percent(4) - PIU_all_percent(3)) ;
PIU_all_std_product_percent = (PIU_all_percent(3) - PIU_all_percent(1)) ;
PIU_all_std_pulseq_percent = (PIU_all_percent(4) - PIU_all_percent(2)) ;

SNR_1acq_all_mean = mean(SNR_1acq_all) ;
SNR_1acq_all_std = std(SNR_1acq_all) ;
SNR_1acq_all_percent = SNR_1acq_all_std./SNR_1acq_all_mean*100 ;

SNR_1acq_all_mean_ice_percent = (SNR_1acq_all_mean(2) - SNR_1acq_all_mean(1)) / SNR_1acq_all_mean(1) * 100 ;
SNR_1acq_all_mean_gt_percent = (SNR_1acq_all_mean(4) - SNR_1acq_all_mean(3)) / SNR_1acq_all_mean(3) * 100 ;
SNR_1acq_all_mean_product_percent = (SNR_1acq_all_mean(3) - SNR_1acq_all_mean(1)) / SNR_1acq_all_mean(1) * 100 ;
SNR_1acq_all_mean_pulseq_percent = (SNR_1acq_all_mean(4) - SNR_1acq_all_mean(2)) / SNR_1acq_all_mean(2) * 100 ;

SNR_1acq_all_std_ice_percent = (SNR_1acq_all_percent(2) - SNR_1acq_all_percent(1)) ;
SNR_1acq_all_std_gt_percent = (SNR_1acq_all_percent(4) - SNR_1acq_all_percent(3)) ;
SNR_1acq_all_std_product_percent = (SNR_1acq_all_percent(3) - SNR_1acq_all_percent(1)) ;
SNR_1acq_all_std_pulseq_percent = (SNR_1acq_all_percent(4) - SNR_1acq_all_percent(2)) ;

SNR_2acq_all_mean = mean(SNR_2acq_all) ;
SNR_2acq_all_std = std(SNR_2acq_all) ;
SNR_2acq_all_percent = SNR_2acq_all_std./SNR_2acq_all_mean*100 ;

SNR_2acq_all_mean_ice_percent = (SNR_2acq_all_mean(2) - SNR_2acq_all_mean(1)) / SNR_2acq_all_mean(1) * 100 ;
SNR_2acq_all_mean_gt_percent = (SNR_2acq_all_mean(4) - SNR_2acq_all_mean(3)) / SNR_2acq_all_mean(3) * 100 ;
SNR_2acq_all_mean_product_percent = (SNR_2acq_all_mean(3) - SNR_2acq_all_mean(1)) / SNR_2acq_all_mean(1) * 100 ;
SNR_2acq_all_mean_pulseq_percent = (SNR_2acq_all_mean(4) - SNR_2acq_all_mean(2)) / SNR_2acq_all_mean(2) * 100 ;

SNR_2acq_all_std_ice_percent = (SNR_2acq_all_percent(2) - SNR_2acq_all_percent(1)) ;
SNR_2acq_all_std_gt_percent = (SNR_2acq_all_percent(4) - SNR_2acq_all_percent(3)) ;
SNR_2acq_all_std_product_percent = (SNR_2acq_all_percent(3) - SNR_2acq_all_percent(1)) ;
SNR_2acq_all_std_pulseq_percent = (SNR_2acq_all_percent(4) - SNR_2acq_all_percent(2)) ;

SNR_ratio_all_mean = mean(SNR_ratio) ;
SNR_ratio_all_std = std(SNR_ratio) ;
SNR_ratio_all_percent = round(SNR_ratio_all_std./SNR_ratio_all_mean*100) ;

SNR_ratio_all_mean_ice_percent = (SNR_ratio_all_mean(2) - SNR_ratio_all_mean(1)) / SNR_ratio_all_mean(1) * 100 ;
SNR_ratio_all_mean_gt_percent = (SNR_ratio_all_mean(4) - SNR_ratio_all_mean(3)) / SNR_ratio_all_mean(3) * 100 ;
SNR_ratio_all_mean_product_percent = (SNR_ratio_all_mean(3) - SNR_ratio_all_mean(1)) / SNR_ratio_all_mean(1) * 100 ;
SNR_ratio_all_mean_pulseq_percent = (SNR_ratio_all_mean(4) - SNR_ratio_all_mean(2)) / SNR_ratio_all_mean(2) * 100 ;

save cimax_metrics_se_mean_std PSG_all PIU_all SNR_1acq_all SNR_2acq_all PSG_all_mean_ice_percent PSG_all_mean_gt_percent PSG_all_mean_product_percent PSG_all_mean_pulseq_percent...
PSG_all_std_ice_percent PSG_all_std_gt_percent PSG_all_std_product_percent PSG_all_std_pulseq_percent...
PIU_all_mean_ice_percent PIU_all_mean_gt_percent PIU_all_mean_product_percent PIU_all_mean_pulseq_percent...
PIU_all_std_ice_percent PIU_all_std_gt_percent PIU_all_std_product_percent PIU_all_std_pulseq_percent...
SNR_1acq_all_mean_ice_percent SNR_1acq_all_mean_gt_percent SNR_1acq_all_mean_product_percent SNR_1acq_all_mean_pulseq_percent...
SNR_1acq_all_std_ice_percent SNR_1acq_all_std_gt_percent SNR_1acq_all_std_product_percent SNR_1acq_all_std_pulseq_percent...
SNR_2acq_all_mean_ice_percent SNR_2acq_all_mean_gt_percent SNR_2acq_all_mean_product_percent SNR_2acq_all_mean_pulseq_percent...
SNR_2acq_all_std_ice_percent SNR_2acq_all_std_gt_percent SNR_2acq_all_std_product_percent SNR_2acq_all_std_pulseq_percent ;

col = {'PSG (‰)', 'PIU (%)', 'SNR1', 'SNR2', 'SNR1/SNR2'} ;
row = {'Product+ICE', 'Pulseq+ICE', 'Product+GT', 'Pulseq+GT'} ;
data = {[num2str(PSG_all_mean(1), '%.2f') '±' num2str(PSG_all_std(1), '%.2f') ' (' num2str(round(PSG_all_percent(1))) '%)'], [num2str(PIU_all_mean(1), '%.1f') '±' num2str(PIU_all_std(1), '%.1f') ' (' num2str(round(PIU_all_percent(1))) '%)'],...
    [num2str(round(SNR_1acq_all_mean(1))) '±' num2str(round(SNR_1acq_all_std(1))) ' (' num2str(round(SNR_1acq_all_percent(1))) '%)'], [num2str(round(SNR_2acq_all_mean(1))) '±' num2str(round(SNR_2acq_all_std(1))) ' (' num2str(round(SNR_2acq_all_percent(1))) '%)'],...
    [num2str(SNR_ratio_all_mean(1), '%.2f') '±' num2str(SNR_ratio_all_std(1), '%.2f') ' (' num2str(round(SNR_ratio_all_percent(1))) '%)']; ...
    
    [num2str(PSG_all_mean(2), '%.2f') '±' num2str(PSG_all_std(2), '%.2f') ' (' num2str(round(PSG_all_percent(2))) '%)'], [num2str(PIU_all_mean(2), '%.1f') '±' num2str(PIU_all_std(2), '%.1f') ' (' num2str(round(PIU_all_percent(2))) '%)'],...
    [num2str(round(SNR_1acq_all_mean(2))) '±' num2str(round(SNR_1acq_all_std(2))) ' (' num2str(SNR_1acq_all_percent(2)) '%)'], [num2str(round(SNR_2acq_all_mean(2))) '±' num2str(round(SNR_2acq_all_std(2))) ' (' num2str(round(SNR_2acq_all_percent(2))) '%)'],...
    [num2str(SNR_ratio_all_mean(2), '%.2f') '±' num2str(SNR_ratio_all_std(2), '%.2f') ' (' num2str(round(SNR_ratio_all_percent(2))) '%)']; ...
    
    [num2str(PSG_all_mean(3), '%.2f') '±' num2str(PSG_all_std(3), '%.2f') ' (' num2str(round(PSG_all_percent(3))) '%)'], [num2str(PIU_all_mean(3), '%.1f') '±' num2str(PIU_all_std(3), '%.1f') ' (' num2str(round(PIU_all_percent(3))) '%)'],...
    [num2str(round(SNR_1acq_all_mean(3))) '±' num2str(round(SNR_1acq_all_std(3))) ' (' num2str(round(SNR_1acq_all_percent(3))) '%)'], [num2str(round(SNR_2acq_all_mean(3))) '±' num2str(round(SNR_2acq_all_std(3))) ' (' num2str(round(SNR_2acq_all_percent(3))) '%)'],...
    [num2str(SNR_ratio_all_mean(3), '%.2f') '±' num2str(SNR_ratio_all_std(3), '%.2f') ' (' num2str(round(SNR_ratio_all_percent(3))) '%)']; ...
    
    [num2str(PSG_all_mean(4), '%.2f') '±' num2str(PSG_all_std(4), '%.2f') ' (' num2str(round(PSG_all_percent(4))) '%)'], [num2str(PIU_all_mean(4), '%.1f') '±' num2str(PIU_all_std(4), '%.1f') ' (' num2str(round(PIU_all_percent(4))) '%)'],...
    [num2str(round(SNR_1acq_all_mean(4))) '±' num2str(round(SNR_1acq_all_std(4))) ' (' num2str(round(SNR_1acq_all_percent(4))) '%)'], [num2str(round(SNR_2acq_all_mean(4))) '±' num2str(round(SNR_2acq_all_std(4))) ' (' num2str(round(SNR_2acq_all_percent(4))) '%)'],...
    [num2str(SNR_ratio_all_mean(4), '%.2f') '±' num2str(SNR_ratio_all_std(4), '%.2f') ' (' num2str(round(SNR_ratio_all_percent(4))) '%)']; } ;
h = figure ;
u = uitable('columnname', col, 'rowname', row,'Position', [20 20 500 70], 'data', data) ;
set(u, 'ColumnWidth', {150,150,150,150,150}) ;
table_extent = get(u,'Extent');
set(u,'Position',[1 1 table_extent(3) table_extent(4)])
figure_size = get(h,'outerposition');
desired_fig_size = [figure_size(1) figure_size(2) table_extent(3) table_extent(4)+80];
set(h,'outerposition', desired_fig_size+20);
% set(gcf, 'Position', get(0, 'Screensize'));
print('-dpng', 'Structural_quality_metrics', '-r0'); 

[p_psg_ice, h_psg_ice] = ranksum(PSG_all(:,1), PSG_all(:,2)) ;
[p_psg_gt, h_psg_gt] = ranksum(PSG_all(:,3), PSG_all(:,4)) ;
[p_psg_product, h_psg_product] = ranksum(PSG_all(:,1), PSG_all(:,3)) ;
[p_psg_pulseq, h_psg_pulseq] = ranksum(PSG_all(:,2), PSG_all(:,4)) ;

[p_piu_ice, h_piu_ice] = ranksum(PIU_all(:,1), PIU_all(:,2)) ;
[p_piu_gt, h_piu_gt] = ranksum(PIU_all(:,3), PIU_all(:,4)) ;
[p_piu_product, h_piu_product] = ranksum(PIU_all(:,1), PIU_all(:,3)) ;
[p_piu_pulseq, h_piu_pulseq] = ranksum(PIU_all(:,2), PIU_all(:,4)) ;

[p_snr1_ice, h_snr1_ice] = ranksum(SNR_1acq_all(:,1), SNR_1acq_all(:,2)) ;
[p_snr1_gt, h_snr1_gt] = ranksum(SNR_1acq_all(:,3), SNR_1acq_all(:,4)) ;
[p_snr1_product, h_snr1_product] = ranksum(SNR_1acq_all(:,1), SNR_1acq_all(:,3)) ;
[p_snr1_pulseq, h_snr1_pulseq] = ranksum(SNR_1acq_all(:,2), SNR_1acq_all(:,4)) ;

[p_snr2_ice, h_snr2_ice] = ranksum(SNR_2acq_all(:,1), SNR_2acq_all(:,2)) ;
[p_snr2_gt, h_snr2_gt] = ranksum(SNR_2acq_all(:,3), SNR_2acq_all(:,4)) ;
[p_snr2_product, h_snr2_product] = ranksum(SNR_2acq_all(:,1), SNR_2acq_all(:,3)) ;
[p_snr2_pulseq, h_snr2_pulseq] = ranksum(SNR_2acq_all(:,2), SNR_2acq_all(:,4)) ;

%%
figure ;
subplot (2,2,1) ;
boxplot(PSG_all, 'Labels', {'product+ICE', 'Pulseq+ICE', 'Product+GT', 'Pulseq+GT'}) ;
title('Percent signal ghosting') ; % meanI = mean(sub(:));
subplot (2,2,2) ;
boxplot(PIU_all, 'Labels', {'product+ICE', 'Pulseq+ICE', 'Product+GT', 'Pulseq+GT'}) ;
title('Percent intensity uniformity') ; % STD of detrended signal residual sd=std(avg_signal_roi - yfit);
subplot (2,2,3) ;
boxplot(SNR_1acq_all, 'Labels', {'product+ICE', 'Pulseq+ICE', 'Product+GT', 'Pulseq+GT'}) ;
title ('SNR (1 acq.)') ; % snr = meanI/sqrt(varI/N); varI = var(Idiff_noise(:)*(fix(N/2)));         %Variance of Isub within ROI only
subplot (2,2,4) ;
boxplot(SNR_2acq_all, 'Labels', {'product+ICE', 'Pulseq+ICE', 'Product+GT', 'Pulseq+GT'}) ;
title ('SNR (2 acq.)') ; % sfnr = Iave./(Isd + eps); sfnr_mean = mean(sfnr)

figure ;
% subplot (2,2,1) ;
boxplot(SNR_ratio, 'Labels', {'product+ICE', 'Pulseq+ICE', 'Product+GT', 'Pulseq+GT'}) ;
title('SNR ratio') ; % meanI = mean(sub(:));

%% plot signal image
im = squeeze(im_all_4recon(:,:,6,1,round(Ndata/2),:)) ;
% I_min = min(im(:)) ;
% I_max = max(im(:)) ;
lineWidth = 2 ;fontSize = 24 ;
figure ;
t = tiledlayout(1, 4) ;
nexttile
imshow(rot90(im(:,:,1)),[]) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
% axis normal
% axis on
set(ax,'XTick',[], 'YTick', [],'FontSize', fontSize, 'fontWeight', 'bold');
colormap gray ;
% colorbar ;
% clim([I_min, I_max]) ;

nexttile
imshow(rot90(im(:,:,2)),[]) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
% axis normal
% axis on
set(ax,'XTick',[], 'YTick', [],'FontSize', fontSize, 'fontWeight', 'bold');
colormap gray ;
% colorbar ;
% clim([I_min, I_max]) ;

nexttile
imshow(rot90(im(:,:,3)),[]) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
% axis normal
% axis on
set(ax,'XTick',[], 'YTick', [],'FontSize', fontSize, 'fontWeight', 'bold');
colormap gray ;
% colorbar ;
% clim([I_min, I_max]) ;

nexttile
imshow(rot90(im(:,:,4)),[]) ;
hold on ;
% line([Nx Nx],[0 Ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
% axis normal
% axis on
set(ax,'XTick',[], 'YTick', [],'FontSize', fontSize, 'fontWeight', 'bold');
colormap gray ;
% colorbar ;
% clim([I_min, I_max]) ;
% cbar = colorbar ;
% cbar.Label.String = 'Intensity';
t.Padding = 'none' ;
t.TileSpacing = 'none' ;
set(gcf, 'Position', get(0, 'Screensize'));
exportgraphics(gcf,'se_image.png') ;

%%
im_diff = squeeze(im_all_4recon(:,:,6,1,round(Ndata/2),:) - im_all_4recon(:,:,6,2,round(Ndata/2),:)) ;
figure ;
t = tiledlayout(1, 4) ;
nexttile
imshow(rot90(im_diff(:,:,1)),[]) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
% axis normal
% axis on
set(ax,'XTick',[], 'YTick', [],'FontSize', fontSize, 'fontWeight', 'bold');
colormap gray ;
% colorbar ;
% clim([I_min, I_max]) ;

nexttile
imshow(rot90(im_diff(:,:,2)),[]) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
% axis normal
% axis on
set(ax,'XTick',[], 'YTick', [],'FontSize', fontSize, 'fontWeight', 'bold');
colormap gray ;
% colorbar ;
% clim([I_min, I_max]) ;

nexttile
imshow(rot90(im_diff(:,:,3)),[]) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
% axis normal
% axis on
set(ax,'XTick',[], 'YTick', [],'FontSize', fontSize, 'fontWeight', 'bold');
colormap gray ;
% colorbar ;
% clim([I_min, I_max]) ;

nexttile
imshow(rot90(im_diff(:,:,4)),[]) ;
hold on ;
% line([Nx Nx],[0 Ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
% axis normal
% axis on
set(ax,'XTick',[], 'YTick', [],'FontSize', fontSize, 'fontWeight', 'bold');
colormap gray ;
% colorbar ;
% clim([I_min, I_max]) ;
% cbar = colorbar ;
% cbar.Label.String = 'Intensity';
t.Padding = 'none' ;
t.TileSpacing = 'none' ;
set(gcf, 'Position', get(0, 'Screensize'));
exportgraphics(gcf,'se_image_diff.png') ;

%% mask plotting

magFactor = 400 ;
allmask = mask_all.bkg_left+mask_all.bkg_right+mask_all.bkg_top+mask_all.bkg_btm ;
green = zeros(size(allmask,1), size(allmask,2), 3) ;
green(:,:,2) = 0.6 ;
red = zeros(size(green)) ;
red(:,:,1) = 0.8 ;
blue = zeros(size(green)) ;
blue(:,:,3) = 1 ;
yellow = zeros(size(green)) ;
yellow(:,:,2) = 1 ;
yellow(:,:,1) = 1 ;
figure;
imshow(rot90(im_all_4recon(:,:,6,1,end)) , [],'InitialMagnification', magFactor) ;
colormap(gray);
hold all;
h=imshow(blue,'InitialMagnification', magFactor) ;
set(h,'AlphaData', rot90(edge(mask_all.largeROI_mask))) ;
h=imshow(red, 'InitialMagnification', magFactor) ;
set(h,'AlphaData', rot90(allmask)) ;
h=imshow(green, 'InitialMagnification', magFactor) ;
set(h,'AlphaData', rot90(mask_all.smallROI_min_mask)) ;
h=imshow(yellow, 'InitialMagnification', magFactor) ;
set(h,'AlphaData', rot90(mask_all.smallROI_max_mask)) ;
% set(gcf, 'Position', get(0, 'Screensize'));
exportgraphics(gcf, 'se_mask.png') ;

