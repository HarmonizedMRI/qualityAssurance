clear all ; close all; clc ;
imfile = dir([pwd '/product_epi_ice/*.nii']) ;
Ndata = size(imfile, 1) ;
Nx = 64 ; Ny = 64 ;
Nrecon = 4 ; ROIsize = 15 ;
Idiff_noise_all = zeros(Nx, Ny, Ndata, Nrecon) ;
Iave_all = zeros(Nx, Ny, Ndata, Nrecon) ;
Isd_all = zeros(Nx, Ny, Ndata, Nrecon) ;
sfnr_all = zeros(Nx, Ny, Ndata, Nrecon) ;
nf_all = fix(198/2+1) ;
spectrum_all = zeros(nf_all, Ndata, Nrecon) ;
F_weis_all = zeros(ROIsize, Ndata, Nrecon) ;
meanI_all = zeros(Ndata, Nrecon) ;
fluct_sd_all = zeros(Ndata, Nrecon) ;
snr_all = zeros(Ndata, Nrecon) ;
sfnr_mean_all = zeros(Ndata, Nrecon) ;
rms_all = zeros(Ndata, Nrecon) ;
temp_drift_all = zeros(Ndata, Nrecon) ;
temp_drift_per_minute_all = zeros(Ndata, Nrecon) ;
max_temp_drift_all = zeros(Ndata, Nrecon) ;
rdc_all = zeros(Ndata, Nrecon) ;

load([pwd '/product_epi_ice/results.mat']) ;
Idiff_noise_all(:,:,:,1) = Idiff_noise ;
Iave_all(:,:,:,1) = Iave ;
Isd_all(:,:,:,1) = Isd ;
sfnr_all(:,:,:,1) = sfnr ;
spectrum_all(:,:,1) = spectrum ;
F_weis_all(:,:,1) = F_weis ;
meanI_all(:,1) = meanI ;
fluct_sd_all(:,1) = fluct_sd ;
snr_all(:,1) = snr ;
sfnr_mean_all(:,1) = sfnr_mean ;
rms_all(:,1) = rms ;
temp_drift_all(:,1) = temp_drift ;
temp_drift_per_minute_all(:,1) = temp_drift_per_minute ;
max_temp_drift_all(:,1) = max_temp_drift ;
rdc_all(:,1) = rdc ;
mask{1} = mask_all;

load([pwd '/pulseq_epi_ice/results.mat']) ;
Idiff_noise_all(:,:,:,2) = Idiff_noise ;
Iave_all(:,:,:,2) = Iave ;
Isd_all(:,:,:,2) = Isd ;
sfnr_all(:,:,:,2) = sfnr ;
spectrum_all(:,:,2) = spectrum ;
F_weis_all(:,:,2) = F_weis ;
meanI_all(:,2) = meanI ;
fluct_sd_all(:,2) = fluct_sd ;
snr_all(:,2) = snr ;
sfnr_mean_all(:,2) = sfnr_mean ;
rms_all(:,2) = rms ;
temp_drift_all(:,2) = temp_drift ;
temp_drift_per_minute_all(:,2) = temp_drift_per_minute ;
max_temp_drift_all(:,2) = max_temp_drift ;
rdc_all(:,2) = rdc ;
mask{2} = mask_all;

load([pwd '/product_epi_gt/results.mat']) ;
Idiff_noise_all(:,:,:,3) = Idiff_noise ;
Iave_all(:,:,:,3) = Iave ;
Isd_all(:,:,:,3) = Isd ;
sfnr_all(:,:,:,3) = sfnr ;
spectrum_all(:,:,3) = spectrum ;
F_weis_all(:,:,3) = F_weis ;
meanI_all(:,3) = meanI ;
fluct_sd_all(:,3) = fluct_sd ;
snr_all(:,3) = snr ;
sfnr_mean_all(:,3) = sfnr_mean ;
rms_all(:,3) = rms ;
temp_drift_all(:,3) = temp_drift ;
temp_drift_per_minute_all(:,3) = temp_drift_per_minute ;
max_temp_drift_all(:,3) = max_temp_drift ;
rdc_all(:,3) = rdc ;
mask{3} = mask_all;

load([pwd '/pulseq_epi_gt/results.mat']) ;
Idiff_noise_all(:,:,:,4) = Idiff_noise ;
Iave_all(:,:,:,4) = Iave ;
Isd_all(:,:,:,4) = Isd ;
sfnr_all(:,:,:,4) = sfnr ;
spectrum_all(:,:,4) = spectrum ;
F_weis_all(:,:,4) = F_weis ;
meanI_all(:,4) = meanI ;
fluct_sd_all(:,4) = fluct_sd ;
snr_all(:,4) = snr ;
sfnr_mean_all(:,4) = sfnr_mean ;
rms_all(:,4) = rms ;
temp_drift_all(:,4) = temp_drift ;
temp_drift_per_minute_all(:,4) = temp_drift_per_minute ;
max_temp_drift_all(:,4) = max_temp_drift ;
rdc_all(:,4) = rdc ;
mask{4} = mask_all;

%% statistical analysis
% mean and std
sfnr_mean_all_mean = mean(sfnr_mean_all) ;
sfnr_mean_all_std = std(sfnr_mean_all) ;
sfnr_mean_all_percent = sfnr_mean_all_std./sfnr_mean_all_mean*100 ;
sfnr_mean_all_mean_ice_percent = (sfnr_mean_all_mean(2) - sfnr_mean_all_mean(1)) / sfnr_mean_all_mean(1) * 100 ;
sfnr_mean_all_mean_gt_percent = (sfnr_mean_all_mean(4) - sfnr_mean_all_mean(3)) / sfnr_mean_all_mean(3) * 100 ;
sfnr_mean_all_mean_product_percent = (sfnr_mean_all_mean(3) - sfnr_mean_all_mean(1)) / sfnr_mean_all_mean(1) * 100 ;
sfnr_mean_all_mean_pulseq_percent = (sfnr_mean_all_mean(4) - sfnr_mean_all_mean(2)) / sfnr_mean_all_mean(2) * 100 ;

sfnr_mean_all_std_ice_percent = (sfnr_mean_all_percent(2) - sfnr_mean_all_percent(1))  ;
sfnr_mean_all_std_gt_percent = (sfnr_mean_all_percent(4) - sfnr_mean_all_percent(3))  ;
sfnr_mean_all_std_product_percent = (sfnr_mean_all_percent(3) - sfnr_mean_all_percent(1))  ;
sfnr_mean_all_std_pulseq_percent = (sfnr_mean_all_percent(4) - sfnr_mean_all_percent(2))  ;

rdc_all_mean = mean(rdc_all) ;
rdc_all_std = std(rdc_all) ;
rdc_all_percent = rdc_all_std./rdc_all_mean*100 ;
rdc_all_mean_ice_percent = (rdc_all_mean(2) - rdc_all_mean(1)) / rdc_all_mean(1) * 100 ;
rdc_all_mean_gt_percent = (rdc_all_mean(4) - rdc_all_mean(3)) / rdc_all_mean(3) * 100 ;
rdc_all_mean_product_percent = (rdc_all_mean(3) - rdc_all_mean(1)) / rdc_all_mean(1) * 100 ;
rdc_all_mean_pulseq_percent = (rdc_all_mean(4) - rdc_all_mean(2)) / rdc_all_mean(2) * 100 ;

rdc_all_std_ice_percent = (rdc_all_percent(2) - rdc_all_percent(1)) ;
rdc_all_std_gt_percent = (rdc_all_percent(4) - rdc_all_percent(3)) ;
rdc_all_std_product_percent = (rdc_all_percent(3) - rdc_all_percent(1)) ;
rdc_all_std_pulseq_percent = (rdc_all_percent(4) - rdc_all_percent(2)) ;

rms_all_mean = 10*mean(rms_all) ;
rms_all_std = 10*std(rms_all) ;
rms_all_percent = rms_all_std./rms_all_mean*100 ;

rms_all_mean_ice_percent = (rms_all_mean(2) - rms_all_mean(1)) / rms_all_mean(1) * 100 ;
rms_all_mean_gt_percent = (rms_all_mean(4) - rms_all_mean(3)) / rms_all_mean(3) * 100 ;
rms_all_mean_product_percent = (rms_all_mean(3) - rms_all_mean(1)) / rms_all_mean(1) * 100 ;
rms_all_mean_pulseq_percent = (rms_all_mean(4) - rms_all_mean(2)) / rms_all_mean(2) * 100 ;

rms_all_std_ice_percent = (rms_all_percent(2) - rms_all_percent(1)) ;
rms_all_std_gt_percent = (rms_all_percent(4) - rms_all_percent(3)) ;
rms_all_std_product_percent = (rms_all_percent(3) - rms_all_percent(1)) ;
rms_all_std_pulseq_percent = (rms_all_percent(4) - rms_all_percent(2)) ;

max_temp_drift_all_mean = mean(max_temp_drift_all) ;
max_temp_drift_all_std = std(max_temp_drift_all) ;
max_temp_drift_all_percent = max_temp_drift_all_std./max_temp_drift_all_mean*100 ;

max_temp_drift_all_mean_ice_percent = (max_temp_drift_all_mean(2) - max_temp_drift_all_mean(1)) / max_temp_drift_all_mean(1) * 100 ;
max_temp_drift_all_mean_gt_percent = (max_temp_drift_all_mean(4) - max_temp_drift_all_mean(3)) / max_temp_drift_all_mean(3) * 100 ;
max_temp_drift_all_mean_product_percent = (max_temp_drift_all_mean(3) - max_temp_drift_all_mean(1)) / max_temp_drift_all_mean(1) * 100 ;
max_temp_drift_all_mean_pulseq_percent = (max_temp_drift_all_mean(4) - max_temp_drift_all_mean(2)) / max_temp_drift_all_mean(2) * 100 ;

max_temp_drift_all_std_ice_percent = (max_temp_drift_all_percent(2) - max_temp_drift_all_percent(1)) ;
max_temp_drift_all_std_gt_percent = (max_temp_drift_all_percent(4) - max_temp_drift_all_percent(3)) ;
max_temp_drift_all_std_product_percent = (max_temp_drift_all_percent(3) - max_temp_drift_all_percent(1)) ;
max_temp_drift_all_std_pulseq_percent = (max_temp_drift_all_percent(4) - max_temp_drift_all_percent(2)) ;

save cimax_metrics_epi_mean_std sfnr_mean_all rdc_all rms_all max_temp_drift_all  sfnr_mean_all_mean_ice_percent sfnr_mean_all_mean_gt_percent sfnr_mean_all_mean_product_percent sfnr_mean_all_mean_pulseq_percent...
sfnr_mean_all_std_ice_percent sfnr_mean_all_std_gt_percent sfnr_mean_all_std_product_percent sfnr_mean_all_std_pulseq_percent...
rdc_all_mean_ice_percent rdc_all_mean_gt_percent rdc_all_mean_product_percent rdc_all_mean_pulseq_percent...
rdc_all_std_ice_percent rdc_all_std_gt_percent rdc_all_std_product_percent rdc_all_std_pulseq_percent...
rms_all_mean_ice_percent rms_all_mean_gt_percent rms_all_mean_product_percent rms_all_mean_pulseq_percent...
rms_all_std_ice_percent rms_all_std_gt_percent rms_all_std_product_percent rms_all_std_pulseq_percent...
max_temp_drift_all_mean_ice_percent max_temp_drift_all_mean_gt_percent max_temp_drift_all_mean_product_percent max_temp_drift_all_mean_pulseq_percent...
max_temp_drift_all_std_ice_percent max_temp_drift_all_std_gt_percent max_temp_drift_all_std_product_percent max_temp_drift_all_std_pulseq_percent ;

col = {'SFNR', 'RDC', 'RMS (‰)', 'Drift (%)'} ;
row = {'Product+ICE', 'Pulseq+ICE', 'Product+GT', 'Pulseq+GT'} ;
data = {[num2str(round(sfnr_mean_all_mean(1))) '±' num2str(round(sfnr_mean_all_std(1))) ' (' num2str(round(sfnr_mean_all_percent(1))) '%)'], [num2str(rdc_all_mean(1), '%.2f') '±' num2str(rdc_all_std(1), '%.2f') ' (' num2str(round(rdc_all_percent(1))) '%)'],...
    [num2str(rms_all_mean(1), '%.2f') '±' num2str(rms_all_std(1), '%.2f') ' (' num2str(round(rms_all_percent(1))) '%)'], ...
    [num2str(max_temp_drift_all_mean(1), '%.2f') '±' num2str(max_temp_drift_all_std(1), '%.2f') ' (' num2str(round(max_temp_drift_all_percent(1))) '%)']; ...
    
    [num2str(round(sfnr_mean_all_mean(2))) '±' num2str(round(sfnr_mean_all_std(2))) ' (' num2str(round(sfnr_mean_all_percent(2))) '%)'], [num2str(rdc_all_mean(2), '%.2f') '±' num2str(rdc_all_std(2), '%.2f') ' (' num2str(round(rdc_all_percent(2))) '%)'],...
    [num2str(rms_all_mean(2), '%.2f') '±' num2str(rms_all_std(2), '%.2f') ' (' num2str(round(rms_all_percent(2))) '%)'], ...
    [num2str(max_temp_drift_all_mean(2), '%.2f') '±' num2str(max_temp_drift_all_std(2), '%.2f') ' (' num2str(round(max_temp_drift_all_percent(2))) '%)']; ...
    
    [num2str(round(sfnr_mean_all_mean(3))) '±' num2str(round(sfnr_mean_all_std(3))) ' (' num2str(round(sfnr_mean_all_percent(3))) '%)'], [num2str(rdc_all_mean(3), '%.2f') '±' num2str(rdc_all_std(3), '%.2f') ' (' num2str(round(rdc_all_percent(3))) '%)'],...
    [num2str(rms_all_mean(3), '%.2f') '±' num2str(rms_all_std(3), '%.2f') ' (' num2str(round(rms_all_percent(3))) '%)'], ...
    [num2str(max_temp_drift_all_mean(3), '%.2f') '±' num2str(max_temp_drift_all_std(3), '%.2f') ' (' num2str(round(max_temp_drift_all_percent(3))) '%)']; ...
    
    [num2str(round(sfnr_mean_all_mean(4))) '±' num2str(round(sfnr_mean_all_std(4))) ' (' num2str(round(sfnr_mean_all_percent(4))) '%)'], [num2str(rdc_all_mean(1), '%.2f') '±' num2str(rdc_all_std(4), '%.2f') ' (' num2str(round(rdc_all_percent(4))) '%)'],...
    [num2str(rms_all_mean(4), '%.2f') '±' num2str(rms_all_std(4), '%.2f') ' (' num2str(round(rms_all_percent(4))) '%)'], ...
    [num2str(max_temp_drift_all_mean(4), '%.2f') '±' num2str(max_temp_drift_all_std(4), '%.2f') ' (' num2str(round(max_temp_drift_all_percent(4))) '%)']; } ;
h = figure ;
u = uitable('columnname', col, 'rowname', row,'Position', [20 20 500 70], 'data', data) ;
set(u, 'ColumnWidth', {150,150,150,150,150}) ;
table_extent = get(u,'Extent');
set(u,'Position',[1 1 table_extent(3) table_extent(4)])
figure_size = get(h,'outerposition');
desired_fig_size = [figure_size(1) figure_size(2) table_extent(3) table_extent(4)+80];
set(h,'outerposition', desired_fig_size+20);
% set(gcf, 'Position', get(0, 'Screensize'));
print('-dpng', 'Temporal_quality_metrics', '-r0'); 

% Wilcoxon rank sum test
[p_sfnr_ice, h_sfnr_ice] = ranksum(sfnr_mean_all(:,1), sfnr_mean_all(:,2)) ; % product vs pulseq on ice
[p_sfnr_gt, h_sfnr_gt] = ranksum(sfnr_mean_all(:,3), sfnr_mean_all(:,4)) ; % product vs pulseq on gt
[p_sfnr_product, h_sfnr_product] = ranksum(sfnr_mean_all(:,1), sfnr_mean_all(:,3)) ; % ice vs gt on product
[p_sfnr_pulseq, h_sfnr_pulseq] = ranksum(sfnr_mean_all(:,2), sfnr_mean_all(:,4)) ; % ice vs gt on pulseq

[p_rms_ice, h_rms_ice] = ranksum(rms_all(:,1), rms_all(:,2)) ;
[p_rms_gt, h_rms_gt] = ranksum(rms_all(:,3), rms_all(:,4)) ;
[p_rms_product, h_rms_product] = ranksum(rms_all(:,1), rms_all(:,3)) ;
[p_rms_pulseq, h_rms_pulseq] = ranksum(rms_all(:,2), rms_all(:,4)) ;

[p_rdc_ice, h_rdc_ice] = ranksum(rdc_all(:,1), rdc_all(:,2)) ;
[p_rdc_gt, h_rdc_gt] = ranksum(rdc_all(:,3), rdc_all(:,4)) ;
[p_rdc_product, h_rdc_product] = ranksum(rdc_all(:,1), rdc_all(:,3)) ;
[p_rdc_pulseq, h_rdc_pulseq] = ranksum(rdc_all(:,2), rdc_all(:,4)) ;

[p_max_drift_ice, h_max_drift_ice] = ranksum(max_temp_drift_all(:,1), max_temp_drift_all(:,2)) ;
[p_max_drift_gt, h_max_drift_gt] = ranksum(max_temp_drift_all(:,3), max_temp_drift_all(:,4)) ;
[p_max_drift_product, h_max_drift_product] = ranksum(max_temp_drift_all(:,1), max_temp_drift_all(:,3)) ;
[p_max_drift_pulseq, h_max_drift_pulseq] = ranksum(max_temp_drift_all(:,2), max_temp_drift_all(:,4)) ;

%% plot signal image
lineWidth = 2 ;fontSize = 24 ;
figure ;
t = tiledlayout(1, 4) ;
nexttile
imshow(rot90(Iave_all(:,:,round(Ndata/2),1)),[]) ;
hold on ;
line([Nx Nx],[0 Ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
% axis normal
% axis on
set(ax,'XTick',[], 'YTick', [],'FontSize', fontSize, 'fontWeight', 'bold');
colormap gray ;

nexttile
imshow(rot90(Iave_all(:,:,round(Ndata/2),2)),[]) ;
hold on ;
line([Nx Nx],[0 Ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;

set(ax,'XTick',[], 'YTick', [],'FontSize', fontSize, 'fontWeight', 'bold');
colormap gray ;

nexttile
imshow(rot90(Iave_all(:,:,round(Ndata/2),3)),[]) ;
hold on ;
line([Nx Nx],[0 Ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', [],'FontSize', fontSize, 'fontWeight', 'bold');
colormap gray ;

nexttile
imshow(rot90(Iave_all(:,:,round(Ndata/2),4)),[]) ;
hold on ;
ax = gca;
set(ax,'XTick',[], 'YTick', [],'FontSize', fontSize, 'fontWeight', 'bold');
colormap gray ;
t.Padding = 'none' ;
t.TileSpacing = 'none' ;
set(gcf, 'Position', get(0, 'Screensize'));
exportgraphics(gcf,'signalImage.png') ;

%%
figure ;
t = tiledlayout(1, 4) ;
std_factor = 3 ;
nexttile
Isd_temp = squeeze(Isd_all(:,:,round(Ndata/2),1)) .* mask{1}.ROImask ;
Isd_temp(Isd_temp==0) = NaN ;
Isd_temp_mean = std_factor*mean(Isd_temp(:), 'omitnan') ;
imshow(rot90(Isd_all(:,:,round(Ndata/2),1)),[0 Isd_temp_mean]) ;
hold on ;
line([Nx Nx],[0 Ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
% axis normal
% axis on
set(ax,'XTick',[], 'YTick', [],'FontSize', fontSize, 'fontWeight', 'bold');
colormap gray ;
% colorbar ;
% clim([Iave_min, Iave_max]) ;

nexttile
Isd_temp = squeeze(Isd_all(:,:,round(Ndata/2),2)) .* mask{2}.ROImask ;
Isd_temp(Isd_temp==0) = NaN ;
Isd_temp_mean = std_factor*mean(Isd_temp(:), 'omitnan') ;
imshow(rot90(Isd_all(:,:,round(Ndata/2),2)),[0 Isd_temp_mean]) ;
hold on ;
line([Nx Nx],[0 Ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', [],'FontSize', fontSize, 'fontWeight', 'bold');
colormap gray ;

nexttile
Isd_temp = squeeze(Isd_all(:,:,round(Ndata/2),3)) .* mask{3}.ROImask ;
Isd_temp(Isd_temp==0) = NaN ;
Isd_temp_mean = std_factor*mean(Isd_temp(:), 'omitnan') ;
imshow(rot90(Isd_all(:,:,round(Ndata/2),3)),[0 Isd_temp_mean]) ;
hold on ;
line([Nx Nx],[0 Ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', [],'FontSize', fontSize, 'fontWeight', 'bold');
colormap gray ;

nexttile
Isd_temp = squeeze(Isd_all(:,:,round(Ndata/2),4)) .* mask{4}.ROImask ;
Isd_temp(Isd_temp==0) = NaN ;
Isd_temp_mean = std_factor*mean(Isd_temp(:), 'omitnan') ;
imshow(rot90(Isd_all(:,:,round(Ndata/2),4)),[0 Isd_temp_mean]) ;
hold on ;
ax = gca;
set(ax,'XTick',[], 'YTick', [],'FontSize', fontSize, 'fontWeight', 'bold');
colormap gray ;
t.Padding = 'none' ;
t.TileSpacing = 'none' ;
set(gcf, 'Position', get(0, 'Screensize'));
exportgraphics(gcf,'noise_image.png') ;
%% Weisskoff Graph
r1 = 1 ; r2 = ROIsize ;
rr = (r1:r2) ;
fcalc_all = zeros(ROIsize, Ndata, Nrecon) ;
for i = 1:Ndata
    for j = 1:Nrecon
    fcalc_all(:,i,j) = F_weis_all(1,i,j) ./ rr' ;
    end
end

MarkerSize = 50 ;

figure ;
loglog(rr, squeeze(F_weis_all(:,round(Ndata/2),1)), 'r-', rr, squeeze(fcalc_all(:,round(Ndata/2),1)), 'r--', 'LineWidth', lineWidth) ;
hold on ;
loglog(rr, squeeze(F_weis_all(:,round(Ndata/2),1)), 'r.', 'MarkerSize', MarkerSize) ;

loglog(rr, squeeze(F_weis_all(:,round(Ndata/2),2)), 'b-', rr, squeeze(fcalc_all(:,round(Ndata/2),2)), 'b--', 'LineWidth', lineWidth) ;
loglog(rr, squeeze(F_weis_all(:,round(Ndata/2),2)), 'b.', 'MarkerSize', MarkerSize) ;

loglog(rr, squeeze(F_weis_all(:,round(Ndata/2),3)),'Color', [0 0.5 0], 'LineStyle','-', 'LineWidth', lineWidth) ;
loglog(rr, squeeze(fcalc_all(:,round(Ndata/2),3)),'Color', [0 0.5 0],'LineStyle', '--', 'LineWidth', lineWidth) ;
loglog(rr, squeeze(F_weis_all(:,round(Ndata/2),3)),'Color', [0 0.5 0], 'Marker', '.', 'MarkerSize', MarkerSize) ;

loglog(rr, squeeze(F_weis_all(:,round(Ndata/2),4)), 'k-', rr, squeeze(fcalc_all(:,round(Ndata/2),4)), 'k--', 'LineWidth', lineWidth) ;
loglog(rr, squeeze(F_weis_all(:,round(Ndata/2),4)), 'k.', 'MarkerSize', MarkerSize) ;

grid
xlabel('ROI full width, voxels');
ylabel('Relative std, %');
legend('\color{red} Product+ICE', '','','\color{blue} Pulseq+ICE', '','','\color[rgb]{0,0.5,0} Product+GT','','', '\color{black} Pulseq+GT',  'Location', 'southwest') ;
title(['RDC = ', num2str(squeeze(rdc_all(round(Ndata/2),:)))]);
xlim([r1 r2]);
% ylim([0 1]) ;
set(gca, 'FontSize', 2*fontSize, 'FontWeight', 'bold') ;
set(gcf, 'Position', get(0, 'Screensize'));
exportgraphics(gcf,'Weisskoff_Graph.png') ;


%% plot masks
magFactor = 400 ;
red=zeros(64,64,3);
red(:,:,1)=0.8;

figure;
imshow(rot90(Iave(:,:,end)),[],'InitialMagnification', magFactor);

colormap(gray);
hold all;
h=imshow(red,'InitialMagnification', magFactor);
set(h,'AlphaData',rot90(mask{end}.ROImask)) ;
exportgraphics(gcf,'epi_mask.png') ;