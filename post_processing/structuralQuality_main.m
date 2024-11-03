%% This is structural quality analysis for quality assurance.
% transversal T1-weighted multi-slice spin-echo sequence:
% write_QA_Tran_T1.m
clear all ; close all; clc ;
magFactor = 100 ;
lineWidth = 1.5 ;
fontSize = 24 ;
imfile = dir('*.nii') ;
im = double(niftiread(imfile(1).name)) ;
Ndata = size(imfile, 1) ;
nx = size(im, 1) ; ny = size(im, 2) ; nslice = 11 ; nacq = 2 ;

im_all = zeros(nx,ny,nslice,nacq,Ndata) ;
PSG = zeros(1, Ndata) ;
PIU = zeros(1, Ndata) ;
SNR_1acq = zeros(1, Ndata) ;
SNR_2acq = zeros(1, Ndata) ;

for i = 1:Ndata
    im = double(niftiread(imfile(i).name)) ;
    im_all(:,:,:,:,i) = reshape(im, [nx, ny, nslice, nacq]) ;
    [PSG(i), PIU(i), SNR_1acq(i), SNR_2acq(i), mask_all] = structuralQuality(im_all(:,:,:,:,i)) ;
end
save results im_all PSG PIU SNR_1acq SNR_2acq mask_all ;

%%
figure ;
montage(mat2gray(im_all(:,:,:,1,end))) ;
title('11 slices of the first acquisition') ;
colormap default ;
figure ;
montage(mat2gray(im_all(:,:,:,1,end)-im_all(:,:,:,2,end))) ;
title('The difference between the first and second acquisitions') ;
colormap default ;

figure ;
t = tiledlayout(1, 3) ;
nexttile
im_edge_fuse = labeloverlay(mat2gray(squeeze(im_all(:,:,6,1,end))), mask_all.im_edge,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_edge_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Edge mask','FontSize', fontSize) ;
nexttile
im_phantom_mask_fuse = labeloverlay(mat2gray(squeeze(im_all(:,:,6,1,end))), mask_all.phantom_mask,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_phantom_mask_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Phantom mask','FontSize', fontSize) ;

nexttile
im_largeROI_mask_fuse = labeloverlay(mat2gray(squeeze(im_all(:,:,6,1,end))), mask_all.largeROI_mask,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_largeROI_mask_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Large ROI mask','FontSize', fontSize) ;
t.Padding = 'none';
t.TileSpacing = 'none';


figure ;
t = tiledlayout(1, 4) ;
nexttile
im_bkg_left_fuse = labeloverlay(mat2gray(squeeze(im_all(:,:,6,1,end))), mask_all.bkg_left,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_bkg_left_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Left bkg ROI','FontSize', fontSize) ;
nexttile
im_bkg_right_fuse = labeloverlay(mat2gray(squeeze(im_all(:,:,6,1,end))), mask_all.bkg_right,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_bkg_right_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Right bkg ROI','FontSize', fontSize) ;

nexttile
im_bkg_btm_fuse = labeloverlay(mat2gray(squeeze(im_all(:,:,6,1,end))), mask_all.bkg_btm,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_bkg_btm_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Bottm bkg ROI','FontSize', fontSize) ;
nexttile
im_bkg_top_fuse = labeloverlay(mat2gray(squeeze(im_all(:,:,6,1,end))), mask_all.bkg_top,'Colormap','autumn','Transparency', 0.25) ;
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

figure ;
t = tiledlayout(1, 2) ;
nexttile
mask_temp = mask_all.largeROI_mask_edge+mask_all.smallROI_max_mask ;
mask_temp(mask_temp>0) = 1 ;
im_smallROI_max_fuse = labeloverlay(mat2gray(squeeze(im_all(:,:,6,1,end))), mask_temp,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_smallROI_max_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('small ROI with maximum mean intensity','FontSize', fontSize) ;
nexttile
mask_temp = mask_all.largeROI_mask_edge + mask_all.smallROI_min_mask ;
mask_temp(mask_temp>0) = 1 ;
im_smallROI_min_fuse = labeloverlay(mat2gray(squeeze(im_all(:,:,6,1,end))), mask_temp,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_smallROI_min_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([nx nx],[0 ny+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('small ROI with minimum mean intensity','FontSize', fontSize) ;
t.Padding = 'none';
t.TileSpacing = 'none';

%%
allmask = 2*mask_all.bkg_left+2*mask_all.bkg_right+2*mask_all.bkg_top+2*mask_all.bkg_btm + edge(mask_all.largeROI_mask) + mask_all.im_edge ;
figure ;
im_allMask_fuse = labeloverlay(mat2gray(squeeze(im_all(:,:,6,1,end))), allmask,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_allMask_fuse),[],'InitialMagnification', magFactor) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('allMask','FontSize', fontSize) ;

%%
smallROImask = 1*edge(mask_all.largeROI_mask) + 1*mask_all.smallROI_min_mask + 2*mask_all.smallROI_max_mask ;
figure ;
im_allMask_fuse = labeloverlay(mat2gray(squeeze(im_all(:,:,6,1,end))), smallROImask,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_allMask_fuse),[],'InitialMagnification', magFactor) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('allMask','FontSize', fontSize) ;

%%
allmask = 2*mask_all.bkg_left+2*mask_all.bkg_right+2*mask_all.bkg_top+2*mask_all.bkg_btm + edge(mask_all.largeROI_mask) + mask_all.im_edge ;
smallROImask = 1*edge(mask_all.largeROI_mask) + 1*mask_all.smallROI_min_mask + 2*mask_all.smallROI_max_mask ;
figure ;
im_allMask_fuse = labeloverlay(mat2gray(squeeze(im_all(:,:,6,1,end))), smallROImask+allmask,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_allMask_fuse),[],'InitialMagnification', magFactor) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('allMask','FontSize', fontSize) ;