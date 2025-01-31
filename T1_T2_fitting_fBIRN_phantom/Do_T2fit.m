%clean matlab for a  fresh run of this program
clear all; close all; clc ;
magFactor = 100 ;
lineWidth = 1.5 ;
fontSize = 24 ;

%% T1 fitting
% Where to find the DICOM data 
currentPath = pwd ;
T2path = [currentPath '\T2'] ;
files = dir(T2path) ;
n = 0 ;
% Get list of all filenames in folder
for i = 3:length(files)
    if files(i).isdir
%         info = dicominfo( fullfile(T2path, files(i).name ) ) ;   % Load Dicom tag parameters to 'info'
%         if(strcmp(info.SeriesDescription,'localizer') == 0 && strcmp(info.SeriesDescription,'PhoenixZIPReport') == 0)
            n = n + 1 ;
            dicomFiles{n} = [files(i).name, '\image0001.dcm'] ;
%         end
    end
end
% load the first image for data acquisition information and masking
info_first = dicominfo( fullfile(T2path, dicomFiles{1} ) ) ;
Nread = info_first.Height ;
Nphase = info_first.Width ;
image_first = dicomread(info_first) ;
im_edge = double(edge(image_first) ) ; % detect circular edge
[im_X, im_Y] = ind2sub(size(im_edge), find(im_edge) ) ;
% fit for circular center (Xc, Yc), and radius: R
[Xc, Yc, R] = circfit(im_X, im_Y) ;
phantom_mask = makeCircleMask([Nread Nphase], Xc, Yc, R) ; % phantom mask
% large ROI mask, whose radius is 0.8 * R
largeROI_mask = makeCircleMask([size(phantom_mask,1) size(phantom_mask,2)], Xc, Yc, R*0.8) ; % large ROI mask
largeROI_size = nnz(largeROI_mask) ;

figure ;
t = tiledlayout(1, 3) ;

nexttile
im_edge_fuse = labeloverlay(mat2gray(image_first), im_edge,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_edge_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([Nread Nread],[0 Nphase+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Edge mask','FontSize', fontSize) ;

nexttile
im_phantom_mask_fuse = labeloverlay(mat2gray(image_first), phantom_mask,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_phantom_mask_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
line([Nread Nread],[0 Nphase+1], 'color','y', 'lineWidth',lineWidth) ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Phantom mask','FontSize', fontSize) ;

nexttile
im_largeROI_mask_fuse = labeloverlay(mat2gray(image_first), largeROI_mask,'Colormap','autumn','Transparency', 0.25) ;
imshow(rot90(im_largeROI_mask_fuse),[],'InitialMagnification', magFactor) ;
hold on ;
ax = gca;
set(ax,'XTick',[], 'YTick', []);
colormap default ;
title('Large ROI mask','FontSize', fontSize) ;
t.Padding = 'none';
t.TileSpacing = 'none';
set(gcf, 'Position', get(0, 'Screensize'));

% Browse through all files in folder
imageNum = length(dicomFiles) ;
TEs = zeros(1, imageNum) ;
pixel_value = zeros(largeROI_size, imageNum) ;
image = zeros(Nread, Nphase, imageNum) ;
for i = 1:imageNum
    info = dicominfo( fullfile(T2path, dicomFiles{i} ) ) ;   % Load Dicom tag parameters to 'info'
    TEs(i) = info.PerFrameFunctionalGroupsSequence.Item_1.MREchoSequence.Item_1.EffectiveEchoTime ; % Inversion time TI, this might not be the right parameter
    image(:,:,i) = dicomread(info) ;
    image_temp = image(:,:,i) ;
    pixel_value(:, i) = image_temp(largeROI_mask==1) ;
end
figure ;
montage(mat2gray(image)) ;
title('All images') ;
set(gca, 'FontSize', fontSize, 'FontWeight', 'bold') ;
t = TEs' ;   % from row to column vectors, ' means transpose
y = pixel_value' ;

% Negate the values before zero-crossing (if inversion recovery is used)
% [minval, zeropoint] = min(y) ;
% y_cor = [-y(1:zeropoint); y(zeropoint+1:size(y))] ;

% Do function fitting to the data stored in y_cor
% Use e.g the non-linear fitting function 'lsqcurvefit'
% If for example you want to fit the data to the function:  y = At+Bt^2)
% this is done by first setting initial values for the A and B parameters,
Mxy_init = 1400 ;
T2_init = 80 ;
x0 = [Mxy_init T2_init] ; 
% and then define the model function:
% Mxy(t) = Mxy(t=0) * exp(-t/T2)
F = @(x,t)(x(1)*exp(-t./x(2)) ) ; % Notice: this should be changed to the right function

% Then do the non-linear fitting for every pixel within the large ROI
T2_fit = zeros(largeROI_size, 1) ;
options = optimset('Display','off');
for i = 1:largeROI_size
    [x, resnorm, ~, exitflag, output] = lsqcurvefit(F, x0, t, y(:,i),[],[],options) ;
    T2_fit(i) = x(2) ;
end

%%
markerSize = 36 ;
fontSize = 28 ;
lineWidth = 4 ;
figure ;
hold on
plot(t, y(:,end), 'r.','lineWidth',lineWidth, 'markerSize', markerSize) ; % acquired data
% and plot the model function together with the data
plot(t, F(x, t), 'blue','lineWidth',lineWidth) ; % fitted curve
hold off ;
grid on ;
title('T2 decay of a pixel') ;
xlabel('TE (ms)') ;
ylabel('Pixel intensity') ;
legend('measured data', 'fitted T2 curve') ;
set(gca, 'FontSize', fontSize, 'FontWeight', 'bold') ;
set(gcf, 'Position', get(0, 'Screensize'));
exportgraphics(gcf,'T2decay_example.png') ;

% Print out resulting T2 value
disp(['The fitted T2 value is ' num2str(mean(T2_fit), '%.1f') '±' num2str(std(T2_fit), '%.1f') ' ms'] ) ;

%%
function [circleMask] = makeCircleMask(im_size, Xc, Yc, Radius)
circleMask = zeros(im_size) ;
I = im_size(1) ; J = im_size(2) ;
for i = 1:I
    for j = 1:J
        distance = (double(i) - Xc).^2 + (double(j) - Yc).^2 ;
        if distance <= Radius.^2
            circleMask(i, j) = 1 ;
        end
    end
end
end

function [xc,yc,R,a] = circfit(x,y)
%
%   [xc yx R] = circfit(x,y)
%
%   fits a circle  in x,y plane in a more accurate
%   (less prone to ill condition )
%  procedure than circfit2 but using more memory
%  x,y are column vector where (x(i),y(i)) is a measured point
%
%  result is center point (yc,xc) and radius R
%  an optional output is the vector of coeficient a
% describing the circle's equation
%
%   x^2+y^2+a(1)*x+a(2)*y+a(3)=0
%
%  By:  Izhak bucher 25/oct /1991, 
    x=x(:); y=y(:);
   a=[x y ones(size(x))]\[-(x.^2+y.^2)];
   xc = -.5*a(1);
   yc = -.5*a(2);
   R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
end