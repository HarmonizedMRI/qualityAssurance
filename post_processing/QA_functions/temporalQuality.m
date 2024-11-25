%% This is the temporal quality analysis for the quality assurance protocol.
% transversal multi-slice EPI sequence with ramp smapling and navigator:
% write_QA_Tran_EPIrs.m
function [Idiff_noise, Iave, Isd, sfnr, spectrum, F_weis,qa_metrics,mask_all] = EPI_QA_ABCD(vol4D, meta)
%% load images
Nx = size(vol4D, 1) ; Ny = size(vol4D, 2) ; % matrix size
% meta information
FOVx = 220 ; FOVy = 220 ;
meta.TR = 2000 ;
meta.sx = FOVx/Nx ; % pixel spacing
meta.sy = FOVy/Ny ; % pixel spacing
meta.sz = 3 ; % slice thickness
TR = meta.TR/1000 ;

selectedSlice = ceil(size(vol4D,3)/2) ;
% The first two volumes are discarded to allow NMR and eddy current
% equilibrium to be achieved.
data = squeeze(vol4D(:, :, selectedSlice, 3:end) ) ;

%% mask making
% detect edge of the 2D epi image
im_edge = double(edge(data(:,:,selectedSlice)) ) ;
mask_all.im_edge = im_edge ;

[im_X, im_Y] = ind2sub(size(im_edge), find (im_edge)) ;
% fit for the center and radius of the circular phantom
[Xc, Yc, R_circ] = circfit(im_X, im_Y) ;
% phantom_mask = makeCircleMask([Nx, Ny], Xc, Yc, R) ;
Xc = round(Xc) ;
Yc = round(Yc) ;
pxsize = FOVx/Nx * FOVy/Ny ;
ROI_area = 8 * 100 ; % backgorund ROI area = 10 cm2
R_circ = round(R_circ) ;
R1_circ = ceil(R_circ * 0.1) ;
% left background ROI
Wl = size(R1_circ:Xc-R_circ-R1_circ, 2) ;
Hl2 = round(ROI_area / Wl / pxsize / 2) ;
bkg_mask = zeros(Nx,Ny) ;
bkg_mask(R1_circ:Xc-R_circ-R1_circ,Yc-Hl2+1:Yc+Hl2) = 1 ;
% right background ROI
Wr = size(Xc+R_circ+R1_circ-1:Nx-R1_circ+1, 2) ;
Hr2 = round(ROI_area / Wr / pxsize / 2) ;
bkg_mask(Xc+R_circ+R1_circ:Nx-R1_circ+1,Yc-Hr2+1:Yc+Hr2) = 1 ;
mask_all.bkg_mask = bkg_mask ;

ROImask = zeros(Nx, Ny) ;
% 15*15 voxel ROI centered in the phantom for matrix size = 64*64
R = 15 ;
ro2 = fix(R/2);             %half ROI size
X1 = Xc-ro2 ;            %Beginging of masked image in X
X2 = X1 + R - 1;            %End of masked image X
Y1 = Yc-ro2;                    %Beginning of masked image in Y
Y2 = Y1 + R - 1;                    %Beginning of masked image in Y
r1 = 1;                     %First ROI voxel
r2 = R;                     %Last ROI voxel

ROImask(X1:X2, Y1:Y2) = 1 ;
mask_all.ROImask = ROImask ;

%  set up input and ROI mask

npx_roi = sum(ROImask(:));       %Total number of voxels in ROI
i1=1;                           % initial time frame
i2=size(data,3);                % final time frame

N = i2 - i1 + 1;                % num time frames
M = R - 1 + 1;                % num ROI's
weisskoff = zeros(N,M);         % matrix with frames by ROIs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  begin loop through images

Iodd = zeros(Nx^2,1);     %Images in odd positions as vector
Ieven = zeros(Nx^2,1);    %Images in even positions as vector
Sy = zeros(Nx^2,1);
Syy = zeros(Nx^2,1);      % Sum t=1,nFrame; It*It
Syt = zeros(Nx^2,1);      % Sum t=1,nFrame; t*It
St = 0;                     % Sum t=1,nFrame; t
Stt = 0;                    % Sum t=1,nFrame; t^2 
S0 = 0;                     % Time counter 1,2,3... nFrame
img = zeros(Nx) ;
avg_signal_roi = zeros(1, i2) ;

if(mod(N,2)==1) 
    even_tf_flag=0;
else
    even_tf_flag=1;
end

for j = i1:i2                       %For each tFrame
    workingSlice = data(:,:,j);     %Get appropiate slice
    I = workingSlice(:);            %Image slice as vector
    clear workingSlice;
    if(mod(j,2)==1)
        if even_tf_flag
            Iodd = Iodd + I; %Add odd images together
        else
            even_tf_flag=1;
        end
    else
      Ieven = Ieven + I;        %Add even images together
    end
    Sy = Sy + I;             % Add all time frames
    Syt = Syt + I*j;         % Calc Sum t=1,nFrame; t*It
    Syy = Syy + I.*I;        % Calc Sum t=1,nFrame; It*It
    S0 = S0 + 1;             % Update counter
    St = St + j;             % Calc Sum t=1,nFrame; t
    Stt = Stt + j*j;         % Calc Sum t=1,nFrame; t^2 
    img(:) = I;
    sub = img(X1:X2,Y1:Y2);         % masked part of image phantom
    avg_signal_roi(S0) = sum(sum(sub))/npx_roi;    % average signal intensity trough ROI
    for r = r1:r2                   % each roi size
      ro2 = fix(r/2);
      x1 = Xc - ro2;
      x2 = x1 + r - 1;
      y1 = Yc - ro2;
      y2 = y1 + r - 1;
      sub = img(x1:x2,y1:y2);
      weisskoff(j-i1+1, r) = mean(sub(:));   %mean of masked image through time and roi size
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  write out diff image

Isub = Iodd - Ieven;        %Isub in vector form
Isub = reshape(Isub, [Nx, Ny]);
sub = Isub(X1:X2,Y1:Y2);     %Isub within ROI
varI = var(sub(:));         %Variance of Isub within ROI only
Idiff_noise=Isub/(fix(N/2)) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  write out ave image

Iave = Sy/N;                %Iave in vector form
Iave = reshape(Iave, [Nx,Ny]) ;
sub = Iave(X1:X2,Y1:Y2);     %Iave with ROI masked
meanI = mean(sub(:));       %Mean of Iave within ROI

% second-order polynominal detrending
im_detrended = zeros(size(data)) ;
Isd = zeros(Nx, Ny) ;
x = (1:198)' ;
for i = 1:Nx
    for j = 1:Ny
        p = polyfit(x, squeeze(data(i,j,:)), 2) ;
        predicted = polyval(p, x) ;
        im_detrended(i,j,:) = squeeze(data(i,j,:)) - predicted ;
        Isd(i,j) = std(im_detrended(i,j,:)) ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make sfnr image
sfnr = Iave./Isd;     %sfnr image in vector
sfnr(sfnr>10000) = NaN ;
sfnr = reshape(sfnr, [Nx,Ny]) ;
sub = sfnr(X1:X2,Y1:Y2);       %sfnr image within ROI
sfnrI = mean(sub(:), 'omitnan');         %mean sfnr value within ROI


snr = meanI/sqrt(varI/N); %SNR value


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%  Do fluctation analysis

% First plot: Percent fluctuation and percent drift

x=(1:N);
p=polyfit(x,avg_signal_roi,2);
yfit = polyval(p, x);
y = avg_signal_roi - yfit;
time_minutes = x*TR/60;

m=mean(avg_signal_roi);
fluct_sd=std(y);
drift = ((yfit(N)-yfit(1))/m);
drift_per_minute = ((yfit(N)-yfit(1))/m)/time_minutes(end); %Drift using trend line. Last intensity-initial divided by mean of trend line (averaged by minute)
maxabsdrift = (max(avg_signal_roi) - min(avg_signal_roi))/m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Second plot frequency analysis of fluctuations (detrended)

z = fft(y);

% fs = 1/TR;
nf = N/2+1;
% freq_axis = 0.5*(1:nf)*fs/nf;
spectrum = z(1:fix(nf)) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  now do analysis for each roi size
F_weis = zeros(1, R) ;
t = (1:N);
for r = r1:r2
  y = weisskoff(:, r)';
  yfit = polyval(polyfit(t, y, 2), t);  % 2nd order trend
  F_weis(r) = std(y - yfit)/mean(yfit); % coefficient of variation (CV, the SD of a time-series divided by the mean of the time-series)
end
% rr = (r1:r2);
F_weis = 100*F_weis;              % percent
% fcalc = F_weis(1)./rr;
rdc = F_weis(1)/F_weis(r2);	% decorrelation distance, radius of decorrelation (RDC) = CV(1)/CV(Nmax)

qa_metrics=struct('meanI',meanI,'fluct_sd',fluct_sd,'snr',snr,'sfnr_mean',sfnrI,'rms',100*fluct_sd/m,'temp_drift',100*drift, 'temp_drift_per_minute', 100*drift_per_minute,...
          'max_temp_drift',100*maxabsdrift,'rdc',rdc);
end