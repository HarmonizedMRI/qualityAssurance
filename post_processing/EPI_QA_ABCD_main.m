clear all; close all ; clc ;
magFactor = 100 ;
lineWidth = 1.5 ;
fontSize = 24 ;
path = pwd ;
TR = 2 ;
%%
imfile = dir('*.nii') ;
Ndata = size(imfile, 1) ;
im = niftiread(imfile(1).name) ;
Nx = size(im, 1) ; Ny = size(im, 2) ; Nrun = size(im, 4) - 2 ; ROIsize = 15 ;
Idiff_noise = zeros(Nx, Ny, Ndata) ;
Iave = zeros(Nx, Ny, Ndata) ;
Isd = zeros(Nx, Ny, Ndata) ;
sfnr = zeros(Nx, Ny, Ndata) ;
nf = fix(Nrun/2+1) ;
spectrum = zeros(nf, Ndata) ;
F_weis = zeros(ROIsize, Ndata) ;
meanI = zeros(1, Ndata) ;
fluct_sd = zeros(1, Ndata) ;
snr = zeros(1, Ndata) ;
sfnr_mean = zeros(1, Ndata) ;
rms = zeros(1, Ndata) ;
temp_drift = zeros(1, Ndata) ;
temp_drift_per_minute = zeros(1, Ndata) ;
max_temp_drift = zeros(1, Ndata) ;
rdc = zeros(1, Ndata) ;

for i = 1:Ndata
    disp(['start processing the ', num2str(i), 'th dataset...']) ;
    vol4D = double(niftiread(imfile(i).name)) ; % [Nx, Ny, Nslice, Nrun]
    if max(vol4D(:))<1 % rescale it.
        vol4D = vol4D / max(vol4D(:)) * 4096 ;
    end
    [Idiff_noise(:,:,i), Iave(:,:,i), Isd(:,:,i), sfnr(:,:,i), spectrum(:,i), F_weis(:,i), qa_metrics, mask_all] = EPI_QA_ABCD(vol4D) ;
    meanI(:,i) = qa_metrics.meanI ;
    fluct_sd(:,i) = qa_metrics.fluct_sd ;
    snr(:,i) = qa_metrics.snr ;
    sfnr_mean(:,i) = qa_metrics.sfnr_mean ;
    rms(:,i) = qa_metrics.rms ;
    temp_drift(:,i) = qa_metrics.temp_drift ;
    temp_drift_per_minute(:,i) = qa_metrics.temp_drift_per_minute ;
    max_temp_drift(:,i) = qa_metrics.max_temp_drift ;
    rdc(:,i) = qa_metrics.rdc ;
    clear qa_metrics ;
end

figure ;
montage(mat2gray(Idiff_noise)); % error image between even and odd frames

figure ;
montage(mat2gray(Iave)); % average iamge across all frames

figure ;
montage(mat2gray(Isd)); % sfnr = Iave./Isd ;

figure ;
montage(mat2gray(sfnr)) ;

% plot frequency analysis of fluctuations (detrended)
figure ;
fs = 1/TR;
f = 0.5*(1:nf)*fs/nf;
plot(f, abs(spectrum(:,:)));
grid
ylabel('spectrum');
xlabel('frequency, Hz');
r1 = 1 ; r2 = ROIsize ;
rr = (r1:r2) ;
figure ;
for i = 1:Ndata
loglog(rr, F_weis(:,i), '-x', rr, F_weis(1,i)./rr', '--', 'LineWidth', lineWidth) ;
hold on ;
end
grid
title(['RDC = ', num2str(rdc)]);
xlabel('ROI full width, voxels');
ylabel('Relative std, %');
xlim([r1 r2]);
legend('measured', 'theoretical') ;
set(gca, 'FontSize', fontSize, 'FontWeight', 'bold') ;
set(gcf, 'Position', get(0, 'Screensize'));
exportgraphics(gcf,'Weisskoff_Graph.png') ;

save results Idiff_noise Iave Isd sfnr spectrum F_weis meanI fluct_sd snr sfnr_mean rms temp_drift temp_drift_per_minute max_temp_drift rdc mask_all ;


