%% load Siemens .dat data
datafile = dir('meas_MID00245_FID24189_pulseq_QA_fmri_sli27_rep200_run4.dat') ;
twix_obj = mapVBVD(datafile.name) ; % load Siemens .dat data
phasecorData = twix_obj{end}.phasecor.unsorted() ; % navigator data
imageData = twix_obj{end}.image.unsorted() ;
Nsample = size(imageData, 1) ; % number of ADC
Ncoil = size(imageData, 2) ; % number of coils
hdr = twix_obj{end}.hdr ;
Npe = twix_obj{end}.image.NLin ;
Nslice = twix_obj{end}.image.NSli ;
Nrep = twix_obj{1, 2}.image.NRep ;
Nnav = 3 ;
REP_image = reshape(twix_obj{end}.image.Rep, [Npe, Nslice, Nrep]) ;
REP_phasecor = reshape(twix_obj{end}.phasecor.Rep, [Nnav, Nslice, Nrep]) ;
SLC_image = reshape(twix_obj{end}.image.Sli, [Npe, Nslice, Nrep]) ;
SLC_phasecor = reshape(twix_obj{end}.phasecor.Sli, [Nnav, Nslice, Nrep]) ;
LIN_image = reshape(twix_obj{end}.image.Lin, [Npe, Nslice, Nrep]) ;
LIN_phasecor = reshape(twix_obj{end}.phasecor.Lin, [Nnav, Nslice, Nrep]) ;
NAV_image = reshape(zeros(1,size(imageData,3)), [Npe, Nslice, Nrep]) ;
NAV_phasecor = reshape(ones(1,size(phasecorData,Nnav)), [Nnav, Nslice, Nrep]) ;
AVG_image = reshape(twix_obj{end}.image.Ave, [Npe, Nslice, Nrep]) ;
AVG_phasecor = reshape(twix_obj{end}.phasecor.Ave, [Nnav, Nslice, Nrep]) ;
SEG_image = reshape(twix_obj{end}.image.Seg, [Npe, Nslice, Nrep]) ;
SEG_phasecor = reshape(twix_obj{end}.phasecor.Seg, [Nnav, Nslice, Nrep]) ;
REV_image = reshape(twix_obj{end}.image.IsReflected, [Npe, Nslice, Nrep]) ;
REV_phasecor = reshape(twix_obj{end}.phasecor.IsReflected, [Nnav, Nslice, Nrep]) ;

REPlbl = [REP_phasecor; REP_image] ; REPlbl = REPlbl(:)-1 ;
SLClbl = [SLC_phasecor; SLC_image] ; SLClbl = SLClbl(:)-1 ;
LINlbl = [LIN_phasecor; LIN_image] ; LINlbl = LINlbl(:)-1 ;
AVGlbl = [AVG_phasecor; AVG_image] ; AVGlbl = AVGlbl(:)-1 ;
SEGlbl = [SEG_phasecor; SEG_image] ; SEGlbl = SEGlbl(:)-1 ;
REVlbl = [REV_phasecor; REV_image] ; REVlbl = REVlbl(:) ;
NAVlbl = [NAV_phasecor; NAV_image] ; NAVlbl = NAVlbl(:) ;
LINlbl_min = min(LINlbl) ; LINlbl_max = max(LINlbl) ; % min and max LIN
SLClbl_min = min(SLClbl) ; SLClbl_max = max(SLClbl) ;
AVGlbl_min = min(AVGlbl) ; AVGlbl_max = max(AVGlbl) ;
SEGlbl_min = min(SEGlbl) ; SEGlbl_max = max(SEGlbl) ;
REPlbl_min = min(REPlbl) ; REPlbl_max = max(REPlbl) ;
Nscan = size(LINlbl, 1) ; % total scan number
% other sequence parameters
% encoded space fov
fov = 1e3 * [220, 220, 144] *1e-3 ; %seq.getDefinition('FOV') ;
% x: readout; y: phase encoding; z: partition/slice

e_fov_x = fov(1) ; e_fov_y = fov(2) ; e_fov_z = 3;%fov(3) ;
readout_os = 2;%seq.getDefinition('ReadoutOversamplingFactor') ;
sliceThickness = 5e-3 ;%seq.getDefinition('SliceThickness') ;

% encoded space matrix size
e_matrixSize_x = Nsample/readout_os ;
e_matrixSize_y = Npe ;
e_matrixSize_z = 1;%SLClbl_max + 1 ;
% reconspace fov and matrix size
r_fov_x = e_fov_x ; r_fov_y = e_fov_y ; r_fov_z = e_fov_z ;
r_matrixSize_x = e_matrixSize_x ;
r_matrixSize_y = e_matrixSize_y ;
r_matrixSize_z = e_matrixSize_z ;

%% combine noise data, ref data, and image data together (typical for Siemens data)
totalData = zeros(Nsample, Ncoil, Nscan) ; % initiliaze combined data
phasecor_count = 1 ; % initialize navigator counter
image_count = 1 ; % initialize image scan counter
for i = 1:Nscan
    if NAVlbl(i) == 1 % navigator scan
        totalData(:,:,i) = phasecorData(:,:,phasecor_count) ;
        phasecor_count = phasecor_count + 1 ;
    else % image scan
        totalData(:,:,i) = imageData(:,:,image_count) ;
        image_count = image_count + 1 ;
    end
end

%%
% Output file Name
filename = 'pulseq_epi_data.h5';
dset = ismrmrd.Dataset(filename);
acqblock = ismrmrd.Acquisition(Nscan) ;
% Set the header elements that don't change
acqblock.head.version(:) = 1 ;
acqblock.head.number_of_samples(:) = Nsample ;
acqblock.head.center_sample(:) = floor(Nsample/2) ;
acqblock.head.active_channels(:) = Ncoil ;
acqblock.head.read_dir  = repmat([0 0 1]', [1 Nscan]) ;
acqblock.head.phase_dir = repmat([0 1 0]', [1 Nscan]) ;
acqblock.head.slice_dir = repmat([1 0 0]', [1 Nscan]) ;
%%
% Loop over the acquisitions, set the header, set the data and append
for i = 1:Nscan
    % Set the header elements that change from acquisition to the next
    % c-style counting
    acqblock.head.scan_counter(i) = i-1 ;
    % Note next entry is k-space encoded line number (not i which
    % is just the sequential acquisition number)
    acqblock.head.idx.kspace_encode_step_1(i) = LINlbl(i) ;
    acqblock.head.idx.slice(i) = SLClbl(i) ;
    acqblock.head.idx.repetition(i) = REPlbl(i) ;
    acqblock.head.idx.average(i) = AVGlbl(i) ;
    acqblock.head.idx.segment(i) = SEGlbl(i) ;

    % Set the flags
    acqblock.head.flagClearAll(i) ;
    if NAVlbl(i) == 1
        acqblock.head.flagSet('ACQ_IS_PHASECORR_DATA', i) ;
%         acqblock.head.flagSet('ACQ_IS_RTFEEDBACK_DATA', i) ;
    end
    if REVlbl(i) == 1
        acqblock.head.flagSet('ACQ_IS_REVERSE', i) ;
        % flip and fill the data
        acqblock.data{i} = flip(squeeze(totalData(:,:,i))) ;
    else
        % fill the data
        acqblock.data{i} = squeeze(totalData(:,:,i)) ;
    end
end

% Append the acquisition block
dset.appendAcquisition(acqblock) ;

%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill the xml header %
%%%%%%%%%%%%%%%%%%%%%%%%
% We create a matlab struct and then serialize it to xml.
% Look at the xml schema to see what the field names should be

header = [];

% Experimental Conditions (Required)
header.experimentalConditions.H1resonanceFrequency_Hz = 123194105 ; % 3T

% Acquisition System Information (Optional)
header.acquisitionSystemInformation.systemVendor = 'ISMRMRD Labs' ;
header.acquisitionSystemInformation.systemModel = 'Virtual Scanner' ;
header.acquisitionSystemInformation.receiverChannels = Ncoil ;

% The Encoding (Required)
header.encoding.trajectory = 'cartesian' ;
header.encoding.encodedSpace.fieldOfView_mm.x = e_fov_x ;
header.encoding.encodedSpace.fieldOfView_mm.y = e_fov_y ;
header.encoding.encodedSpace.fieldOfView_mm.z = sliceThickness * 1e3 ;
header.encoding.encodedSpace.matrixSize.x = e_matrixSize_x ;
header.encoding.encodedSpace.matrixSize.y = e_matrixSize_y ;
header.encoding.encodedSpace.matrixSize.z = e_matrixSize_z ;
% Recon Space
header.encoding.reconSpace.fieldOfView_mm.x = r_fov_x ;
header.encoding.reconSpace.fieldOfView_mm.y = r_fov_y ;
header.encoding.reconSpace.fieldOfView_mm.z = sliceThickness * 1e3 ;
header.encoding.reconSpace.matrixSize.x = r_matrixSize_x ;
header.encoding.reconSpace.matrixSize.y = r_matrixSize_x ;
header.encoding.reconSpace.matrixSize.z = r_matrixSize_z ;
% Encoding Limits
header.encoding.encodingLimits.kspace_encoding_step_0.minimum = 0 ;
header.encoding.encodingLimits.kspace_encoding_step_0.maximum = Nsample-1 ;
header.encoding.encodingLimits.kspace_encoding_step_0.center = floor(Nsample/2) ;
header.encoding.encodingLimits.kspace_encoding_step_1.minimum = LINlbl_min ;
header.encoding.encodingLimits.kspace_encoding_step_1.maximum = LINlbl_max ;
header.encoding.encodingLimits.kspace_encoding_step_1.center = floor((LINlbl_max+1)/2) ;
header.encoding.encodingLimits.kspace_encoding_step_2.minimum = 0 ;
header.encoding.encodingLimits.kspace_encoding_step_2.maximum = 0 ;
header.encoding.encodingLimits.kspace_encoding_step_2.center = 0 ;
header.encoding.encodingLimits.average.minimum = 0 ;
header.encoding.encodingLimits.average.maximum = 0 ;
header.encoding.encodingLimits.average.center = 0 ;
header.encoding.encodingLimits.slice.minimum = SLClbl_min ;
header.encoding.encodingLimits.slice.maximum = SLClbl_max ;
header.encoding.encodingLimits.slice.center = 0 ;
header.encoding.encodingLimits.repetition.minimum = REPlbl_min ;
header.encoding.encodingLimits.repetition.maximum = REPlbl_max ;
header.encoding.encodingLimits.repetition.center = 0 ;

rampUpTime = hdr.Meas.alRegridRampupTime(1) ;
flatTopTime = hdr.Meas.alRegridFlattopTime(1) ;
rampDownTime = hdr.Meas.alRegridRampdownTime(1) ;
acqDelayTime = hdr.Meas.alRegridDelaySamplesTime(1) ;
numSamples = hdr.Meas.alRegridDestSamples(1) ;
numberOfNavigator = 3 ;
dwellTime = 1e-3 * hdr.Meas.alDwellTime(1);
etl = numSamples ;
header.encoding.trajectoryDescription.identifier = 'ConventionalEPI' ;
header.encoding.trajectoryDescription.userParameterLong   = ...
    struct('value', {etl,numberOfNavigator,rampUpTime,rampDownTime,flatTopTime,acqDelayTime,numSamples},...
    'name',{'etl','numberOfNavigators','rampUpTime','rampDownTime','flatTopTime','acqDelayTime','numSamples'}) ;
header.encoding.trajectoryDescription.userParameterDouble.value   = dwellTime ;
header.encoding.trajectoryDescription.userParameterDouble.name   = 'dwellTime' ;
header.encoding.trajectoryDescription.comment = 'Conventional EPI sequence' ;
header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_0 = 1 ;
header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1 = 1 ;
header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_2 = 1 ;
header.encoding.parallelImaging.calibrationMode = 'other' ;

%% Serialize and write to the data set
xmlstring = ismrmrd.xml.serialize(header) ;
dset.writexml(xmlstring) ;

%% Write the dataset
dset.close() ;