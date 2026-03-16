% Convert .seq files to .pge binary files for execution on GE

% get Pulseq toolbox
system('git clone --branch master git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab

% get toolbox to convert .seq file to a PulSeg sequence (psq) object
system('git clone git@github.com:HarmonizedMRI/pulseg.git');
addpath pulseg/matlab
addpath(genpath('pulseg/matlab/third_party'));

% get toolbox for plotting psq object and exporting to binary file for GE 
system('git clone git@github.com:HarmonizedMRI/pge2.git');
addpath pge2/matlab

% To load the ScanArchive raw data files you will need the Orchestra toolbox
% which is available for download at http://weconnect.gehealthcare.com/ 
% addpath ~/Programs/orchestra-sdk-2.1-1.matlab/

% Options
scanner_pge_location = '/export/home/sdc/PulseqQA/';
opuser1 = [21];     % Determines .entry file number, i.e., pge<opuser1>.entry

% Loop over the .seq files
scans = {'QA_epi_final.seq'};
pislquant = [100];  % number of ADC events to use for receive gain calibration in Auto Prescan

for s = 1 : length(scans)
    seq_name = erase(scans{s}, {'.seq'});

    % Convert .seq file to a PulSeg sequence (psq) object
    psq = pulseg.fromSeq(strcat(seq_name, '.seq'));

    % Check PNS, timing, and b1/gradient limits
    PNSwt = [0.8 1 0.7];   % directional PNS weights, see pge2.pns()
    params = pge2.check(psq, sys_ge, 'PNSwt', PNSwt);

    % Write to .pge file
    pge2.serialize(psq, strcat(seq_name, '.pge'), 'pislquant', pislquant(s), 'params', params);

    % Write the corresponding .entry file.
    pge2.writeentryfile(opuser1(s), seq_name, 'path', scanner_pge_location);

    % (Optional) Validate psq representation against the original .seq file
    %seq = mr.Sequence();
    seq.read(strcat(seq_name, '.seq'));
    pge2.validate(psq, sys_ge, seq, [], 'row', [], 'plot', false);

    % (Optional) Save psq object as .mat file for Matlab runtime based scanner workflow,
    % see https://github.com/HarmonizedMRI/pge2/tree/main/scanner/fov_prescription for details.
    save(seq_name, 'psq', 'params', 'pislquant');  
end
    
