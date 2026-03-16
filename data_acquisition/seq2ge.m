% Convert .seq files to .pge binary files for execution on GE

scanner_pge_location = '/srv/nfs/psd/usr/psd/pulseq/v7/QA/';

% Loop over the .seq files
scans = {'QA_epi_final.seq', 'QA_T1_final.seq'};
pislquant = [100 100];  % number of ADC events to use for receive gain calibration in Auto Prescan
opuser1 = [21 22];     % Determines .entry file number, i.e., pge<opuser1>.entry

for s = 1 %: length(scans)
    fprintf('\nConverting %s.seq to .pge file:\n', scans{s});

    seq_name = erase(scans{s}, {'.seq'});

    % Get system limits
    seq = mr.Sequence();
    seq.read(strcat(seq_name, '.seq'));

    % Convert .seq file to a PulSeg sequence (psq) object
    psq = pulseg.fromSeq(strcat(seq_name, '.seq'));

    % GE system limits
    b1_max   = seq.sys.maxB1/sys.gamma/1e-4;     % Gauss
    g_max    = seq.sys.maxGrad/sys.gamma*100;    % Gauss/cm
    slew_max = seq.sys.maxSlew/sys.gamma/10;     % Gauss/cm/ms
    sys_ge_tmp = pge2.opts(sys_ge.psd_rf_wait, sys_ge.psd_grd_wait, ...
        b1_max, g_max, slew_max, coil);

    % Check PNS, timing, and b1/gradient limits
    PNSwt = [0.8 1 0.7];   % directional PNS weights, see pge2.pns()
    params = pge2.check(psq, sys_ge_tmp, 'PNSwt', PNSwt);

    % Write to .pge file
    pge2.serialize(psq, strcat(seq_name, '.pge'), 'pislquant', pislquant(s), 'params', params);

    % Write the corresponding .entry file.
    pge2.writeentryfile(opuser1(s), seq_name, 'path', scanner_pge_location);

    % (Optional) Validate psq representation against the original .seq file
    pge2.validate(psq, sys_ge_tmp, seq, [], 'row', [], 'plot', false);

    % (Optional) Save psq object as .mat file for Matlab runtime based scanner workflow,
    % see https://github.com/HarmonizedMRI/pge2/tree/main/scanner/fov_prescription for details.
    save(seq_name, 'psq', 'params', 'pislquant');  
end
    
