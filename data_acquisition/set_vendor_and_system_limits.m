
vendor = 'g' ;

switch lower(vendor(1))
    case 's'
        sys = mr.opts('maxGrad', max_grad, 'gradUnit','mT/m', ...
              'maxSlew', max_slew, 'slewUnit', 'T/m/s', ...
              'rfDeadTime', 100e-6, ... 
              'rfRingdownTime', 20e-6, ...
              'adcDeadTime', 20e-6, ... 
              'flag_trid', false, ...
              'B0', 2.89);                  % this is Siemens' 3T

    case 'g'
        % System limits used in design.
        % On GE, block boundaries disappear inside segments, so it may be ok
        % to set dead/ringdown times to 0 in practice here.
        sys = mr.opts('maxGrad', max_grad, 'gradUnit','mT/m', ...
              'maxSlew', max_slew, 'slewUnit', 'T/m/s', ...
              'rfDeadTime', 100e-6, ...     % or 0
              'rfRingdownTime', 0e-6, ...  % or 0
              'adcDeadTime', 0e-6, ...     % or 0
              'adcRasterTime', 2e-6, ...    % GE dwell time must be a multiple of 2us
              'rfRasterTime', 4e-6, ...     % 2e-6, or any integer multiple thereof
              'gradRasterTime', 4e-6, ...   % 4e-6, or any integer multiple thereof
              'blockDurationRaster', 4e-6, ... % 4e-6, or any integer multiple thereof
              'flag_trid', true, ...
              'B0', 3.0);

        % additional system limits for GE
        psd_rf_wait  = 50e-6;   % RF–gradient delay (s), scanner-specific
        psd_grd_wait = 50e-6;   % ADC–gradient delay (s), scanner-specific
        b1_max   = sys.maxB1/sys.gamma/1e-4;  % Gauss
        g_max    = max_grad/10;           % Gauss/cm
        slew_max = max_slew/10;           % Gauss/cm/ms
        coil     = 'xrm';        % See pge2.opts(). 'xrm' (MR750), 'hrmw' (Premier), 'magnus', ...
        sys_ge = pge2.opts(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, coil);

    otherwise
        error("Vendor must be 'GE' or 'Siemens' (for now)");
end
