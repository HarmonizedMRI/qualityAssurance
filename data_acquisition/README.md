# This is the instruction for data acquisition.
* `write_QA_Tran_EPIrs.m`: to generate `QA_epi.seq` file for EPI scans.
* `write_QA_Tran_T1.m`: to generate QA_T1.seq file for spin-echo scans.
* `20241122_QA_protocol_instruction_siemens.docx`: standard operating procedure for QA measurements.
* `QA_record.xlsx`: Excel sheet for the record of QA measurements.
* `QA_MAGMA.pdf`: product and Pulseq-based sequence protocols on Cima.X for quality assurance.

## Setup

### Set up Python environment 

Need for `mr.makeSLRpulse` call.

On Linux command line:
```bash
sudo apt install python3.12-venv   # or whichever python version you have available
python3 -m venv myvenv
source myvenv/bin/activate
pip install matplotlib
pip install scipy
pip install sigpy
```

Then start MATLAB from within that environment:
```bash
$ source myvenv/bin/active
$ matlab
```

### Set MATLAB paths

From inside MATLAB:
```matlab
>> setup
```


## Usage

1. Edit `set_vendor_and_system_limits.m`:
   1. Set `vendor` to `'Siemens'` or `'GE'`
   3. Set any other scanner-specific system settings if needed.

2. Set paths and execute the `write*.m` scripts:
   ```matlab
   >> write_QA_Tran_EPIrs
   >> write_QA_Tran_T1
   ```

This will create the Pulseq (`.seq`) files needed to execute this QA protocol on the scanner.


### Additional steps for GE users

1. In `set_vendor_and_system_limits.m`, set `'coil'` to the appropriate for your scanner (used to check PNS), e.g., 
   ```matlab
   coil = 'xrm';   % MR750
   ```
1. In `seq2ge.m`, set the desired options:
   1. Set `scanner_pge_location` to the directory on the scanner where you
   wish to put the `.pge` files.
   2. Set `opuser1` to control the desired `.entry` file names. 
   3. Run it:
       ```matlab
       >> seq2ge;
       ```
   This will create `.pge` and `.entry` files that can be executed on GE with the 
   [pge2](https://github.com/HarmonizedMRI/SequenceExamples-GE/tree/main/pge2)
   GE interpreter.

On the scanner:
1. Copy the `.entry` files to `/srv/nfs/psd/usr/psd/pulseq/v7/`.
   This path is hardcoded in the `pge2` interpreter.
  As usual, take care to not overwrite existing `.entry` files in that folder.
2. Copy the `.pge` files to `scanner_pge_location`.
3. Run the `pge2` sequences. 

