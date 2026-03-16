# This is the instruction for data acquisition.
* `write_QA_Tran_EPIrs.m`: to generate `QA_epi.seq` file for EPI scans.
* `write_QA_Tran_T1.m`: to generate QA_T1.seq file for spin-echo scans.
* `20241122_QA_protocol_instruction_siemens.docx`: standard operating procedure for QA measurements.
* `QA_record.xlsx`: Excel sheet for the record of QA measurements.
* `QA_MAGMA.pdf`: product and Pulseq-based sequence protocols on Cima.X for quality assurance.

## GE users

In MATLAB:
1. Set `vendor` to `'GE'` in `setvendor.m`
2. Execute the `write*.m` scripts.
3. Set the desired options in `seq2ge.m`:
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
3. Run the `pge2` sequences. Choose `Axial` scan orientation.

