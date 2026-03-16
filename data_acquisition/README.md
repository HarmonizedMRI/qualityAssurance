# This is the instruction for data acquisition.
* `write_QA_Tran_EPIrs.m`: to generate `QA_epi.seq` file for EPI scans.
* `write_QA_Tran_T1.m`: to generate QA_T1.seq file for spin-echo scans.
* `20241122_QA_protocol_instruction_siemens.docx`: standard operating procedure for QA measurements.
* `QA_record.xlsx`: Excel sheet for the record of QA measurements.
* `QA_MAGMA.pdf`: product and Pulseq-based sequence protocols on Cima.X for quality assurance.

## GE users

1. Set `vendor` to `'GE'` in `setvendor.m`
2. Execute the `write*.m` scripts.
3. Set the desired options in `seg2ge.m` and run it.
   This will create `.pge` and `.entry` files that can be executed on GE with the 
   [pge2](https://github.com/HarmonizedMRI/SequenceExamples-GE/tree/main/pge2)
   GE interpreter.

