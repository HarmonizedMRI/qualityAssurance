# QA Data Acquisition

This directory contains scripts and documentation for generating the Pulseq
sequences used in the **HarmonizedMRI quality assurance (QA) protocol**.

The MATLAB scripts generate Pulseq (`.seq`) files for phantom QA scans.
These sequences can be:

* executed directly on **Siemens scanners**, or
* converted to **GE format** using the `pge2` interpreter.

---

# Contents

| File                                            | Description                                                          |
| ----------------------------------------------- | -------------------------------------------------------------------- |
| `write_QA_Tran_EPIrs.m`                         | Generates `QA_epi_final.seq` for EPI QA scans                        |
| `write_QA_Tran_T1.m`                            | Generates `QA_T1_final.seq` for spin-echo QA scans                   |
| `20241122_QA_protocol_instruction_siemens.docx` | Standard operating procedure for QA measurements on Siemens scanners |
| `QA_record.xlsx`                                | Spreadsheet for recording QA measurements                            |
| `QA_MAGMA.pdf`                                  | QA protocol documentation for Siemens Cima.X                         |

---

# Workflow Overview

1. Generate Pulseq sequence files using MATLAB.
2. Run the sequences on the MRI scanner.

   * Siemens: run `.seq` files directly.
   * GE: convert `.seq` → `.pge` using `seq2ge`.
3. Acquire QA scans.
4. Record measurements in `QA_record.xlsx`.

---

# Setup

## Python Environment

The scripts call `mr.makeSLRpulse`, which requires Python.

Example setup on Linux:

```bash
sudo apt install python3.12-venv   # or whichever python version you have available
python3 -m venv qa_env
source qa_env/bin/activate
pip install matplotlib scipy sigpy
```

Start MATLAB from the same environment:

```bash
matlab
```

---

## MATLAB Paths

From inside MATLAB (in the repository root):

```matlab
setup
```

This adds the required directories to the MATLAB path.

---

# Generating Pulseq Sequences

1. Configure scanner settings in

    ```
    set_vendor_and_system_limits.m
    ```

    Set the scanner vendor:

    ```matlab
    vendor = 'Siemens';   % or 'GE'
    ```

    Adjust other scanner-specific limits if necessary.

2. Generate the sequences:

    ```matlab
    write_QA_Tran_EPIrs
    write_QA_Tran_T1
    ```

    This produces the Pulseq files required for the QA scans:

    ```
    QA_epi_final.seq
    QA_T1_final.seq
    ```

---

# Running on Siemens Scanners

The generated `.seq` files can be executed directly using the Pulseq interpreter.

Follow the instructions in:

```
20241122_QA_protocol_instruction_siemens.docx
```

for the complete scanner procedure.

---

# Running on GE Scanners

GE scanners require conversion to the `pge2` format.

## 1. Configure Coil Model

Edit `set_vendor_and_system_limits.m`:

```matlab
coil = 'xrm';   % example for MR750
```

This is used for PNS checking.

---

## 2. Convert Pulseq Files

Edit options in `seq2ge.m`:

* `scanner_pge_location`: location on the scanner where `.pge` files will be stored
* `opuser1`: controls the generated `.entry` filenames

Run:

```matlab
seq2ge
```

This generates:

```
.pge
.entry
```

files that can be executed on GE scanners.

---

## 3. Copy Files to the Scanner

```
.entry → /srv/nfs/psd/usr/psd/pulseq/v7/
.pge   → scanner_pge_location
```

⚠️ **Important**

Do **not overwrite existing `.entry` files** in the Pulseq directory.

---

## 4. Run the Sequences

Execute the sequences on the scanner using the **pge2 interpreter**, with the following settings:

* EPI scan: 200 runs (opnex = 200)    
* T1 scan: 2 runs (opnex = 2)

---

# Dependencies

* MATLAB
* Python (for SLR pulse design)
* Pulseq
* `pge2` GE interpreter
  https://github.com/HarmonizedMRI/SequenceExamples-GE/tree/main/pge2

---

# Tested Systems

The protocol has been tested on:

* Siemens Cima.X
* GE MR750 (with `pge2` interpreter)

---

# QA Data Recording

QA measurements should be recorded in:

```
QA_record.xlsx
```

following each acquisition session.
