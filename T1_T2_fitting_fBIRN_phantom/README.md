# Images and scripts for T1/T2 fitting
This package includes DICOM images and T1/T2 fitting scripts for the fBIRN phantom. Images for T1 fitting were acquired using a product turbo spin echo sequence with an inversion recovery pulse (repetition time = 4000 ms, echo train length = 4). Images for T2 fitting were obtained using a product SE sequence (repetition time = 3500 ms). Both measurements were conducted on the Siemens Prisma.Fit 3T scanner on 05.06.2024.                       

The package contains the following sub-folders/scripts:
* T1 sub-folder: contains all DICOM images for T1 fitting with inversion recovery times of {50, 150, 300, 450, 600, 750, 900, 1050, 1200, 1350, 1500, 2200, 3000} ms.
* T2 sub-folder: contains all DICOM images for T2 fitting with echo times of {7.5, 15, 30, 45, 60, 75, 90, 130, 200, 250} ms.
* Do_T1fit.m: Matlab script for T1 fitting.
* Do_T2fit.m: Matlab script for T2 fitting.
