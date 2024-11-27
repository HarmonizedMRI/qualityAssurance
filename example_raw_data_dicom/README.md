# Example raw data and DICOM images
The data and DICOM images were acquired from Cima.X on the fBIRN phantom on 06.08.2024.              
**Note**: The package is too big to upload here and can be accessed via [Zenodo](https://doi.org/10.5281/zenodo.14217778) instead.
* DICOM folder: contains the DICOM images for four EPI scans (the first two scans for warm-up) and two spin-echo scans.
* `meas*.dat`: Siemens raw data of two EPI scans for temporal quality analysis and two spin-echo scans for structural quality analysis.
* `*data.h5` files: the ISMRMRD data of the four raw datasets.
* `*out.h5` files: the images reconstructed by Gadgetron.
* `*.nii`: the NIFTI-format reconstructed images.
* `siemens2mrd_epi.m`: to convert the Siemens EPI raw data to ISMRMRD data.
* `read_image.m`: to convert the Gadgetron-reconstructed h5-format images to NIFTI-format images.
