# Example raw data and DICOM images from Cima.X.
**Note**: The package is too big to upload here and can be accessed via this link instead.
* DICOM folder: contains the DICOM images for four EPI scans (the first two scans for warm-up) and two spin-echo scans.
* `meas*.dat`: Siemens raw data of two EPI scans for temporal quality analysis and two spin-echo scans for structural quality analysis.
* `*data.h5` files: the ISMRMRD data of the four raw datasets.
* `*out.h5` files: the images reconstructed by Gadgetron.
* `*.nii`: the nifti-format reconstructed images.
* `siemens2mrd_epi.m`: to convert the Siemens EPI raw data to ISMRMRD data.
* `read_image.m`: to convert the Gadgetron-reconstructed h5-format images to nifti-format images.