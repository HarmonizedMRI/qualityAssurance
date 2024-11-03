# This directory contains information for image reconstruction
* `pulseq2mrd_epi.m`: the Matlab script to convert the EPI rawdata to ISMRMRD data.
* `pulseq2mrd_se.m`: the Matlab script to convert the spin-echo rawdata to ISMRMRD data.
* The EPI images were reconstructed with the configuration file `qc_epi.xml` (modified from the `epi.xml`) using Gadgetron. The command is `gadgetron_ismrmrd_client -f epi_data.h5 -c qc_epi.xml -o epi_out.h5`.
* The spin-echo images were reconstructed with the configuration file `default.xml` using Gadgetron. The command is `gadgetron_ismrmrd_client -f se_data.h5 -c default.xml -o se_out.h5`.

For detailed tutorials of the image reconstruction using Siemens ICE and Gadgetron, please visit [this link](https://github.com/pulseq/Pulseq-Rocks-2023-24-ISMRM-Reproducibility-Challenge/tree/main/image_reconstruction_tutorial).
