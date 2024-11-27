# This is the instruction for post-processing.
The example post-processing is based on the reconstructed images from Cima.X over five days.
## Reconstructed images from Cima.X
**Note**: The images are too big to upload here, and can be accessed via [Zenodo](https://doi.org/10.5281/zenodo.14217778) as an alternative. All `se` folders contain a `structuralQuality_main.m` to call the `structuralQuality.m` function for structural quality analysis. All `epi` folders contain a `temporalQuality_main.m` to call the `temporalQuality.m` function for temporal quality analysis.
* `product_epi_ice`: ICE-reconstructed product EPI images.
* `product_epi_gt`: Gadgetron-reconstructed product EPI images.
* `pulseq_epi_ice`: ICE-reconstructed Pulseq EPI images.
* `pulseq_epi_gt`: Gadgetron-reconstructed Pulseq EPI images.
* `product_se_ice`: ICE-reconstructed product SE images.
* `product_se_gt`: Gadgetron-reconstructed product SE images.
* `pulseq_se_ice`: ICE-reconstructed Pulseq SE images.
* `pulseq_se_gt`: Gadgetron-reconstructed Pulseq SE images.
## QA analysis Matlab package: `QA_functions`
* `circfit.m`: to find the center point and radius of the phantom.
* `makeCircleMask.m`: to make a circular mask based on the center point and radius.
* `structuralQuality.m`: to analyze the structural quality of the SE images.
* `temporalQuality.m`: to analyze the temporal quality of the EPI images.
### Post-processing procedures
* Step 1: Add the `QA_functions` folder to your Matlab Path.
* Step 2: Run the `temporalQuality_main.m` or `structuralQuality_main.m` script in each folder to produce the QA results of all reconstructed images inside the folder.
* Step 3: Run the `make_figure_epi.m` and `make_figure_se.m` to produce some of the tables and figures used in the manuscript.
