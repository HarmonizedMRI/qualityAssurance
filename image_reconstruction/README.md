# This is the instruction for image reconstruction.

## Documents
### MRD conversion
* `pulseq2mrd_epi.m`: convert GE Pulseq EPI raw data (`.mat`) to MRD raw data (`.h5`) using the LABEL information in the `QA_epi.seq` file.
* `pulseq2mrd_se.m`: convert GE Pulseq spin-echo raw data (`.mat`) to MRD raw data (`.h5`) using the LABEL information in the `QA_T1.seq` file.
* `siemens2mrd_epi.m`: convert Siemens Pulseq EPI raw data (`.dat`) to MRD raw data (`.h5`) using the information in the `.dat` raw data itself.

### Gadgetron configuration files
* `default.xml`: Gadgetron configuration file for spin-echo image reconstruction. (This document is already in the Gadgetron container: `/opt/conda/envs/gadgetron/share/gadgetron/config/default.xml`)
* `qc_epi.xml`: Gadgetron configuration file for EPI image reconstruction, which is modified from the default `epi.xml` located in the Gadgetron container: `/opt/conda/envs/gadgetron/share/gadgetron/config/`.

### ICE online reconstruction
* `specialCard_ICE.png`: Special card setting for ICE online reconstruction.

## Procedures for Gadgetron offline reconstruction
### Step 1: Gadgetron installation ([for more details](https://gadgetron.github.io/tutorial/))
* Download and install [Docker](https://www.docker.com/) software.
* Open your terminal and run: `docker run -t --name gt_lastest --detach --volume %cd%:/opt/data ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_dev_nocuda:latest`. This will download and then launch the [latest-version Gadgetron](https://gadgetron.readthedocs.io/en/latest/building.html) in a Docker container. It will also mount your current folder as a data folder inside the container.
* Run this command: `docker exec -ti lastest /bin/bash`. This will execute your Gadgetron container.
### Step 2
* Place your SE/EPI .dat/.h5 data in the mounted folder.
* Run the command in Terminal: `cd /opt/data` to enter the mounted folder.
### Step 3: MRD conversion
* for Siemens data, you can convert the `.dat` data to MRD data by using Gadgetron. If Gsdgetron doesn't work (e.g. for XA EPI data), you can then use the Matlab script, `siemens2mrd_epi.m`.
* Command for spin-echo data conversion: `siemens_to_ismrmrd -f meas_MID*.dat -z 2 -o se_data.h5`.
* Command for EPI data conversion: `siemens_to_ismrmrd -f meas_MID*.dat -z 2 -m IsmrmrdParameterMap_Siemens.xml -x IsmrmrdParameterMap_Siemens_EPI.xsl -o epi_data.h5`.
* for GE data, you can convert the `.mat` raw data to MRD data by using the Matlab scripts with the corresponding `.seq` files. spin-echo conversion: use `pulseq2mrd_se.m` with `QA_T1.seq`. EPI conversion: use `pulseq2mrd_epi.m` with `QA_epi.seq`.
### Step 4: Gadgetron reconstruction
* SE reconstruction: `gadgetron_ismrmrd_client -f se_data.h5 -c default.xml -o se_out.h5`.
* EPI reconstruction: first, put qc_epi.xml to the mounted folder and then copy it to Gadgetron container: `cp /opt/data/qc_epi.xml /opt/conda/envs/gadgetron/share/gadgetron/config/`. Then, run reconstruction: `gadgetron_ismrmrd_client -f epi_data.h5 -c qc_epi.xml -o epi_out.h5`.
### Step 5: Load Gadgetron-reconstructed images (`.h5`)
* Load SE `.h5` images             
`filename = 'pulseq_se_out.h5' ;`        
`info = hdf5info(filename) ;`        
`address_data_1 = info.GroupHierarchy.Groups(1).Groups.Datasets(2).Name ;`     
`pulseq_se_im = squeeze(double( hdf5read(filename, address_data_1) ) ) ;`      
`pulseq_se_im = reshape(pulseq_se_im, [256, 256, 11, 2]) ;`
* Load EPI `.h5` images:             
`filename = 'pulseq_epi_out.h5' ;`
`info = hdf5info(filename) ;`         
`address_data_1 = info.GroupHierarchy.Groups(1).Groups.Datasets(2).Name ;`         
`pulseq_epi_im = squeeze(double( hdf5read(filename, address_data_1) ) ) ;`            
`pulseq_epi_im = reshape(pulseq_epi_im, [64, 64, 27, 200]) ;`         
### Procedure for ICE online reconstruction
Before executing the Pulseq-based sequences, you can enable ICE online Reconstruction following the procedure below:
* Navigate to the Special Card (`specialCard_ICE.png`), set `Data handling` to `ICE STD` for NUMARIS/X (e.g. XA60A, XA61A), and `ICE 2D` for NUMARIS/4 (e.g. VB, VD, VE).
* Select `Sum-of-Square` for coil combination.
* Be sure that the maximal pixel intensity does not violate the intensity threshold of **4096**.
