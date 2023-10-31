# Pipeline-GATK_CNV

![image](https://github.com/cmasotti-Lab/Pipeline_Gatk-CNV/assets/11162991/71bfc5e6-c512-4db2-bc8d-5e0e06231686)
[Detailed Extended Pipelines](https://drive.google.com/file/d/100eEe_oiofVWKpySCTJfEtS2OvwRBM9m/view?usp=sharing)

This pipeline has been developed to identify somatic copy number variations (CNV) in tumor-only exome data.

The samples used are from the locally advanced rectal cancer project.

The pipeline is based on the GATK CNV approach for identifying somatic mutations in tumor samples in the presence of normal controls.

GATK CNV uses a Panel of Normals (PoN) with unrelated samples to separate germline variants from somatic ones.

To develop the PoN, we used samples from 100 non-cancer individuals from the SARS-CoV2-Brasil project (SECOLIN, et al. 2021).


## Pipeline Stages

This pipeline is divided into two stages:

- Creation of the Panel of Normals (PoN).
- Identification of CNV.

[Detailed Extended Pipelines](https://drive.google.com/file/d/100eEe_oiofVWKpySCTJfEtS2OvwRBM9m/view?usp=sharing)
