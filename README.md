# INSPR (*in situ* PSF retrieval)
*In situ* point spread function retrieval (INSPR) toolbox is distributed as accompanying software for manuscript: Fan Xu, Donghan Ma, Kathryn P. MacPherson, Sheng Liu, Ye Bu, Yu Wang, Yu Tang, Cheng Bi, Tim Kwok, Alexander A. Chubykin, Peng Yin, Sarah Calve, Gary E. Landreth, and Fang Huang, "Three-dimensional nanoscopy of whole cells and tissues with in situ point spread function retrieval" (2020) **Nature Methods**, 17(5): 531-540, https://www.nature.com/articles/s41592-020-0816-x

INSPR toolbox is developed for both biplane and astigmatism-based setups. It constructs an *in situ* 3D point spread function (PSF) directly from the obtained single molecule dataset and features an easy-to-use user interface including all steps of 3D single molecule localization from INSPR model generation, pupil-based 3D localization (supporting both GPU with cubic spline implementation and CPU versions), drift correction, volume alignment, to super-resolution image reconstruction. It also contains a small single molecule dataset for users to run as an example.

Due to the limited space of Github (<100 MB), we shrink the single-molecule dataset. The data that support the findings of this study are available from the corresponding authors upon request.  

Note: The target of INSPR is to deal with whole cell and tissue specimens. If the specimen is very thin (typically less than 1 µm), the range of localization may not be enough for reliable model generation.

# Installation environment
* Windows 7 or later, 64 bit.
* MATLAB R2016b, 64 bit (downloadable at http://www.mathworks.com).
* CUDA 7.5 compatible graphics driver (downloadable at https://developer.nvidia.com/cuda-75-downloads-archive).

# INSPR toolbox for biplane setup
## Installation of INSPR toolbox
``` Installation
1) Go to the ‘INSPR for biplane setup/INSPR toolbox’ folder and open the ‘main.m’ file.
2) Check the ‘support_path’ in the ‘main.m’ file to make sure that the path of the ‘Support’ folder is correct.
3) Run the ‘main.m’ file.
4) Check the detail manual "INSPR software user manual".
```

## Dataset and source codes
### Demonstration dataset 
* Data\rawData.mat: Single molecule dataset
* Data\ config.ma: General setting parameters
* Data\tform.mat: Alignment calibration file
* Data\subregions.mat: Cropped sub-regions
* Data\probj.mat: In situ 3D PSF model
* Data\sCMOS_calibration.mat: sCMOS calibration parameters
* Data\recon3D.mat: 3D super-resolution reconstruction results


### INSPR source codes
**‘Main’ folder**
* main.m: Main script for running INSPR
* INSPR_GUI.m: Script for INSPR GUI
* INSPR_GUI.fig: INSPR GUI
* default_cfg.mat: Default configuration for INSPR GUI
* genPupilfigs.m: Script for generating figures of retrieved pupil
* export2csv.m: Script for exporting to ‘csv’ format
* srhist_color.m: Script for generating color-coded super-resolution image


**‘Biplane registration’ folder**
* biplane_registration.m: Script for biplane registration


**‘Segmentation’ folder**
* crop_subregion.m: Script for segmentation
* cMakeSubregions.mexw64: Mex function for cropping sub-regions


**‘INSPR model generation’ folder**
* INSPR_model_generation.m: Script for estimating in situ 3D PSF
* gen_initPupil.m: Script for generating initial pupil
* classify_twoPlanes_par.m: Script for classification and 2D alignment, when XY_shift_mode is ‘Separate shift’
* classify_twoPlanes_together_par.m: Script for classification and 2D alignment, when XY_shift_mode is ‘Together shift’
* registration_in_each_channel.m: Script for 2D alignment, when XY_shift_mode is ‘Separate shift’
* registration_in_biplane.m: Script for 2D alignment, when XY_shift_mode is ‘Together shift’
* cc2.m: Script for calculating 2D cross correlation
* PRPSF_aber_fromAveZ.m: Script for estimating pupil
* merge_Z_ave_img.m: Script for merging same axial position
* realign_Z_ave_img.m: Script for realigning axial position
* subregion_normalization.m: Script for normalizing sub-regions


**‘3D localization’ folder**
* analysis3D_fromPupil.m: Script for 3D reconstruction
* crop_subregion_without_transData.m: Script for segmentation
* loc_channel_specific_model.m: Script for pupil-based 3D localization
* loc_channel_specific_model_CPU.m: Script for pupil-based 3D localization: CPU version
* genIniguess.m: Script for estimating initial lateral position
* geniniBiplane_z_mat_parfor.m: Script for estimating initial axial position
* cal_model_affine.m: Script for calculating affine matrix in model
* gensamplepsf_biplane.m: Script for pre-generating channel specific model
* genpsf_biplane_real.m: Script for generating model from pupil
* cuda_channel_specific_model.mexw64: Mex function for 3D localization
* CalDevBi.m: Script for calculating image derivatives: CPU version
* gen_calCRLB_bi.m: Script for calling CRLB generation: CPU version
* CalCRLB_bi.m: Script for calculating CRLB: CPU version
* gen_LLR_bi.m: Script for calculating log-likelihood ratio: CPU version



# INSPR toolbox for astigmatism-based setup
## Installation of INSPR astigmatism toolbox
``` Installation
1) Go to the ‘INSPR for astigmatism-based setup/INSPR astigmatism toolbox’ folder and open the ‘main.m’ file.
2) Check the ‘support_path’ in the ‘main.m’ file to make sure that the path of the ‘Support’ folder is correct.
3) Run the ‘main.m’ file.
4) Check the detail manual "INSPR astigmatism user manual".
```

## Dataset and source codes
### Demonstration dataset 
* Data\rawData.mat: Single molecule dataset
* Data\config.mat: General setting parameters
* Data\subregions.mat: Cropped sub-regions
* Data\probj.mat: In situ 3D PSF model
* Data\sCMOS_calibration_ast.mat: sCMOS calibration parameters
* Data\recon3D.mat: 3D super-resolution reconstruction results


### INSPR source codes
**‘Main’ folder**
* main.m: Main script for running INSPR
* INSPR_ast_GUI.m: Script for INSPR astigmatism GUI
* INSPR_ast_GUI.fig: INSPR astigmatism GUI
* default_cfg.mat: Default configuration for INSPR astigmatism GUI
* genPupilfigs.m: Script for generating figures of retrieved pupil
* export2csv.m: Script for exporting to ‘csv’ format
* srhist_color.m: Script for generating color-coded super-resolution image


**‘Segmentation’ folder**
* crop_subregion_ast.m: Script for segmentation
* cMakeSubregions.mexw64: Mex function for cropping sub-regions


**‘INSPR model generation’ folder**
* INSPR_model_generation_ast.m: Script for estimating in situ 3D PSF
* gen_initPupil.m: Script for generating initial pupil
* classify_onePlane_par.m: Script for classification and 2D alignment
* registration_in_each_channel.m: Script for 2D alignment
* cc2.m: Script for calculating 2D cross correlation
* PRPSF_aber_fromAveZ_ast.m: Script for estimating pupil
* subregion_normalization.m: Script for normalizing sub-regions


**‘3D localization’ folder**
* analysis3D_fromPupil_ast.m: Script for 3D reconstruction
* crop_subregion_var_ast.m: Script for segmentation
* loc_ast_model.m: Script for pupil-based 3D localization
* loc_ast_model_CPU.m: Script for pupil-based 3D localization: CPU version
* genIniguess.m: Script for estimating initial lateral position
* genini_z_mat_parfor.m: Script for estimating initial axial position
* gensamplepsf.m: Script for pre-generating 3D model
* genpsf_real.m: Script for generating model from pupil
* genpsfstruct.m: Script for calculating image gradient
* cuda_ast_model.mexw64: Mex function for 3D localization
* CalDev.m: Script for calculating image derivatives: CPU version
* gen_calCRLB.m: Script for calling CRLB generation: CPU version
* CalCRLB.m: Script for calculating CRLB: CPU version
* gen_LLR.m: Script for calculating log-likelihood ratio: CPU version

# Mex function for 3D localization
Mex functions include pupil-based 3D localization for both biplane and astigmatism-based setup (running on GPU with cubic spline implementation). 

The detail codes are in 'Mex function for 3D localization' folder.

# Citation

Please cite INSPR in your publications if it helps your research:

``` Citation
@article{xu2020three,
  title={Three-dimensional nanoscopy of whole cells and tissues with in situ point spread function retrieval},
  author={Xu, Fan and Ma, Donghan and MacPherson, Kathryn P and Liu, Sheng and Bu, Ye and Wang, Yu and Tang, Yu and Bi, Cheng and Kwok, Tim and Chubykin, Alexander A and others},
  journal={Nature Methods},
  volume={17},
  number={5},
  pages={531--540},
  year={2020},
  publisher={Nature Publishing Group}
}
```

# Updated versions

``` INSPR versions
INSPR 1.1: Add background subtraction using temporal median filter. Note: INSPR supports the background subtraction option in cases with high background. During background subtraction, the statistical properties of the raw detected camera counts will be no longer maintained, it may decrease localization precisions.

```



# Acknowledgements
We would like to thank Karthigeyan Dhanasekaran and Patrick Lusk (Yale University) for sharing the labeling protocol of Nup98 and interpretation of the resolved Nup98 structures. We thank Michael J. Mlodzianoski for initial instrument design, Sha An for help in instrument alignment and sample preparation, and David A. Miller for providing labeling protocols of mitochondria and microtubules. F.X., D.M., S.L., C.B., and F.H. were supported by grants from the NIH (GM119785) and DARPA (D16AP00093). K.P.M. and G.E.L were supported by grants from the NIH (AG051495 and AG050597). Y.B. and S.C. were supported by a grant from the NIH (R01AR071359). Y.W. and P.Y. were supported by a grant from the NIH (1R01EB018659) and Harvard Medical School Dean’s Initiative Grant. 

