# INSPR (in situ PSF retrieval)
In situ point spread function retrieval (INSPR) toolbox is distributed as accompanying software for manuscript: Fan Xu, Donghan Ma, Kathryn P. MacPherson, Sheng Liu, Ye Bu, Yu Wang, Yu Tang, Cheng Bi, Tim Kwok, Alexander A. Chubykin, Peng Yin, Sarah Calve, Gary E. Landreth, and Fang Huang, "Three dimensional nanoscopy of whole cells and tissues with in situ point spread function retrieval" (2020) **Nature Methods**, accept.
INSPR toolbox is developed for both biplane and astigmatism-based setups. It constructs an in situ 3D point spread function (PSF) directly from the obtained single molecule dataset and features an easy-to-use user interface including all steps of 3D single molecule localization from INSPR model generation, pupil-based 3D localization (supporting both GPU with cubic spline implementation and CPU versions), drift correction, volume alignment, to super-resolution image reconstruction. It also contains a small single molecule dataset for users to run as an example.


## Files included in this package
### Content of smNet software (Pytorch)
**Sample Data:**
* train/data.mat: A small training dataset containing simulated PSFs.
* train/label.mat: Labels of the training dataset (ground truth of the parameters)
* train/CRLB.mat: Calculated CRLB used in training
* test/data.mat: A small test dataset
* test/label.mat: Underlying true position to compare with smNet results
* test/CRLB.mat: Underlying true position to compare with smNet results
