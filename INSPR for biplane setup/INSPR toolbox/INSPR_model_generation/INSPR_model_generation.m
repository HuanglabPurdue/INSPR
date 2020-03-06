%%
% Script for estimating in situ 3D PSF  
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, August 2019
%           
function probj = INSPR_model_generation(subregion_ch1,subregion_ch2,setup,pupil_para)

global pupil_stop

empupil = [];                            
%% setup parameters     
empupil.NA = setup.NA;
empupil.Lambda = setup.Lambda;
empupil.nMed = setup.nMed;
empupil.Pixelsize = setup.Pixelsize;   
%When index mismatch: Dist' = Dist * sample_med / immersion_med
empupil.biplane_dist = setup.biplane_dist * setup.nMed / setup.RefractiveIndex;   %um 
                            
%% pupil parameters
empupil.imsz = size(subregion_ch1,1);
empupil.Z_pos = pupil_para.Z_pos;
empupil.bin_lowerBound = pupil_para.bin_lowerBound;
empupil.min_similarity = pupil_para.min_similarity;
empupil.iter = pupil_para.iter;
empupil.blur_sigma = pupil_para.blur_sigma;
empupil.init_z = pupil_para.init_z;     %initial zernike value in pupil phase
empupil.ZernikeorderN = pupil_para.ZernikeorderN; 
empupil.Zernike_sz = (empupil.ZernikeorderN+1)^2;
empupil.Zshift_mode = pupil_para.Zshift_mode; 
empupil.XYshift_mode = pupil_para.XYshift_mode;     
empupil.iter_mode = 0;  


% initial Z 
empupil.zshift = 0;

%% Generate initial pupil
probj = gen_initPupil(empupil);

%% Iteration 

for iter = 1 : empupil.iter
display(['Iteration: ' num2str(iter) ' ...']);


%% Generate reference Z-postions' PSFs in two channels from pupil function
[ref_plane1, ref_plane2] = gen_aberPSF_fromPR(probj, empupil);



%% classify single molecules of two channels to giving reference Z-position images, and 
%  average these Z-postion images in a certain group 
%  calculate similartiy between single molecules and reference images, do
%  X-Y registration
drawnow
if pupil_stop == 1
    return;
end

if empupil.XYshift_mode == 0
    [ims_Zcal_ave_plane1, ims_Zcal_ave_plane2, index_record_Zplanes] = classify_twoPlanes_par(subregion_ch1, subregion_ch2, ...
        ref_plane1, ref_plane2, empupil);
else
    [ims_Zcal_ave_plane1, ims_Zcal_ave_plane2, index_record_Zplanes] = classify_twoPlanes_together_par(subregion_ch1, subregion_ch2, ...
        ref_plane1, ref_plane2, empupil);
end

%% Refine the pupil function and estimate the aberration using average 
%  images of two channels in different Z postions

drawnow
if pupil_stop == 1
   return;
end

[probj,empupil] = PRPSF_aber_fromAveZ(ims_Zcal_ave_plane1,ims_Zcal_ave_plane2,index_record_Zplanes,empupil);

end

index_number = length(find(index_record_Zplanes==1));

if index_number <= 5
    msgbox('Warning! The range of localization isn’t enough for reliable model generation!');
end

