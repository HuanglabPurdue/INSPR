% Script for generating templates from a given pupil
%
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, August 2019
%       
%
function [ref_plane1, ref_plane2] = gen_aberPSF_fromPR(probj, empupil)


%% reference Z image in plane1 and plane2

disp('Generate reference PSFs in different Z positions of two channels from pupil function');
imsz = empupil.imsz;  
Numimage = size(empupil.Z_pos,2);   
ref_plane1 = zeros(imsz,imsz,Numimage);
ref_plane2 = zeros(imsz,imsz,Numimage);


num = 1;
for ii = 1 : Numimage
    
    %% calculate real Z image in plane1 and plane2
    Ztrue_plane1 = empupil.Z_pos(ii) - empupil.biplane_dist/2 + empupil.zshift;
    Ztrue_plane2 = empupil.Z_pos(ii) + empupil.biplane_dist/2 + empupil.zshift;
    
    
    %% get pupil from Phase retrieval
    PRstruct = probj.PRstruct;
    pxsz = probj.Pixelsize;
    R = 128;
    PRstruct.SigmaX = empupil.blur_sigma;
    PRstruct.SigmaY = empupil.blur_sigma;
    %% generate reference Z images in Plane1
    psfobj = PSF_zernike(PRstruct);
    psfobj.Xpos = 0;
    psfobj.Ypos = 0;
    psfobj.Zpos = Ztrue_plane1;   

    psfobj.Boxsize = imsz;
    psfobj.Pixelsize = pxsz; % micron
    psfobj.PSFsize = R;
    psfobj.nMed = probj.PRstruct.RefractiveIndex;
    
    psfobj.precomputeParam();
    psfobj.setPupil();
    psfobj.genPSF_2();
    psfobj.scalePSF('normal');
    
    %% generate reference Z images in Plane2
    psfobj1 = PSF_zernike(PRstruct);
    psfobj1.Xpos = 0;
    psfobj1.Ypos = 0;
    psfobj1.Zpos = Ztrue_plane2;  

    psfobj1.Boxsize = imsz;
    psfobj1.Pixelsize = pxsz; % micron
    psfobj1.PSFsize = R;
    psfobj1.nMed = probj.PRstruct.RefractiveIndex;
    
    psfobj1.precomputeParam();
    psfobj1.setPupil();
    psfobj1.genPSF_2();
    psfobj1.scalePSF('normal');

    I = permute(repmat(5000, [1, imsz, imsz]), [2,3,1]);
  
    psf = psfobj.ScaledPSFs;

    normf = sum(sum(psfobj.Pupil.mag.^2,1),2);
    img = (psf./normf).*I;
    
    psf1 = psfobj1.ScaledPSFs;

    normf1 = sum(sum(psfobj1.Pupil.mag.^2,1),2);
    img1 = psf1./normf1.*I;
   
    bg = 0;
    
    ref_plane1(:,:,num) = sum(img, 3)+ bg;
    ref_plane2(:,:,num) = sum(img1, 3)+ bg;
    num = num + 1;
end



