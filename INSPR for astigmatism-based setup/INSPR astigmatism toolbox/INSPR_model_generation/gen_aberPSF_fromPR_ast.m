% Script for generating templates from a given pupil
%
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, October 2019
%       
%
function [ref_plane1] = gen_aberPSF_fromPR_ast(probj, empupil)


%% reference Z image

disp('Generate reference PSFs in different Z positions from pupil function');
imsz = empupil.imsz;  
Numimage = size(empupil.Z_pos,2);   
ref_plane1 = zeros(imsz,imsz,Numimage);
label = probj.PRstruct.Zernike_phase(5:25);


num = 1;
for ii = 1 : Numimage
    
    %% calculate template Z image
    Ztrue_plane1 = empupil.Z_pos(ii) + empupil.zshift;
    
    %% get pupil from Phase retrieval
    PRstruct = probj.PRstruct;
    pxsz = probj.Pixelsize;
    R = 128;
    PRstruct.SigmaX = empupil.blur_sigma;
    PRstruct.SigmaY = empupil.blur_sigma;
    
%     if abs(label(1)-empupil.init_z(1)) > 0.4
%         PRstruct.Pupil.phase = zeros(R,R);
%         PRstruct.Pupil.mag = zeros(R,R);
%         phaseZ = zeros(1,25);
% %         phaseZ([5:9]) = label(1:5);%Zernike mode
%         phaseZ([5:25]) = label;%Zernike mode
%         
%         phaseZ(5) = empupil.init_z(1);
%         
%         magZ = zeros(1,25);
%         magZ(1) = 1;
%         PRstruct.Zernike_phase = phaseZ;
%         PRstruct.Zernike_mag = magZ;
%     end
    PRstruct.Pupil.phase = zeros(R,R);
    PRstruct.Pupil.mag = zeros(R,R);
    phaseZ = zeros(1,25);
    
    phaseZ([5:25]) = label;%Zernike mode
    magZ = zeros(1,25);
    magZ(1) = 1;
    PRstruct.Zernike_phase = phaseZ;% this generate aberration
    PRstruct.Zernike_mag = magZ;
    
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

%     if abs(label(1)-empupil.init_z(1)) > 0.4
%         psfobj.genZernike();
%         psfobj.genPupil();
%         psfobj.genPSF_2();
%     else
%         psfobj.setPupil();
%         psfobj.genPSF_2();
%     end
    psfobj.genZernike();
    psfobj.genPupil();
    psfobj.genPSF_2();
    
    psfobj.scalePSF('normal');
    %% add I, bg
    I = permute(repmat(5000, [1, imsz, imsz]), [2,3,1]);
    psf = psfobj.ScaledPSFs;

    normf = sum(sum(psfobj.Pupil.mag.^2,1),2);
    img = (psf./normf).*I;
    
    bg = 0;
    ref_plane1(:,:,num) = sum(img, 3)+ bg;
    num = num + 1;
end



