% Script for realigning axial position
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, August 2019
% 
function [ims_Zcal_align,Zpos] = realign_Z_ave_img(ims_Zcal_ave,Zpos_in,probj)


imsz = size(ims_Zcal_ave,1);
Numimage = size(ims_Zcal_ave,3);  
ims_Zref = zeros(imsz,imsz,Numimage);

%% Re-arrange Z position 
num = 1;
Zpos = [];
while num <= Numimage
    if num <= Numimage -1 && abs(Zpos_in(num)-Zpos_in(num+1)) < 0.04
       Zpos(num) = (Zpos_in(num)+Zpos_in(num+1)) / 2;
       Zpos(num+1) = Zpos(num);
       num = num + 1;
    else 
       Zpos(num) = Zpos_in(num);
    end
    num = num + 1;
end

%% generate Reference image
for ii = 1 : Numimage
    
    % get pupil from Phase retrieval
    PRstruct = probj.PRstruct;
    pxsz = probj.Pixelsize;
    R = 128;

    % generate reference Z images
    psfobj = PSF_zernike(PRstruct);
    psfobj.Xpos = 0;
    psfobj.Ypos = 0;
    psfobj.Zpos = Zpos(ii);    %micron, imporant

    psfobj.Boxsize = imsz;
    psfobj.Pixelsize = pxsz; % micron
    psfobj.PSFsize = R;
    psfobj.nMed = probj.PRstruct.RefractiveIndex;
    
    psfobj.precomputeParam();
    psfobj.setPupil();
    psfobj.genZernikeMag();
    psfobj.genPSF_2();
    psfobj.scalePSF('normal');
    
    % add I, bg
    I = permute(repmat(5000, [1, imsz, imsz]), [2,3,1]);
  
    psf = psfobj.ScaledPSFs;
    normf = sum(sum(psfobj.Pupil.mag.^2,1),2);
    img = (psf./normf).*I;
    
    bg = 0;
    
    ims_Zref(:,:,ii) = sum(img, 3)+ bg;
end

%% Re-alignment Z average imagei
shift_backup = [];
ims_Zcal_align = zeros(imsz,imsz,Numimage);
for ii = 1 : Numimage
    ref_tmp = ims_Zref(:,:,ii);
    reg_tmp = ims_Zcal_ave(:,:,ii);
    
    [output, ~] = dftregistration(fft2(ref_tmp),fft2(reg_tmp),20);
    
    shift_backup = cat(1,shift_backup,output(3:4));
    ims_Zcal_align(:,:,ii) = FourierShift2D(reg_tmp, [output(3) output(4)]);
end

%% Merge same Z group
num = 1;
while num < length(Zpos)    
    if Zpos(num) == Zpos(num+1)
        Zpos(num+1) = [];
        ims_Zcal_align(:,:,num) = (ims_Zcal_align(:,:,num) + ims_Zcal_align(:,:,num+1)) /2; 
        ims_Zcal_align(:,:,num+1) = [];
    end
    num = num + 1;
end


