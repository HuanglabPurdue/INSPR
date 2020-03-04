% Script for estimating initial axial position
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, October 2019
%%  
function [z_guess,mcc_val,psf_model] = genini_z_mat_parfor(data,psfobj_input,psftype,zLim,Nz)

R = size(data,1);
Nfit = size(data,3);
zi = linspace(zLim(1),zLim(2),Nz);
psf_model = zeros(R,R,Nz);  

% generate template PSFs along axial position
psfobj = psfobj_input;
psfobj.Xpos = zeros(Nz,1);% pixel
psfobj.Ypos = zeros(Nz,1);% pixel
switch psftype
    case 'normal'
        psfobj.Zpos = zi;% micron
        psfobj.genPSF_2();
        psfobj.scalePSF('normal');
      
    case 'IMM'
        psfobj.ZposMed = zi;% micron
        psfobj.genIMMPSF();
        psfobj.scalePSF('IMM');
end
psf_model(:,:,:) = psfobj.ScaledPSFs;

   
z_guess = zeros(Nfit,1);
mcc_val = zeros(Nfit,1);
img_fft2 = zeros(R,R,Nfit);
ref_fft2 = zeros(R,R,Nz);

for ii = 1: Nz
    %cal ref fft
    ref = psf_model(:,:,ii);
    ref = ref-mean(ref(:));
    ref = ref/std(ref(:));
    ref_fft2(:,:,ii) = fft2(ref);
end

parfor nn = 1:Nfit
    %cal img fft
    img = data(:,:,nn);
    img = img-mean(img(:));
    img = img/std(img(:));
    img_fft2(:,:,nn) = fft2(img);
end

parfor nn = 1:Nfit
    mcc_max = 0;
    ind = 1;
    for ii = 1:Nz
        cc_value = abs(ifft2(ref_fft2(:,:,ii) .*conj(img_fft2(:,:,nn))));
        mcc = 1/(R*R)*max(cc_value(:));            

        if mcc>mcc_max
            ind = ii;
            mcc_max = mcc;
        end
    end
    z_guess(nn) = zi(ind);
    mcc_val(nn) = mcc_max;
end

end