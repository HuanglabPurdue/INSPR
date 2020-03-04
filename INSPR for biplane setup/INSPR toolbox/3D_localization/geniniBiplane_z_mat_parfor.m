% Script for estimating initial axial position
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, August 2019
%%  
function [z_guess,mcc_val,psf_model] = geniniBiplane_z_mat_parfor(data,psfobj_all,psftype,zLim,Nz,planedis,tform)

R = size(data,1);
Nfit = size(data,3);
zi = linspace(zLim(1),zLim(2),Nz);
Nplane = size(data,4);
psf_model = zeros(R,R,Nz,Nplane);  

for ss = 1:Nplane
    psfobj = psfobj_all{ss};
    psfobj.Xpos = zeros(Nz,1);% pixel
    psfobj.Ypos = zeros(Nz,1);% pixel
    switch psftype
        case 'normal'
            psfobj.Zpos = zi+planedis(ss);% micron
            psfobj.genPSF_2();  %FX
            psfobj.scalePSF('normal');
        case 'IMM'
            psfobj.ZposMed = zi;% micron
            psfobj.genIMMPSF();
            psfobj.scalePSF('IMM');
    end
    
    if ss == 2
        edge_1 = mean(squeeze(psfobj.ScaledPSFs(1,:,:)),1)';
        edge_2 = mean(squeeze(psfobj.ScaledPSFs(:,1,:)),1)';
        edge_3 = mean(squeeze(psfobj.ScaledPSFs(end,:,:)),1)';
        edge_4 = mean(squeeze(psfobj.ScaledPSFs(:,end,:)),1)';
        edge_fill = min([edge_1 edge_2 edge_3 edge_4],[],2);
        
        psf_model(:,:,:,ss) = imwarp(psfobj.ScaledPSFs,tform,'cubic','OutputView',imref2d(size(psfobj.ScaledPSFs)),...
            'FillValues',edge_fill);
    else 
        psf_model(:,:,:,ss) = psfobj.ScaledPSFs;
    end
end
   

z_guess = zeros(Nfit,1);
mcc_val = zeros(Nfit,1);
img_fft2 = zeros(R,R,Nfit,Nplane);
ref_fft2 = zeros(R,R,Nz,Nplane);

for ii = 1: Nz
    %cal ref fft
    for ss = 1:Nplane
        ref = psf_model(:,:,ii,ss);
        ref = ref-mean(ref(:));
        ref = ref/std(ref(:));
        ref_fft2(:,:,ii,ss) = fft2(ref);
    end
end

parfor nn = 1:Nfit
   %cal img fft
   for ss = 1:Nplane
       img = data(:,:,nn,ss);
       img = img-mean(img(:));
       img = img/std(img(:));
       img_fft2(:,:,nn,ss) = fft2(img);
   end 
end

parfor nn = 1:Nfit
    mcc_max = 0;
    ind = 1;
    for ii = 1:Nz
        mcc = 0;
        for ss = 1:Nplane
            cc_value = abs(ifft2(ref_fft2(:,:,ii,ss) .*conj(img_fft2(:,:,nn,ss))));
            maxa = 1/(R*R)*max(cc_value(:));
            mcc = mcc + maxa/Nplane;
        end
        if mcc>mcc_max
            ind = ii;
            mcc_max = mcc;
        end
    end
    z_guess(nn) = zi(ind);
    mcc_val(nn) = mcc_max;
end


end