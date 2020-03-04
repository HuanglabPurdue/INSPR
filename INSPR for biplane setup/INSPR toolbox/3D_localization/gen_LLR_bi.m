% Script for calcaluate the log-likelihood ratio (LLR) and error sum of squares (SSE)
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, December 2019
%%
function errM = gen_LLR_bi(data,psfobj_all,PM,planedis,tform_fitting,offset_seg)

imsz = size(data,1);
Nfit = size(data,3);
Nplane = size(data,4);
psf_model = zeros(imsz,imsz,Nfit,Nplane);  

for ss = 1:Nplane
    psfobj = psfobj_all{ss};
    psfobj.Xpos = PM(1:Nfit,1);% pixel
    psfobj.Ypos = PM(1:Nfit,2);% pixel   
    psfobj.Zpos = PM(1:Nfit,3) + planedis(ss);% micron
    psfobj.genPSF_2();  
    psfobj.scalePSF('normal');
            
    psf = psfobj.ScaledPSFs;
    normf = sum(sum(psfobj.Pupil.mag.^2,1),2);
    
    I = permute(repmat(PM(1:Nfit,3+ss),[1, imsz, imsz]),[2,3,1]);
    bg = permute(repmat(PM(1:Nfit,5+ss),[1, imsz, imsz]),[2,3,1]);
    
    img = (psf./normf).*I + bg;
    
    if ss == 2     
        for ii = 1 : Nfit        
            edge_1 = mean(img(1,:,ii));
            edge_2 = mean(img(1,:,ii));
            edge_3 = mean(img(1,:,ii));
            edge_4 = mean(img(1,:,ii));
            edge_fill = min([edge_1 edge_2 edge_3 edge_4]);
            
            tform_tmp = tform_fitting;
            tform_tmp.T(3,1) = tform_fitting.T(3,1) + offset_seg(ii,1);
            tform_tmp.T(3,2) = tform_fitting.T(3,2) + offset_seg(ii,2);
            
            psf_model(:,:,ii,ss) = imwarp(img(:,:,ii),tform_tmp,'cubic','OutputView',imref2d(size(img)),...
                'FillValues',edge_fill);
        end
    else 
        psf_model(:,:,:,ss) = img;
    end
end

LLR = zeros(1,Nfit);
sse = zeros(1,Nfit);
errM = [];
for ss = 1:Nplane
   psf_input = reshape(psf_model(:,:,:,ss), imsz*imsz,Nfit);
   data_input = reshape(data(:,:,:,ss), imsz*imsz,Nfit);
   
   sse = sse + sum((psf_input-data_input).^2);
   LLR = LLR + 2*(sum(psf_input - data_input - data_input.*log(psf_input) + data_input.*log(data_input)));
end

errM(:,1) =  sse';
errM(:,2) =  LLR';






