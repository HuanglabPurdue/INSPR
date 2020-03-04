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
function errM = gen_LLR(data,psfobj_input,PM)

imsz = size(data,1);
Nfit = size(data,3);

psfobj = psfobj_input;
psfobj.Xpos = PM(1:Nfit,1);% pixel
psfobj.Ypos = PM(1:Nfit,2);% pixel
psfobj.Zpos = PM(1:Nfit,3);% micron
psfobj.genPSF_2();
psfobj.scalePSF('normal');

psf = psfobj.ScaledPSFs;
normf = sum(sum(psfobj.Pupil.mag.^2,1),2);
I = permute(repmat(PM(1:Nfit,4),[1, imsz, imsz]),[2,3,1]);
bg = permute(repmat(PM(1:Nfit,5),[1, imsz, imsz]),[2,3,1]);
psf_model = (psf./normf).*I + bg;
        

errM = [];
psf_input = reshape(psf_model, imsz*imsz,Nfit);
data_input = reshape(data, imsz*imsz,Nfit);

sse = sum((psf_input-data_input).^2);
LLR = 2*(sum(psf_input - data_input - data_input.*log(psf_input) + data_input.*log(data_input)));

errM(:,1) =  sse';
errM(:,2) =  LLR';






