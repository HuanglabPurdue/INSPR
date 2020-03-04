% Script for calculating image shift, version 2
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, August 2019
% 
function [shift1, shift2, similarity]=registration_in_biplane(ref1,ref2,img1,img2)
% normalization
ref1 = ref1-mean(ref1(:));
ref1 = ref1/std(ref1(:));

ref2 = ref2-mean(ref2(:));
ref2 = ref2/std(ref2(:));

img1 = img1-mean(img1(:));
img1 = img1/std(img1(:));

img2 = img2-mean(img2(:));
img2 = img2/std(img2(:));

output = dftregistration_biplane(fft2(ref1),fft2(img1),fft2(ref2),fft2(img2),10);
shift1 = output(3);
shift2 = output(4);

img1_shift = FourierShift2D(img1, [shift1 shift2]); 
img2_shift = FourierShift2D(img2, [shift1 shift2]); 


[similarity1,cc_1] = cc2(ref1(3:end-2,3:end-2),img1_shift(3:end-2,3:end-2));
[similarity2,cc_2] = cc2(ref2(3:end-2,3:end-2),img2_shift(3:end-2,3:end-2));

cc_value = cc_1 + cc_2;
similarity = 1/numel(cc_value)*max(cc_value(:))/2;
