% Script for calculating 2D cross correlation
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
function [maxa,cc] = cc2(ref,img)
ref = ref-mean(ref(:));
ref = ref/std(ref(:));
img = img-mean(img(:));
img = img/std(img(:));
cc = abs(ifft2(fft2(ref).*conj(fft2(img))));
maxa = 1/numel(img)*max(cc(:));
end