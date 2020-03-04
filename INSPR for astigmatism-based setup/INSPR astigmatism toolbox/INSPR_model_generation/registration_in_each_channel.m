% Script for calculating image shift
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, October 2019
% 
function [shift1, shift2, similarity]=registration_in_each_channel(im1,im2)


output = dftregistration(fft2(im1),fft2(im2),10);

shift1 = output(3);
shift2 = output(4);

im2_shift = FourierShift2D(im2, [shift1 shift2]); 

[similarity,cc] = cc2(im1(3:end-2,3:end-2),im2_shift(3:end-2,3:end-2));