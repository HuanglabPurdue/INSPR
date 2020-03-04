% Script for calculating affine matrix in model
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, August 2019
%       
%%
function [tform_fitting, tMatrix_noTranslation] = cal_model_affine(tform,ref_point)

invertform = invert(tform);

[x,y] = transformPointsForward(invertform,ref_point(1),ref_point(2));

x_offset = ref_point(1)-x;
y_offset = ref_point(2)-y;



tform_fitting = invertform;
tform_fitting.T(3,1) = tform_fitting.T(3,1) + x_offset;
tform_fitting.T(3,2) = tform_fitting.T(3,2) + y_offset;

tMatrix_noTranslation = single(zeros(3,2));
tMatrix_noTranslation(1,1) = tform_fitting.T(1,1);
tMatrix_noTranslation(1,2) = tform_fitting.T(1,2);
tMatrix_noTranslation(2,1) = tform_fitting.T(2,1);
tMatrix_noTranslation(2,2) = tform_fitting.T(2,2);
tMatrix_noTranslation(3,1) = tform_fitting.T(3,1);
tMatrix_noTranslation(3,2) = tform_fitting.T(3,2);

