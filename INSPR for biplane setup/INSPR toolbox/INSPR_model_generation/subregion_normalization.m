% Script for normalizing sub-regions
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, August 2019
% 
function [subregion_ch1_norm, subregion_ch2_norm] = subregion_normalization(subregion_ch1, subregion_ch2)

subregion_ch1_norm = [];
subregion_ch2_norm = [];

for i = 1 : size(subregion_ch1,3)
    tmp_ch1 = subregion_ch1(:,:,i);
    tmp_ch1 = reshape(zscore(tmp_ch1(:)),size(tmp_ch1,1),size(tmp_ch1,2));  %Standardized z-scores
    subregion_ch1_norm(:,:,i) = tmp_ch1;
    
    tmp_ch2 = subregion_ch2(:,:,i);
    tmp_ch2 = reshape(zscore(tmp_ch2(:)),size(tmp_ch2,1),size(tmp_ch2,2));
    subregion_ch2_norm(:,:,i) = tmp_ch2;
end