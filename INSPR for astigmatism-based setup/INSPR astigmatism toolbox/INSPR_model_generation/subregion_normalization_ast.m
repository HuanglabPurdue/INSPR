% Script for normalizing sub-regions
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, October 2019
% 
function subregion_ch1_norm = subregion_normalization_ast(subregion_ch1)

subregion_ch1_norm = [];

for i = 1 : size(subregion_ch1,3)
    tmp_ch1 = subregion_ch1(:,:,i);
    tmp_ch1 = reshape(zscore(tmp_ch1(:)),size(tmp_ch1,1),size(tmp_ch1,2));  %Standardized z-scores
    subregion_ch1_norm(:,:,i) = tmp_ch1;
end