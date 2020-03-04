%%
% Main script for running INSPR
%
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, December 2019
%
%% set environment paths

support_path = '..\Support'; %check support path

addpath([support_path '\PSF Toolbox']);
addpath([support_path '\SRsCMOS']);
addpath([support_path '\Helpers']);


%% call EMpupil GUI

INSPR_GUI();