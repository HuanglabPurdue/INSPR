%%
% Main script for running INSPR
%
% (C) Copyright 2020                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, December 2020
% INSPR 1.1: Add background subtraction using temporal median filter
% INSPR 1.3: Add the usage information for mouseover events, update GUI interface and related functions, and fix some bugs. 
% 
%% set environment paths

support_path = '..\Support'; %check support path

addpath([support_path '\PSF Toolbox']);
addpath([support_path '\SRsCMOS']);
addpath([support_path '\Helpers']);


%% call EMpupil GUI

INSPR_GUI();