%%
% Main script for running INSPR
% This code is used for astigmatism based setup
%
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, July 2020
% INSPR 1.1: Add background subtraction using temporal median filter
%

%% set environment paths

support_path = '..\Support'; %check support path

addpath([support_path '\PSF Toolbox']);
addpath([support_path '\SRsCMOS']);
addpath([support_path '\Helpers']);


%% call EMpupil GUI

INSPR_ast_GUI();