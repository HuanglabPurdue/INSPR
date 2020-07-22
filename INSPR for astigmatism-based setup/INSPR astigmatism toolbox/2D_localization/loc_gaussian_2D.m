% Script for 2D Gaussian localization
% (C) Copyright 2020                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, July 2020
%%
function [sub_xM,sub_yM,sub_photonM,sub_bgM,crlbM,llM] = loc_gaussian_2D(subregion_ch1,subvar_ch1)

disp('Run 2D Gaussian fitting');                     

                    
%% parameters
PSFSigma = 1.1;
Iterations = 50;
FitType = 1;

%% localization
tic

All_fit = size(subregion_ch1, 3);
interval = 40000;
N_loop = ceil(All_fit/interval);

sub_xM = [];
sub_yM = [];
sub_photonM = [];
sub_bgM = [];
crlbM = [];
llM = [];
for ii = 1 : N_loop
    index_start = (ii-1) * interval + 1;
    if ii == N_loop
        index_end = All_fit;
    else
        index_end = ii * interval;
    end

    data_selection = subregion_ch1(:,:,index_start:index_end);
    var_selection = subvar_ch1(:,:,index_start:index_end);               
    
    
    [sub_x,sub_y,sub_photon,sub_bg,crlb,ll]=SRsCMOS_MLE(single(data_selection),single(var_selection),single(var_selection.*0+1),FitType,PSFSigma,Iterations);

%     clear mex
%     pause(1)
    
    sub_xM = [sub_xM; sub_x];
    sub_yM = [sub_yM; sub_y];
    sub_photonM = [sub_photonM; sub_photon];
    sub_bgM = [sub_bgM; sub_bg];
    crlbM = [crlbM; crlb]; 
    llM = [llM; -2*ll];
end
toc  


