% Script for merging same axial position
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, August 2019
%       
function [ims_Zcal_merge,Zpos] = merge_Z_ave_img(ims_Zcal_ave,Zpos_in)


Numimage = size(ims_Zcal_ave,3);  

%% Re-arrange Z position 
num = 1;
Zpos = [];
while num <= Numimage
    if num <= Numimage -1 && abs(Zpos_in(num)-Zpos_in(num+1)) < 0.04
       Zpos(num) = (Zpos_in(num)+Zpos_in(num+1)) / 2;
       Zpos(num+1) = Zpos(num);
       num = num + 1;
    else 
       Zpos(num) = Zpos_in(num);
    end
    num = num + 1;
end


%% Merge same Z group
ims_Zcal_merge = ims_Zcal_ave;

num = 1;
while num < length(Zpos)    
    if Zpos(num) == Zpos(num+1)
        Zpos(num+1) = [];
        ims_Zcal_merge(:,:,num) = (ims_Zcal_merge(:,:,num) + ims_Zcal_merge(:,:,num+1)) /2; 
        ims_Zcal_merge(:,:,num+1) = [];
    end
    num = num + 1;
end


