% Script for pupil-based 3D localization: CPU version
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, December 2019

%%
function [PM,crlbM,errM] = loc_channel_specific_model_CPU(subregion_ch1,subregion_ch2,probj,probj2,dist_biplane,tform,offset_seg,setup)

disp('Use CPU version');

oprobj = probj;
PRstruct = oprobj.PRstruct;
bxsz = size(subregion_ch1,1);
pxsz = oprobj.Pixelsize;
R = 128;
psfobj = PSF_zernike(PRstruct);
psfobj.Boxsize = bxsz;
psfobj.Pixelsize = pxsz; % micron
psfobj.PSFsize = R;
psfobj.nMed = oprobj.PRstruct.RefractiveIndex;

psfobj.precomputeParam();  
psfobj.setPupil();
if setup.is_imgsz == 1    
    psfobj.genZernikeMag();
end
%%
oprobj2 = probj2;
PRstruct2 = oprobj2.PRstruct;
bxsz = size(subregion_ch2,1);
pxsz = oprobj2.Pixelsize;
R = 128;
psfobj2 = PSF_zernike(PRstruct2);
psfobj2.Boxsize = bxsz;
psfobj2.Pixelsize = pxsz; % micron
psfobj2.PSFsize = R;
psfobj2.nMed = oprobj2.PRstruct.RefractiveIndex;

psfobj2.precomputeParam();  
psfobj2.setPupil();
if setup.is_imgsz == 1    
    psfobj2.genZernikeMag();
end

%% initial parameters
tic

data = cat(4,subregion_ch1,subregion_ch2);

x1 = genIniguess(mean(data,4),'median'); 
zLim = [-0.9,0.9];  
Nz = 81;
psftype = 'normal';
psfobj_all{1} = psfobj;
psfobj_all{2} = psfobj2;
planedis = [-dist_biplane/2,dist_biplane/2]; 

ref_point = [8, 8];   
[tform_fitting, ~] = cal_model_affine(tform,ref_point);
[z_guess,mcc_val,psf_model] = geniniBiplane_z_mat_parfor(data,psfobj_all,psftype,zLim,Nz,planedis,tform_fitting);  
x_next = x1(:,1) - bxsz/2;
y_next = x1(:,2) - bxsz/2;
z_next = z_guess;

x_1 = genIniguess(data(:,:,:,1),'median'); 
x_2 = genIniguess(data(:,:,:,2),'median'); 

I_next = cat(2,x_1(:,3),x_2(:,3));
bg_next = cat(2,x_1(:,4),x_2(:,4));

toc

%% localization
tic

data_tmp = data;
All_fit = size(data_tmp, 3);
interval = 1000;
N_loop = ceil(All_fit/interval);

disp(['The sub-regions were divided into ' num2str(N_loop) ' groups. Each group has around ' num2str(interval) ' sub-regions']);
numcores = feature('numcores');
disp(['The esimated time is around ' num2str(round(All_fit/1000*numcores*0.6)) ' min in this cycle']);


M_par = cell(N_loop,1);
parfor ii = 1 : N_loop
    disp(['Localization group: ' num2str(ii)]);

    fitobj = CalDevBi(PRstruct,'zernike');
    fitobj.Boxsize = bxsz;
    fitobj.Pixelsize = pxsz;
    fitobj.Deltax = 0.1;
    fitobj.Deltaz = 0.01;
    fitobj.PSFobj = psfobj;

    index_start = (ii-1) * interval + 1;
    if ii == N_loop
        index_end = All_fit;
    else
        index_end = ii * interval;
    end

    % input data
    data_selection = data_tmp(:,:,index_start:index_end,:);
    Nfit = size(data_selection,3);    
    offset_seg_tmp = offset_seg(index_start:index_end,:);
    xtmp = x_next(index_start:index_end);
    ytmp = y_next(index_start:index_end);
    ztmp = z_next(index_start:index_end);
    Itmp = I_next(index_start:index_end,:);
    bgtmp = bg_next(index_start:index_end,:);
    
    iterateN = 30;
     
    xa = 0.2;   %moving step size
    za = 0.06;  
    Ia = 400;
    bga = 5;
    %% estimate position
    for jj = 1:iterateN     
        fitobj.Xpos = xtmp;   %in this case, (0,0) is the center
        fitobj.Ypos = ytmp;
        fitobj.Zpos = ztmp;
        step = fitobj.getnextstep(psfobj_all,data_selection,Itmp,bgtmp,planedis,tform_fitting,offset_seg_tmp);
        xtmp = xtmp + min([max([step.x,-xa.*ones(size(xtmp))],[],2),xa.*ones(size(xtmp))],[],2);
        ytmp = ytmp + min([max([step.y,-xa.*ones(size(ytmp))],[],2),xa.*ones(size(ytmp))],[],2);
        ztmp = ztmp + min([max([step.z,-za.*ones(size(ztmp))],[],2),za.*ones(size(ztmp))],[],2);

        Itmp = Itmp + min(cat(3,max(cat(3,step.I,-Ia.*ones(size(Itmp))),[],3),Ia.*ones(size(Itmp))),[],3);
        bgtmp = bgtmp + min(cat(3,max(cat(3,step.bg,-bga.*ones(size(bgtmp))),[],3),bga.*ones(size(bgtmp))),[],3);
           
        xtmp = min(max(xtmp,-4),4);
        ytmp = min(max(ytmp,-4),4);
        ztmp = min(max(ztmp,-1.2),1.2);
        Itmp(Itmp<=100) = 100;
        bgtmp(bgtmp<=0) = 1e-3;
        
    end
    Ptmp = cat(2,xtmp,ytmp,ztmp,Itmp,bgtmp);    
    M_par{ii}.Ptmp = Ptmp;
    
    %% generate CRLB
    crlb_tmp = gen_calCRLB_bi(psfobj_all,Ptmp,dist_biplane,tform_fitting,offset_seg_tmp,bxsz);
    M_par{ii}.crlb_tmp = crlb_tmp;
    
    %% generate LLR    
    err_tmp = gen_LLR_bi(data_selection,psfobj_all,Ptmp,planedis,tform_fitting,offset_seg_tmp);
    M_par{ii}.err_tmp = err_tmp;
end

PM = [];
crlbM = [];
errM = [];
for ii = 1 : N_loop
    PM = cat(1,PM,M_par{ii}.Ptmp);
    crlbM = cat(1,crlbM,M_par{ii}.crlb_tmp);
    errM = cat(1,errM,M_par{ii}.err_tmp);
end

PM(:,1) = PM(:,1) + (bxsz+1)/2;
PM(:,2) = PM(:,2) + (bxsz+1)/2;

toc
