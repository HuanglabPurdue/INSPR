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
function [PM,crlbM,errM] = loc_ast_model_CPU(subregion_ch1,probj,setup)


%%
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
%% initial parameters

tic

data = subregion_ch1;

x1 = genIniguess(data,'median'); 
zLim = [-0.8,0.8];  
Nz = 81;
psftype = 'normal';
psfobj_input = psfobj;

[z_guess,mcc_val,psf_model] = genini_z_mat_parfor(data,psfobj_input,psftype,zLim,Nz);  %FX changed
x_next = x1(:,1) - bxsz/2;
y_next = x1(:,2) - bxsz/2;
z_next = z_guess;

I_next = x1(:,3);
bg_next = x1(:,4);

toc


%% localization
tic

data_tmp = data;
All_fit = size(data_tmp, 3);
interval = 1000;
N_loop = ceil(All_fit/interval);

disp(['The sub-regions were divided into ' num2str(N_loop) ' groups. Each group has around ' num2str(interval) ' sub-regions']);
numcores = feature('numcores');
disp(['The esimated time is around ' num2str(round(All_fit/1000*numcores*0.3)) ' min in this cycle']);

M_par = cell(N_loop,1);
parfor ii = 1 : N_loop
    disp(['Localization group: ' num2str(ii)]);

    fitobj = CalDev(PRstruct,'zernike');
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
        fitobj.Photon = Itmp;
        fitobj.Bg = bgtmp;
        
        fitobj.prepInputparam();
        fitobj.caldev();
        step = fitobj.getnextstep(data_selection);
       
        xtmp = xtmp + min([max([step(:,1),-xa.*ones(size(xtmp))],[],2),xa.*ones(size(xtmp))],[],2);
        ytmp = ytmp + min([max([step(:,2),-xa.*ones(size(ytmp))],[],2),xa.*ones(size(ytmp))],[],2);
        ztmp = ztmp + min([max([step(:,3),-za.*ones(size(ztmp))],[],2),za.*ones(size(ztmp))],[],2);

        Itmp = Itmp + min([max([step(:,4),-Ia.*ones(size(Itmp))],[],2),Ia.*ones(size(Itmp))],[],2);
        bgtmp = bgtmp + min([max([step(:,5),-bga.*ones(size(bgtmp))],[],2),bga.*ones(size(bgtmp))],[],2);
           
        xtmp = min(max(xtmp,-4),4);
        ytmp = min(max(ytmp,-4),4);
        ztmp = min(max(ztmp,-1.2),1.2);
        Itmp(Itmp<=100) = 100;
        bgtmp(bgtmp<=0) = 1e-3;
    end
    Ptmp = cat(2,xtmp,ytmp,ztmp,Itmp,bgtmp);    
    M_par{ii}.Ptmp = Ptmp;
    
    %% generate CRLB
    crlb_tmp = gen_calCRLB(psfobj_input,Ptmp,bxsz);
    M_par{ii}.crlb_tmp = crlb_tmp;
    
    %% generate LLR    
    err_tmp = gen_LLR(data_selection,psfobj_input,Ptmp);
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
                    

