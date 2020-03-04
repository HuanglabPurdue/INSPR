% Script for pupil-based 3D localization
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, August 2019
%%
function [PcudaM,crlbM,errM] = loc_channel_specific_model(subregion_ch1,subregion_ch2,probj,probj2,dist_biplane,tform,offset_seg,subvar_ch1,subvar_ch2,setup)

disp('Using GPU version');
%%
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
[tform_fitting, tMatrix_noTranslation] = cal_model_affine(tform,ref_point);
[z_guess,mcc_val,psf_model] = geniniBiplane_z_mat_parfor(data,psfobj_all,psftype,zLim,Nz,planedis,tform_fitting);  
x_next = x1(:,1)-bxsz/2;
y_next = x1(:,2)-bxsz/2;
z_next = z_guess;

x_1 = genIniguess(data(:,:,:,1),'median'); 
x_2 = genIniguess(data(:,:,:,2),'median'); 

I_next = cat(2,x_1(:,3),x_2(:,3));
bg_next = cat(2,x_1(:,4),x_2(:,4));

toc

%% generate sample PSF
tic

pixelsize = pxsz;
boxsize = 25;

bin = 4;
psfsize = 128 * bin;
Nzs = 701;
ref_point = [49, 49];  
[tform_fitting, ~] = cal_model_affine(tform,ref_point);

[samplepsf,startx,starty,startz,dz,dx] = gensamplepsf_biplane(psfobj_all,pixelsize,psfsize,boxsize,bin,Nzs,dist_biplane,tform_fitting);

samplepsf_4d = [];
samplepsf_cuda = [];
for ii = 1:2
    samplepsf_cuda = cat(3,samplepsf_cuda,permute(reshape(samplepsf{ii},boxsize*bin,boxsize*bin,Nzs),[2,1,3,4]));
    samplepsf_4d = cat(4,samplepsf_4d,permute(reshape(samplepsf{ii},boxsize*bin,boxsize*bin,Nzs),[2,1,3,4]));  %For calculating spline parameters     
end
samplepsf_cuda = single(samplepsf_cuda);

toc      

%% Calculate spline parameters

st1 = genpsfstruct(samplepsf_4d(:,:,:,1),dx,dz,'matrix');
st2 = genpsfstruct(samplepsf_4d(:,:,:,2),dx,dz,'matrix');
st = catstruct(st1,st2,3);

%% localization
tic

data_tmp = data;
var_tmp = cat(4,subvar_ch1,subvar_ch2);
boxsize_data = 16;

All_fit = size(data_tmp, 3);
interval = 50000;
N_loop = ceil(All_fit/interval);

PcudaM = [];
crlbM = [];
errM = [];
for ii = 1 : N_loop
    index_start = (ii-1) * interval + 1;
    if ii == N_loop
        index_end = All_fit;
    else
        index_end = ii * interval;
    end

    % input data
    data_selection = data_tmp(:,:,index_start:index_end,:);
    Nfit = size(data_selection,3);
    data_cuda = single(data_selection(:)); %psfsize x Nfit x 2
    
    % input sCMOS noise parameter
    var_selection = var_tmp(:,:,index_start:index_end,:);
    var_input = single(var_selection(:)); %psfsize x Nfit x 2

    
    coords_cuda = cat(2,ones(Nfit,2),zeros(Nfit,1))';
    coords_cuda = single(coords_cuda(:));
    
    offset_seg_tmp = offset_seg(index_start:index_end,:);
    xtmp = x_next(index_start:index_end) + boxsize_data/2;
    ytmp = y_next(index_start:index_end) + boxsize_data/2;
    x0 = cat(2,xtmp,ytmp,z_next(index_start:index_end),I_next(index_start:index_end,:), bg_next(index_start:index_end,:));
    x0i = single(reshape(x0',Nfit*7,1));
    lambda = 0;
    
    iterateN = 150;
                                
                                
    [P_cuda,CG,crlb,err,PSF_cuda] = cuda_channel_specific_model(data_cuda,coords_cuda,samplepsf_cuda,...
                                    dx,dz,startx,starty,startz,iterateN,Nfit,lambda,x0i,...
                                    single(st.Fx),single(st.Fy),single(st.Fz),...
                                    single(st.Fxy),single(st.Fxz),single(st.Fyz),single(st.Fxyz),...
                                    single(tMatrix_noTranslation),single(offset_seg_tmp),single(var_input));
                                
                                
                                
    
    PcudaM = [PcudaM; reshape(P_cuda,7,Nfit)'];
    crlbM = [crlbM; reshape(crlb,7,Nfit)']; 
    errM = [errM; reshape(err,2,Nfit)'];

end

toc                         
clear cuda_biplane_sCMOS_loc_spline_affine                         

                    
