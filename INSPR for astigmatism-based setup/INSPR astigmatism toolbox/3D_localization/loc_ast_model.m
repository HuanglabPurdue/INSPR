% Script for pupil-based 3D localization
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, October 2019
%%
function [PcudaM,crlbM,errM] = loc_ast_model(subregion_ch1,probj,subvar_ch1,setup)

disp('Use GPU version');
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

%% initial parameters

tic

data = subregion_ch1;

x1 = genIniguess(data,'median'); 
zLim = [-0.8,0.8];  
Nz = 81;
psftype = 'normal';
psfobj_input = psfobj;

[z_guess,mcc_val,psf_model] = genini_z_mat_parfor(data,psfobj_input,psftype,zLim,Nz);  %FX changed
x_next = x1(:,1)-bxsz/2;
y_next = x1(:,2)-bxsz/2;
z_next = z_guess;

I_next = x1(:,3);
bg_next = x1(:,4);

toc

%% generate sample PSF
tic

pixelsize = pxsz;
boxsize = 25;

bin = 4;
psfsize = 128 * bin;
Nzs = 701;
[samplepsf,startx,starty,startz,dz,dx] = gensamplepsf(psfobj,pixelsize,psfsize,boxsize,bin,Nzs);

samplepsf_cuda = permute(reshape(samplepsf,boxsize*bin,boxsize*bin,Nzs),[2,1,3]);
samplepsf_3d = permute(reshape(samplepsf,boxsize*bin,boxsize*bin,Nzs),[2,1,3]);  %For calculating spline parameters     

samplepsf_cuda = single(samplepsf_cuda);

toc      

%% Calculate spline parameters

st = genpsfstruct(samplepsf_3d,dx,dz,'matrix');

%% localization
tic

data_tmp = data;
var_tmp = subvar_ch1;
boxsize_data = 16;

All_fit = size(data_tmp, 3);
interval = 40000;
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
    data_cuda = single(data_selection(:)); %psfsize x Nfit
    
    % input sCMOS noise
    var_selection = var_tmp(:,:,index_start:index_end,:);
    var_input = single(var_selection(:)); %psfsize x Nfit
    
    coords_cuda = cat(2,ones(Nfit,2),zeros(Nfit,1))';
    coords_cuda = single(coords_cuda(:));
    
    xtmp = x_next(index_start:index_end) + boxsize_data/2;
    ytmp = y_next(index_start:index_end) + boxsize_data/2;
    x0 = cat(2,xtmp,ytmp,z_next(index_start:index_end),I_next(index_start:index_end), bg_next(index_start:index_end));
    x0i = single(reshape(x0',Nfit*5,1));
    lambda = 0;% damping factor
    
    iterateN = 150;

    [P_cuda,CG,crlb,err,PSF_cuda] = cuda_ast_model(data_cuda,coords_cuda,samplepsf_cuda,...
                                dx,dz,startx,starty,startz,iterateN,Nfit,lambda,x0i,...
                                single(st.Fx),single(st.Fy),single(st.Fz),...
                                single(st.Fxy),single(st.Fxz),single(st.Fyz),single(st.Fxyz),single(var_input));
    
  
    PcudaM = [PcudaM; reshape(P_cuda,5,Nfit)'];
    crlbM = [crlbM; reshape(crlb,5,Nfit)']; 
    errM = [errM; reshape(err,2,Nfit)'];

end

toc                         
clear cuda_ast_model                         

                    



