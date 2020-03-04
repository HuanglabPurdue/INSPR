% Script for pre-generating channel specific model
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, August 2019
%%
function [samplepsf,startx,starty,startz,dz,dx] = gensamplepsf_biplane(psfobj_all,pixelsize,psfsize,boxsize,bin,Nzs,dist_biplane,tform)
zpos = linspace(-1.3,1.3,Nzs)';
xpos = zeros(Nzs,1);
ypos = zeros(Nzs,1);
I = ones(Nzs,2);
bg = zeros(Nzs,2);
x1 = cat(2,xpos,ypos,zpos,I,bg);
w = [1,1,1,1,1];


for ii = 1:2
    psfobj_all{ii}.Pixelsize = pixelsize/bin;
    psfobj_all{ii}.PSFsize = psfsize;
    psfobj_all{ii}.Boxsize = boxsize*bin;
end

[psf2d_fit] = genpsf_biplane_real(x1,w,psfobj_all,dist_biplane);

edge_1 = mean(squeeze(psf2d_fit(1,:,:,2)),1)';
edge_2 = mean(squeeze(psf2d_fit(:,1,:,2)),1)';
edge_3 = mean(squeeze(psf2d_fit(end,:,:,2)),1)';
edge_4 = mean(squeeze(psf2d_fit(:,end,:,2)),1)';
edge_fill = min([edge_1 edge_2 edge_3 edge_4],[],2);

psf2d_fit(:,:,:,2) = imwarp(psf2d_fit(:,:,:,2),tform,'cubic','OutputView',imref2d(size(psf2d_fit(:,:,:,2))),...
    'FillValues',edge_fill);

psf_fit = [];
for ii = 1:2
    psf_fit = cat(2,psf_fit,squeeze(psf2d_fit(:,:,:,ii)));
end
samplepsf = cell(2,1);
for ii = 1:2
    f0 = psf2d_fit(:,:,:,ii);
    N = size(f0,1);
    Nz = size(f0,3);
    F = reshape(permute(f0,[2,1,3]),N*N*Nz,1);

    samplepsf{ii} = F;
end
startx = -0.5*psfobj_all{1}.Pixelsize*psfobj_all{1}.Boxsize;
starty = -0.5*psfobj_all{1}.Pixelsize*psfobj_all{1}.Boxsize;
startz = zpos(1);
dz = zpos(2)-zpos(1);
dx = psfobj_all{1}.Pixelsize;
end
