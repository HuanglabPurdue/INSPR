% Script for pre-generating channel specific model
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, October 2019
%%
function [samplepsf,startx,starty,startz, dz,dx] = gensamplepsf(psfobj,pixelsize,psfsize,boxsize,bin,Nzs)
zpos = linspace(-1.3,1.3,Nzs)';
xpos = zeros(Nzs,1);
ypos = zeros(Nzs,1);
I = ones(Nzs,1);
bg = zeros(Nzs,1);
x1 = cat(2,xpos,ypos,zpos,I,bg);
w = [1,1,1,1,1];

psfobj.Pixelsize = pixelsize/bin;
psfobj.PSFsize = psfsize;
psfobj.Boxsize = boxsize*bin;

[psf2d_fit] = genpsf_real(x1,w,psfobj);

samplepsf = [];
f0 = psf2d_fit;
N = size(f0,1);
Nz = size(f0,3);
F = reshape(permute(f0,[2,1,3]),N*N*Nz,1);
samplepsf = F;

startx = -0.5*psfobj.Pixelsize*psfobj.Boxsize;
starty = -0.5*psfobj.Pixelsize*psfobj.Boxsize;
startz = zpos(1);
dz = zpos(2)-zpos(1);
dx = psfobj.Pixelsize;
end
