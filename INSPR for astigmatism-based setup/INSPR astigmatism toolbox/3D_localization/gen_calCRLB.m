% Script for generating CRLB
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, December 2019
%%
function crlbM = gen_calCRLB(psfobj_input,PM,bxsz)

oprobj = psfobj_input;

PRstruct = oprobj.PRstruct;
pxsz = oprobj.Pixelsize;
R = 128;
nmed = oprobj.PRstruct.RefractiveIndex;

crobj = CalCRLB(PRstruct,'zernike');
crobj.Pixelsize = pxsz; %micron
crobj.Xpos = PM(:,1);
crobj.Ypos = PM(:,2);
crobj.Zpos = PM(:,3);
crobj.Photon = PM(:,4);
crobj.Bg = PM(:,5);
crobj.Boxsize = bxsz;   %subregion size
crobj.Deltax = 0.1; %pixel
crobj.Deltaz = 0.01;    %micron
crobj.PSFobj.PSFsize = R;
crobj.PSFobj.nMed = nmed;


crobj.prepInputparam();
crobj.calcrlb();
% crobj.genfigs();

crlbM = crobj.CRLB;
crlbM(:,1:2) = crlbM(:,1:2)*pxsz*pxsz;

