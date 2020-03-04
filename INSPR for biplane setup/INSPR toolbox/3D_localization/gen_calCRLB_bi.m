% Script for generating CRLB
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, December 2019
function crlbM = gen_calCRLB_bi(psfobj_all,PM,planedis,tform_fitting,offset_seg,bxsz)


%% biplane

oprobj = psfobj_all{1};

PRstruct = oprobj.PRstruct;
pxsz = oprobj.Pixelsize;
R = 128;
nmed = oprobj.PRstruct.RefractiveIndex;

crobj = CalCRLB_bi(PRstruct,PRstruct);
crobj.Pixelsize = pxsz; %micron
crobj.Xpos = PM(:,1);
crobj.Ypos = PM(:,2);
crobj.Zpos = PM(:,3) - planedis/2;
crobj.Photon = PM(:,4:5);
crobj.Bg = PM(:,6:7);
crobj.Boxsize = bxsz;   %subregion size
crobj.Deltax = 0.1; %pixel
crobj.Deltaz = 0.01;    %micron
crobj.PSFobj1.PSFsize = R;
crobj.PSFobj1.nMed = nmed;
crobj.PSFobj2.PSFsize = R;
crobj.PSFobj2.nMed = nmed;

crobj.Planedis = planedis;

crobj.prepInputparam();
crobj.calcrlb(tform_fitting,offset_seg);
% crobj.genfigs();

crlbM = crobj.CRLB;
crlbM(:,1:2) = crlbM(:,1:2)*pxsz*pxsz;

