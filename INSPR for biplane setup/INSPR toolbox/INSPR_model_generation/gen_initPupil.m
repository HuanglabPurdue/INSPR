% Script for generating initial pupil
%
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, August 2019
%       
function probj = gen_initPupil(empupil)

disp('generate initial pupil');

%generate Z images
R = 128;
PRstruct = [];
PRstruct.NA = empupil.NA; 
PRstruct.Lambda = empupil.Lambda; 
PRstruct.RefractiveIndex = empupil.nMed ;    
PRstruct.Pupil.phase = zeros(R,R);
PRstruct.Pupil.mag = zeros(R,R);

phaseZ = zeros(1,empupil.Zernike_sz);
phaseZ([5:25]) = empupil.init_z; %initial Zernike

magZ = zeros(1,empupil.Zernike_sz);
magZ(1) = 1;
PRstruct.Zernike_phase = phaseZ;
PRstruct.Zernike_mag = magZ;
PRstruct.SigmaX = empupil.blur_sigma;
PRstruct.SigmaY = empupil.blur_sigma;

probj = PSF_zernike(PRstruct);
probj.Pixelsize = empupil.Pixelsize; % micron
probj.PSFsize = R; 
probj.nMed = empupil.nMed;   
probj.precomputeParam();
probj.genZernike();
probj.genPupil();
