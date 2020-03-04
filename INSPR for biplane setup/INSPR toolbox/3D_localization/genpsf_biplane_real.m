% Script for generating model from pupil
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, August 2019
%%
function [psf_biplane] = genpsf_biplane_real(x,w,psfobj_all,dist_biplane)


R = psfobj_all{1}.Boxsize;
Zpos_change = [-dist_biplane/2 dist_biplane/2];
N = numel(x(:,1));

psf_biplane = zeros(R,R,N,2);
for nn = 1:2
    
    psfobj = psfobj_all{nn};
    psfobj.Xpos = x(:,1).*w(1);
    psfobj.Ypos = x(:,2).*w(2);
    psfobj.Zpos = x(:,3).*w(3) + Zpos_change(nn);
    
    psfobj.precomputeParam();
    psfobj.padPupil();
    psfobj.genPSF_2();


    psfobj.scalePSF('normal');
    

    
    psfI = psfobj.ScaledPSFs;
    normf = sum(sum(psfobj.Pupil.mag.^2,1),2);
    psf = psfI./normf;
    
    I = x(:,nn+3);
    tmp = zeros(1,1,N);
    tmp(1,1,:) = I;
    IL = repmat(tmp,[R,R,1]).*w(4);

    bg = x(:,nn+5);
    tmp = zeros(1,1,N);
    tmp(1,1,:) = bg;
    bgL = repmat(tmp,[R,R,1]).*w(5);
    
    psf_biplane(:,:,:,nn) = psf.*IL+bgL;
end

end
