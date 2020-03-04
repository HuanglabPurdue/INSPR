% Script for generating model from pupil
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, October 2019
%%
function [psf_oneplane] = genpsf_real(x,w,obj)
R = obj.Boxsize;
obj.Xpos = x(:,1).*w(1);
obj.Ypos = x(:,2).*w(2);

N = numel(x(:,1));

psf_oneplane = zeros(R,R,N);

obj.Zpos = x(:,3).*w(3);
obj.precomputeParam();
obj.padPupil();
obj.genPSF_2();

obj.scalePSF('normal');
psfI = obj.ScaledPSFs;

normf = sum(sum(obj.Pupil.mag.^2,1),2);
psf = psfI./normf;

I = x(:,4);
tmp = zeros(1,1,N);
tmp(1,1,:) = I;
IL = repmat(tmp,[R,R,1]).*w(4);

bg = x(:,5);
tmp = zeros(1,1,N);
tmp(1,1,:) = bg;
bgL = repmat(tmp,[R,R,1]).*w(5);

psf_oneplane(:,:,:) = psf.*IL+bgL;

end
