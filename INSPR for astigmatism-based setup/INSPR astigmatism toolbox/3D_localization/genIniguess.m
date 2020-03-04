% Script for estimating initial parameters
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, October 2019
%%  
function x0 = genIniguess(data,bgtype,p,s1,s2,zstage)
sz = size(data);
if numel(sz) < 3
    Num=1;
else
    Num = sz(3);
end
% z
if nargin>2
    [P,CRLB,LL] = GPUgaussMLE(single(data),1,100,4);
    sx = P(:,5);
    sy = P(:,6);
    sr = (sx.^2-sy.^2)./(sx+sy);
    if ~isempty(p)
        z = polyval(p,sr);
    else
        p1 = polyval(s1,zstage);
        p2 = polyval(s2,zstage);
        z = polyval([p1,p2],sr);
    end
else
    z = zeros(Num,1);
end
% xy
start = floor(sz(1)/2);
subroi = data(start-2:start+2,start-2:start+2,:);
[x,y] = meshgrid([start-2:start+2],[start-2:start+2]);
xL = repmat(x,[1,1,Num]);
yL = repmat(y,[1,1,Num]);
area = squeeze(sum(sum(subroi,1),2));
comx = squeeze(sum(sum(subroi.*xL,1),2))./area;
comy = squeeze(sum(sum(subroi.*yL,1),2))./area;
% I,bg
switch bgtype
    case 'mean'
        subroi = data(2:end-1,2:end-1,:);
        sz1 = size(subroi);
        bg = (sum(sum(data,1),2)-sum(sum(subroi,1),2))./(sz(1)^2-sz1(1)^2);
    case 'median'
        tmp = data(:,:,1);
        tmp(2:end-1,2:end-1) = 0;
        mask = tmp ~= 0;
        psfxz = reshape(data,sz(1)*sz(1),sz(3));
        bg = ones(1,1,sz(3));
        bg(1,1,:) = median(psfxz(mask(:),:),1)';
end
BG = repmat(bg,[sz(1),sz(2),1]);
I = squeeze(sum(sum(data-BG,1),2));
I(I<200) = 200;
bg = squeeze(bg);
x0 = cat(2,comx,comy,I,bg,ones(Num,1).*z);
end