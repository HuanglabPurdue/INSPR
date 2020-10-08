%% 2D drift correction method, detail see motion correction in cryo-EM
% INPUT
%    x1/2in: input x y postions in pixels 
%    pixelsz: pixel size of orginal image
%    errothresh: error threshold
%    cormask: section index. Single molecules of some sectiion will assign together 
%      N x 1 integer array grouping observations for corellation (eg filenumber)
%    interval: maximum image interval in analyzing drift 
% OUTPUT
%    xyshift: Found shift * pixel size
%    error_final: error parameters
%
% Edited by Fan Xu, 2018-02-05



function [xshift, yshift, error1_final, error2_final]=driftcorrection_Redun2D(x1in,x2in,pixelsz,errothresh,cormask,interval)

if nargin < 6
    interval = 5;
end

if numel(unique(cormask)) == 1
    xshift = 0;
    yshift = 0;
    return
end

%threshold
thresh = errothresh / pixelsz;

x1 = floor((x1in - min(x1in(:))));
x2 = floor((x2in - min(x2in(:))));

rndsz = 2;
x1sz = ceil(max(x1)/rndsz)*rndsz;
x2sz = ceil(max(x2)/rndsz)*rndsz;
sz = [x1sz x2sz];

%zoom in 2
cczoom = 2;
x1sz = (sz(1)*cczoom);
x2sz = (sz(2)*cczoom);
x1 = single(cczoom*(x1in - min(x1in)));
x2 = single(cczoom*(x2in - min(x2in)));

% max section index
segnum = max(cormask);
shift1_tmp = zeros(segnum, segnum);
shift2_tmp = zeros(segnum, segnum);

%% find drift between two image cycle
parfor ii = 1:1:segnum-1
    disp(['Processing drift correction: '  num2str(ii)]);
     
    x1_1 = x1(cormask == ii);
    x2_1 = x2(cormask == ii);
    
    if size(x1_1) < 1000
        disp(['Error in cycle of drift corection: ' num2str(ii) ' so less localization spots, please check the image']);
    else
        im1 = cHistRecon(x1sz, x2sz, single(x1_1), single(x2_1), 1);
        
        %analysis all image, or set min(interval,segnum)
        for jj=1:segnum
            if jj >= ii+1 && jj <=min(ii+1+interval,segnum)
                
                x1_2 = x1(cormask == jj);
                x2_2 = x2(cormask == jj);
                
                if size(x1_2) < 1000
                    disp(['Error in cycle of drift corection: ' num2str(jj) ' so less localization spots, please check the image']);
                else
                    im2 = cHistRecon(x1sz, x2sz, single(x1_2), single(x2_2),1);
                    
                    %calcualte drift between two images
                    [shift1_tmp(ii,jj) shift2_tmp(ii,jj)]=driftcorrection_core2D(imgaussfilt(single(im1),1),imgaussfilt(single(im2),1),[0 0]);
  
                end
            end
        end
    end
end


kk=1;
for ii = 1:1:segnum-1
    for jj=ii+1:1:min(ii+1+interval,segnum)
        
        shift1(kk) = shift1_tmp(ii,jj);
        shift2(kk) = shift2_tmp(ii,jj);
        
        AA(kk,ii:jj-1)=1;
        kk=kk+1;
    end
end

%% calcate R_shift according to redundancy-based method
R_shift1 = pinv(AA)*shift1';
R_shift2 = pinv(AA)*shift2';

error1 = sqrt(((AA*R_shift1)-shift1').^2);
error2 = sqrt(((AA*R_shift2)-shift2').^2);

newAA = AA;
oldAA = newAA;
oldshift = shift1;
while rank(newAA)>=min(size(AA)) && max(error1)>thresh
    oldAA = newAA;
    oldshift = shift1;
    [maxval ind] = max(error1);
    newAA(ind,:) = [];
    error1(ind) = [];
    shift1(ind) = [];
end

if rank(newAA) < min(size(AA))
    newAA = oldAA;
    shift1 = oldshift;
    flag(1) = 1;
end

R_shift1 = pinv(newAA)*shift1';
error1_final = sqrt(((newAA*R_shift1)-shift1').^2);

newAA = AA;
oldAA = newAA;
oldshift = shift2;
while (rank(newAA)>=min(size(AA))) && (max(error2)>thresh)
    oldAA = newAA;
    oldshift = shift2;
    [maxval ind] = max(error2);
    newAA(ind,:) = [];
    error2(ind) = [];
    shift2(ind) = [];
end

if rank(newAA) < min(size(AA))
    newAA = oldAA;
    shift2 = oldshift;
    flag(2) = 1;
end

error2_final = sqrt(((newAA*R_shift2)-shift2').^2);
R_shift2 = pinv(newAA)*shift2';


%% output XY shift
xshift = zeros(segnum,1);
yshift = zeros(segnum,1);

xshift(1) = 0;
yshift(1) = 0;
for ii = 1 : segnum-1
    xshift(ii+1) = xshift(ii) - R_shift1(ii);
    yshift(ii+1) = yshift(ii) - R_shift2(ii);
end

xshift =(xshift/cczoom) * pixelsz;
yshift =(yshift/cczoom) * pixelsz; 

