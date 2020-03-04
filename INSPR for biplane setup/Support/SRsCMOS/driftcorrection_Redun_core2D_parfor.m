
function [R_shift1 R_shift2]=driftcorrection_Redun_core2D_parfor(x1in,x2in,toutin,frmnum,pixelsz,errothresh,cutmeth,currcy)


if nargin<6
    pixelsz=15;
    errothresh=7.5; % in nm
end

rndsz=2;
thresh=errothresh/pixelsz;

x1=floor((x1in-min(x1in))./pixelsz);
x2=floor((x2in-min(x2in))./pixelsz);

x1sz=max(24,ceil(max(x1)/rndsz)*rndsz+1);
x2sz=max(24,ceil(max(x2)/rndsz)*rndsz+1);
display(['drift_correction_reconstruction size: ' num2str(x1sz) ', ' num2str(x2sz)]);

if nargin<=6
    bomask1=(x1>=max(x1)*0.2)&(x1<=max(x1)*0.8);
    bomask2=(x2>=max(x2)*0.2)&(x2<=max(x2)*0.8);
elseif strcmp(cutmeth,'percentile')
    p1_20 = prctile(x1,20);
    p1_80 = prctile(x1,80);
    p2_20 = prctile(x2,20);
    p2_80 = prctile(x2,80);
    bomask1=(x1>=p1_20)&(x1<=p1_80);
    bomask2=(x2>=p2_20)&(x2<=p2_80);
elseif strcmp(cutmeth,'nocut')
    bomask1=x1>=1e-37;
    bomask2=x2>=1e-37;
end

% determine segments
segn=unique(currcy);

segnum = length(segn);
shift1_tmp = zeros(segnum, segnum);
shift2_tmp = zeros(segnum, segnum);

parfor ii=1:(segnum-1)

    disp(['Processing drift correction: '  num2str(ii)]);
    
    maskt=(currcy==segn(ii));
    x1_1=x1(maskt&bomask1&bomask2);
    x2_1=x2(maskt&bomask1&bomask2);
    im1=cHistRecon(x1sz,x2sz,x1_1,x2_1,0);

    for jj=1:segnum

        if jj >= ii+1 && jj <=min(ii+1+10,segnum)
            maskt=(currcy==segn(jj));
            x1_2=x1(maskt&bomask1&bomask2);
            x2_2=x2(maskt&bomask1&bomask2);
            
            im2=cHistRecon(x1sz,x2sz,x1_2,x2_2,0);
            [shift1_tmp(ii,jj) shift2_tmp(ii,jj)]=driftcorrection_core2D(imgaussfilt(single(im1),1),imgaussfilt(single(im2),1),[0 0]);
        end
    end
end


kk=1;
for ii = 1:1:segnum-1
    for jj=ii+1:1:min(ii+1+4,segnum)
        
        shift1(kk) = shift1_tmp(ii,jj);
        shift2(kk) = shift2_tmp(ii,jj);
        
        AA(kk,ii:jj-1)=1;
        kk=kk+1;
    end
end


R_shift1=pinv(AA)*shift1';
R_shift2=pinv(AA)*shift2';

error1=sqrt(((AA*R_shift1)-shift1').^2);
error2=sqrt(((AA*R_shift2)-shift2').^2);

newAA=AA;
oldAA=newAA;
oldshift=shift1;
while (rank(newAA)>=min(size(AA)))&&(max(error1)>thresh)
    oldAA=newAA;
    oldshift=shift1;
    [maxval ind]=max(error1);
    newAA(ind,:)=[];
    error1(ind)=[];
    shift1(ind)=[];
end

if rank(newAA)<min(size(AA))
    newAA=oldAA;
    shift1=oldshift;
end

R_shift1=pinv(newAA)*shift1';


newAA=AA;
oldAA=newAA;
oldshift=shift2;
while (rank(newAA)>=min(size(AA)))&&(max(error2)>thresh)
    oldAA=newAA;
    oldshift=shift2;
    [maxval ind]=max(error2);
    newAA(ind,:)=[];
    error2(ind)=[];
    shift2(ind)=[];
end

if rank(newAA)<min(size(AA))
    newAA=oldAA;
    shift2=oldshift;
end

R_shift2=pinv(newAA)*shift2';

R_shift1=R_shift1*pixelsz;
R_shift2=R_shift2*pixelsz;