function [R_shift1 R_shift2 R_shift3 flag error1_final error2_final error3_final]=stackalignment_Redun3D(stksort,x1in,x2in,x3in,pixelsz,errothresh,cutmeth,iniguess,maskflag,cropsz)

flag=[0 0];
if nargin<6
    pixelsz=15;
    errothresh=15; % in nm
end

if nargin<7
    inipx=[0 0];
else
    inipx=round(iniguess./pixelsz);
end

rndsz=2;
thresh=errothresh/pixelsz;

minx1=min(cat(1,x1in{:}));
minx2=min(cat(1,x2in{:}));
minx3=min(cat(1,x3in{:}));

x1=cellfun(@(x) (x-minx1)./pixelsz, x1in,'UniformOutput',false);
x2=cellfun(@(x) (x-minx2)./pixelsz, x2in,'UniformOutput',false);
x3=cellfun(@(x) (x-minx3)./pixelsz, x3in,'UniformOutput',false);


x1sz=ceil(max((cat(1,x1{:})))./rndsz)*rndsz+1;
x2sz=ceil(max((cat(1,x2{:})))./rndsz)*rndsz+1;
x3sz=ceil(max((cat(1,x3{:})))./rndsz)*rndsz+1;

display(['3D drift_correction_reconstruction size: ' num2str(x1sz) ', ' num2str(x2sz) ', ' num2str(x3sz)]);


if nargin<=5
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
%     bomask1=x1>=1e-37;
%     bomask2=x2>=1e-37;
end



segnum=numel(x1in);
kk=1;
for ii=1:1:segnum
%     ii

    
    x1_1=x1{ii};
    x2_1=x2{ii};
    x3_1=x3{ii};
    im1=cHistRecon3D(x1sz,x2sz,x3sz,single(x1_1),single(x2_1),single(x3_1),0);

    for jj=ii+1:1:min(ii+1,segnum)
 
        x1_2=x1{jj};
        x2_2=x2{jj};
        x3_2=x3{jj};
        
        im2=cHistRecon3D(x1sz,x2sz,x3sz,single(x1_2),single(x2_2),single(x3_2),0);
        
        [shift1(kk) shift2(kk) shift3(kk)]=stackalignment_core3D(single(im1),single(im2),inipx*(jj-ii));
        AA(kk,ii:jj-1)=1;
        kk=kk+1;
    end
end

R_shift1=pinv(AA)*shift1';
R_shift2=pinv(AA)*shift2';
R_shift3=pinv(AA)*shift3';

error1=sqrt(((AA*R_shift1)-shift1').^2);
error2=sqrt(((AA*R_shift2)-shift2').^2);
error3=sqrt(((AA*R_shift3)-shift3').^2);

newAA=AA;
oldAA=newAA;
oldshift=shift1;
while rank(newAA)>=min(size(AA))&&max(error1)>thresh
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
    flag(1)=1;
end

R_shift1=pinv(newAA)*shift1';
error1_final=sqrt(((newAA*R_shift1)-shift1').^2);

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
    flag(2)=1;
end

error2_final=sqrt(((newAA*R_shift2)-shift2').^2);
R_shift2=pinv(newAA)*shift2';

newAA=AA;
oldAA=newAA;
oldshift=shift3;
while rank(newAA)>=min(size(AA))&&max(error3)>thresh
    oldAA=newAA;
    oldshift=shift3;
    [maxval ind]=max(error3);
    newAA(ind,:)=[];
    error3(ind)=[];
    shift3(ind)=[];
end

if rank(newAA)<min(size(AA))
    newAA=oldAA;
    shift3=oldshift;
    flag(1)=1;
end

R_shift3=pinv(newAA)*shift3';
error3_final=sqrt(((newAA*R_shift3)-shift3').^2);

R_shift1=R_shift1*pixelsz;
R_shift2=R_shift2*pixelsz;
R_shift3=R_shift3*pixelsz;