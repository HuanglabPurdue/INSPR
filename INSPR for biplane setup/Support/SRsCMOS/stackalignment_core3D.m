function [shift1 shift2 shift3]=stackalignment_core3D(im1,im2,iniguess)

kernelsize=1;   
if isempty(im1)||isempty(im2)||(sum(im1(:))==0)||(sum(im2(:))==0)
    shift1=1e25;
    shift2=1e25;
    shift3=1e25;
else
    corrim3d=normxcorr3_sparse(single(imgaussfilt3(im1,kernelsize)),single(imgaussfilt3(im2,kernelsize)),'full');
    
    if nargin>=3
        shift_col=iniguess(2);
        shift_row=iniguess(1);
        shift_z=iniguess(3);
        rowval=shift_row+round((size(corrim3d,1)+1)./2);
        colval=shift_col+round((size(corrim3d,2)+1)./2);
        zval=shift_z+round((size(corrim3d,3)+1)./2);
        cropsz2=12;
        ori=[colval-cropsz2/2-1 rowval-cropsz2/2-1 zval-cropsz2/2-1];
%         smallim=cut(corrim3d,cropsz2,ori);
        xLeft = max(rowval-cropsz2/2,1);
        xRight = min(rowval+cropsz2/2-1, size(corrim3d,1));
        yLeft = max(colval-cropsz2/2,1);
        yRight = min(colval+cropsz2/2-1, size(corrim3d,2));
        zLeft = max(zval-cropsz2/2,1);
        zRight = min(zval+cropsz2/2-1, size(corrim3d,3));
        smallim = corrim3d(xLeft:xRight,yLeft:yRight,zLeft:zRight);

        [tmpval colval2 rowval2 zval2]=findmax3d(double(smallim));
        colval=ori(1)+colval2;
        rowval=ori(2)+rowval2;
        zval=ori(3)+zval2;
    else
        [tmpval colval rowval zval]=findmax3d(double(corrim3d));
    end
    
    shift_row=(rowval-(size(corrim3d,1)+1)./2);
    shift_col=(colval-(size(corrim3d,2)+1)./2);
    shift_z=(zval-(size(corrim3d,3)+1)./2);
    cropsz=12;
    exsz=256;
    ori=[colval-cropsz/2-1 rowval-cropsz/2-1 zval-cropsz/2-1];
    
    xLeft = max(rowval-cropsz/2,1);
    xRight = min(rowval+cropsz/2-1, size(corrim3d,1));
    yLeft = max(colval-cropsz/2,1);
    yRight = min(colval+cropsz/2-1, size(corrim3d,2));
    zLeft = max(zval-cropsz/2,1);
    zRight = min(zval+cropsz/2-1, size(corrim3d,3));

    rstmp_fft = fftshift(fftn(corrim3d(xLeft:xRight,yLeft:yRight,zLeft:zRight)));
    rstmp_pad = padarray(rstmp_fft, [(exsz-cropsz)/2 (exsz-cropsz)/2 (exsz-cropsz)/2]);
    rscorrim = ifftn(rstmp_pad);
%     rscorrim=ift(extend(ft(cut(corrim3d,cropsz,ori)),exsz));
    [tmpval rowval colval zval]=findmax3d(double(abs(rscorrim)));
    s2row=(rowval-(exsz/2+1))*cropsz/exsz;
    s2col=(colval-(exsz/2+1))*cropsz/exsz;
    s2z=(zval-(exsz/2+1))*cropsz/exsz;
    shift1=shift_row+s2row;
    shift2=shift_col+s2col;
    shift3=shift_z+s2z;
end