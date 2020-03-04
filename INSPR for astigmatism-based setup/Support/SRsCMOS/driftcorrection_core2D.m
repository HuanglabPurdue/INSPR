
function [shift1 shift2]=driftcorrection_core2D(im1,im2,iniguess)
corrim=imgaussfilt(single(normxcorr2(im1,im2)),1);  %old 1


if nargin>=3
    shift_col=iniguess(2);
    shift_row=iniguess(1);
    rowval=shift_row+(size(corrim,1)+1)./2;
    colval=shift_col+(size(corrim,2)+1)./2;
    cropsz2=24;  
    ori=[colval-cropsz2/2-1 rowval-cropsz2/2-1];
    smallim = imcrop(corrim, [ori(1)+1 ori(2)+1 cropsz2-1 cropsz2-1]);
    
    [tmpval rowval2 colval2]=findmax(double(smallim));
    colval=ori(1)+colval2;
    rowval=ori(2)+rowval2;
else
    [tmpval rowval colval]=findmax(double(corrim));
    
end
shift_row=(rowval-(size(corrim,1)+1)./2);
shift_col=(colval-(size(corrim,2)+1)./2);
cropsz=24;   
exsz2=512;
ori=[colval-cropsz/2-1 rowval-cropsz/2-1];

rstmp_fft = fftshift(fft2(imcrop(corrim, [ori(1)+1 ori(2)+1 cropsz-1 cropsz-1])));
rstmp_pad = padarray(rstmp_fft, [(exsz2-cropsz)/2 (exsz2-cropsz)/2]);
rscorrim = ifft2(rstmp_pad);

xLeft = (exsz2-128)/2+1;
xRight = (exsz2+128)/2;
yLeft = (exsz2-128)/2+1;
yRight = (exsz2+128)/2;
cropex = rscorrim(xLeft:xRight,yLeft:yRight);

exsz=128;
[tmpval rowval colval]=findmax(double(abs(cropex)));
s2row=(rowval-(exsz/2+1))*cropsz/exsz2;
s2col=(colval-(exsz/2+1))*cropsz/exsz2;
shift1=shift_row+s2row;
shift2=shift_col+s2col;