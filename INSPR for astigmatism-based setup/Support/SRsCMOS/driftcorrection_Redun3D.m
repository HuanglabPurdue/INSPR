
function [xout3,yout3,zout2,shifts]=driftcorrection_Redun3D(xout,yout,zout,tout,frmnum,reverseflag,currcy)
%% drift correction

if numel(unique(currcy))==1
    xout3=xout;
    yout3=yout;
    zout2=zout;
    shifts=[NaN NaN NaN];
    return
end
cutmeth='percentile';
pixelsz=50; % nm
thresh=25; % nm
[shiftx,shifty,R_shift3,flag,error1_final,error2_final,error3_final]=driftcorrection_Redun_core3D_parfor(...
    xout,yout,zout,pixelsz,thresh,cutmeth,frmnum,tout,currcy);

[xout2]=shiftcoords(xout,shiftx,tout,frmnum,reverseflag,currcy);
[yout2]=shiftcoords(yout,shifty,tout,frmnum,reverseflag,currcy);
[zout2]=shiftcoords(zout,R_shift3,tout,frmnum,reverseflag,currcy);



pixelsz=25; % nm
thresh=15; % nm

[shiftx_d,shifty_d]=driftcorrection_Redun_core2D_parfor(single(xout2),single(yout2),tout,frmnum,pixelsz,thresh,cutmeth,currcy);
[xout3]=shiftcoords(xout2,shiftx_d,tout,frmnum,reverseflag,currcy);
[yout3]=shiftcoords(yout2,shifty_d,tout,frmnum,reverseflag,currcy);

% save drftcorrresult xout2 yout2 zout2 R_shift1 R_shift2 R_shift3
shifts=[shiftx(:)+shiftx_d(:) shifty(:)+shifty_d(:) R_shift3(:)];