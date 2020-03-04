function OTFparam=findOTFparam(obj)
% findOTFparam - find SigmaX and SigmaY of a Gaussian filter for OTF rescale. 
%   They are found by fitting the ratio OTF with a 2D gaussian, the ratio OTF
%   is the ratio between the measured OTF and the phase retreived OTF
R=obj.PSFsize;
%z=[obj.Zstart:obj.Zstep:obj.Zend];
z = obj.Zpos;
N=obj.DatadimZ;
[val,ind]=min(abs(z));
realsize0=floor(obj.OTFratioSize/2);
realsize1=ceil(obj.OTFratioSize/2);
starty=-realsize0+R/2+1;endy=realsize1+R/2;
startx=-realsize0+R/2+1;endx=realsize1+R/2;

mOTF=fftshift(ifft2(obj.Mpsf_extend(:,:,ind)));% measured OTF
rOTF=fftshift(ifft2(obj.PSFstruct.PRpsf(:,:,ind))); % phase retrieved OTF, changed by FX, from ZKpsf to PRpsf
tmp=abs(mOTF)./abs(rOTF);
tmp1=tmp(startx:endx,starty:endy);
ratio=tmp1;
% fit_param=[1,2,2,0];
% [I,sigmax,sigmay,bg,fit_im]=GaussRfit(obj,fit_param,ratio);
% OTFparam=[I,sigmax,sigmay,bg];

fit_param=[1,2,0];
[I,sigma,bg,fit_im]=GaussRfit(obj,fit_param,ratio);
OTFparam=[I,sigma,bg];

% generate zernike fitted PSF modified by OTF rescale
Mod_psf=zeros(R,R,N);
for ii=1:N
    Fig4=obj.PSFstruct.PRpsf(:,:,ii);   %changed by FX, from ZKpsf to PRpsf
    Fig4=Fig4./sum(sum(Fig4));
    Mod_OTF=fftshift(ifft2(Fig4)).*fit_im;
    Fig5=abs(fft2(Mod_OTF));
    Mod_psf(:,:,ii)=Fig5;
end
% save SigmaX and SigmaY in PRstruct
% Fan test, SigmaX and SigmaY should be same?
obj.PRstruct.SigmaX=sigma;
obj.PRstruct.SigmaY=sigma;
% save modified PSF in PSFstruct
obj.PSFstruct.Modpsf=Mod_psf;
end

% function [I,sigmax,sigmay,bg,fit_im]=GaussRfit(obj,startpoint,input_im)
function [I,sigma,bg,fit_im]=GaussRfit(obj,startpoint,input_im)

R1=obj.OTFratioSize;
R=obj.PSFsize;
estimate=fminsearch(@(x) Gauss2(x,input_im,R1,R,obj.Pixelsize),startpoint,optimset('MaxIter',50,'Display','off'));

I=estimate(:,1);
% restrict sigmax and sigmay to be less than 5
% tmp=estimate(:,2);
% tmp(tmp>8)=5;
% sigmax=tmp;
% tmp=estimate(:,3);
% tmp(tmp>8)=5;
% sigmay=tmp;
tmp=estimate(:,2);
tmp(tmp>5)=5;
sigma = tmp;

bg=0;
scale=R*obj.Pixelsize;
[xx,yy]=meshgrid(-R/2:R/2-1,-R/2:R/2-1);

X=abs(xx)./scale;
Y=abs(yy)./scale;
% fit_im=I.*exp(-X.^2./2./sigmax^2).*exp(-Y.^2./2./sigmay^2)+bg;
fit_im=I.*exp(-X.^2./2./sigma^2).*exp(-Y.^2./2./sigma^2)+bg;
end


function [sse,Model]=Gauss2(x,input_im,R1,R,pixelsize)
I=x(1);
% sigmax=x(2);
% sigmay=x(3);
sigma=x(2);
bg=0;
[xx,yy]=meshgrid(-R1/2:R1/2-1,-R1/2:R1/2-1);

scale=R*pixelsize;
X=abs(xx)./scale;
Y=abs(yy)./scale;
% Model=I.*exp(-X.^2./2./sigmax^2).*exp(-Y.^2./2./sigmay^2)+bg;
Model=I.*exp(-X.^2./2./sigma^2).*exp(-Y.^2./2./sigma^2)+bg;

sse=sum(sum((Model-input_im).^2));
end