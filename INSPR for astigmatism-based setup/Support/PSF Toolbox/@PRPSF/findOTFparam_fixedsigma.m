function OTFparam=findOTFparam_fixedsigma(obj,sigma)
% findOTFparam - fixed SigmaX and SigmaY of a Gaussian filter for OTF rescale. 
R=obj.PSFsize;
%z=[obj.Zstart:obj.Zstep:obj.Zend];
N=obj.DatadimZ;
OTFparam=[1,sigma,0];

% FX test
scale=R*obj.Pixelsize;
[xx,yy]=meshgrid(-R/2:R/2-1,-R/2:R/2-1);
X=abs(xx)./scale;
Y=abs(yy)./scale;
fit_im=1.*exp(-X.^2./2./sigma^2).*exp(-Y.^2./2./sigma^2);

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
obj.PRstruct.SigmaX=sigma;
obj.PRstruct.SigmaY=sigma;
% save modified PSF in PSFstruct
obj.PSFstruct.Modpsf=Mod_psf;
end
