classdef PSF_pupil < handle
    % PSF_pupil class for generating PSFs from a pupil function
    %   create object: obj = PSF_pupil(PRstruct)
    %   input required: PRstruct - define necessary parameters for a PSF model
    %
    % PSF_pupil Properties (Input):
    %   Boxsize - 
    %   PRstruct - 
    %   Pixelsize - 
    %   Xpos - 
    %   Ypos - 
    %   Zpos - 
    %   ZposMed - 
    %   
    % PSF_pupil Properties (Output):
    %   PSFs - 
    %   ScaledPSFs - 
    %   IMMPSFs - 
    %
    % PSF_pupil Methods:
    %   precomputeParam - generate images for k space operation
    %   genPSF - generate PSFs from the given pupil function
    %   scalePSF - generate OTF rescaled PSFs
    %   genIMMPSF - generate PSFs considering refractive index mismatch aberration
    %
    %   see also OTFrescale
    properties
        % PRstruct - define necessary parameters for a PSF model
        %   NA
        %   Lambda
        %   RefractiveIndex
        %   Pupil: phase retrieved pupil function
        %           phase: phase image
        %           mag: magnitude image
        %   Zernike_phase: coefficient of zernike polynomials representing the pupil phase
        %   Zernike_mag: coefficient of zernike polynomials representing the pupil phase
        %   SigmaX: sigmax of Gaussian filter for OTF rescale, unit is 1/micron in k space, the conversion to real space is 1/(2*pi*SigmaX), unit is micron
        %   SigmaY: sigmay of Gaussian filter for OTF rescale, unit is 1/micron in k space, the conversion to real space is 1/(2*pi*SigmaY), unit is micron
        PRstruct;
        Xpos;% x positions of simulated emitters, a vector of N elements, unit is pixel
        Ypos;% y positions of simulated emitters, a vector of N elements, unit is pixel
        Zpos;% z positions of simulated emitters relative to the immersion medium, a vector of N elements, unit is micron
        ZposMed;% z positions of simulated emitters relative to the sample medium, a vector of N elements, unit is micron
        nMed;% refractive index of the sample medium
        PSFsize; % image size of the pupil function defined in PRstruct
        Boxsize; % image size of out put PSF
        Pixelsize;% pixel size at sample plane, unit is micron
    end
    
    properties (SetAccess = private, GetAccess = private)
        % precompute images for k space operation
        Zo;% r coordinates of out put PSF, it's a image of PSFsize x PSFsize
        k_r;% k_r coordinates of out put OTF, it's a image of PSFsize x PSFsize
        k_z;% k_z coordinates of out put OTF, it's a image of PSFsize x PSFsize
        Phi;% phi coordinates out put PSF, it's a image of PSFsize x PSFsize
        NA_constrain;% a circle function defining the limit of k_r, it's a image of PSFsize x PSFsize
        Cos1;% cos(theta1), theta1 is the angle of between the k vector and the optical axis in the sample medium
        Cos3;% cos(theta3), theta3 is the angle of between the k vector and the optical axis in the immersion medium
    end
    
    properties (SetAccess = private, GetAccess = public)
        % PSFs - out put PSFs from Fourier transform of the pupil function,
        % it's a 3D matrix of Boxsize x Boxsize x N, N is the number of
        % elements in Xpos.
        PSFs;
        % IMMPSFs - out put PSFs from Fourier transform of the pupil
        % function modified with the index mismatch aberration, it's a 3D
        % matrix of Boxsize x Boxsize x N, N is the number of elements in
        % Xpos.
        IMMPSFs;
        % ScaledPSFs - out put PSFs after OTF rescale, it's a 3D matrix of
        % Boxsize x Boxsize x N, N is the number of elements in Xpos.
        ScaledPSFs;
    end
    
    methods
        function obj=PSF_pupil(PRstruct)
            sz=size(PRstruct.Pupil.mag);
            obj.PSFsize=sz(1);
            obj.PRstruct=PRstruct;
        end
        
        function precomputeParam(obj)
            % precomputeParam - generate images for k space operation, and saved in
            % precomputed parameters.
            [X,Y]=meshgrid(-obj.PSFsize/2:obj.PSFsize/2-1,-obj.PSFsize/2:obj.PSFsize/2-1);
            obj.Zo=sqrt(X.^2+Y.^2);
            scale=obj.PSFsize*obj.Pixelsize;
            obj.k_r=obj.Zo./scale;
            obj.Phi=atan2(Y,X);
            n=obj.PRstruct.RefractiveIndex;
            Freq_max=obj.PRstruct.NA/obj.PRstruct.Lambda;
            obj.NA_constrain=obj.k_r<Freq_max;
            sin_theta3=obj.k_r.*obj.PRstruct.Lambda./n;
            if isempty(obj.nMed)
                obj.nMed=n;
            end
            sin_theta1=n./obj.nMed.*sin_theta3;
            obj.Cos1=sqrt(1-sin_theta1.^2);
            obj.Cos3=sqrt(1-sin_theta3.^2);

            obj.k_z=sqrt((n/obj.PRstruct.Lambda)^2-obj.k_r.^2).*obj.NA_constrain;
        end
        
        function genPSF(obj)
            % genPSF - generate PSFs from the given pupil function.
            %   The PSFs are directly calculated from the Fourier transform
            %   of pupil functions modified by shift phase in x, y and
            %   defocus phase in z. The out put is 'PSFs'
            N=numel(obj.Xpos);
            R=obj.PSFsize;
            Ri=obj.Boxsize;
            psfs=zeros(Ri,Ri,N);
            for ii=1:N
                shiftphase=-obj.k_r.*cos(obj.Phi).*obj.Xpos(ii).*obj.Pixelsize-obj.k_r.*sin(obj.Phi).*obj.Ypos(ii).*obj.Pixelsize;
                shiftphaseE=exp(-1i.*2.*pi.*shiftphase);
                defocusphaseE=exp(2.*pi.*1i.*obj.Zpos(ii).*obj.k_z);
                pupil_complex=obj.PRstruct.Pupil.mag.*obj.PRstruct.Pupil.phase.*shiftphaseE.*defocusphaseE;
                psfA=abs(fftshift(fft2(pupil_complex)));
                Fig2=psfA.^2;
                realsize0=floor(Ri/2);
                realsize1=ceil(Ri/2);
                startx=-realsize0+R/2+1;endx=realsize1+R/2;
                starty=-realsize0+R/2+1;endy=realsize1+R/2;
                psfs(:,:,ii)=Fig2(startx:endx,starty:endy)./R^2;
            end
            obj.PSFs=psfs;
        end
        
        function genIMMPSF(obj)
            % genIMMPSF - generate PSFs considering refractive index mismatch aberration.
            %   For oil immersion objective lens, when imaging a emitter at
            %   a distance from the cover glass, because of the index
            %   mismatch between the sample medium and the immersion oil,
            %   there will be an additional aberration phase added to the
            %   pupil function. Here for simplification, assume the camera
            %   is at the designed position and the immersion oil has the
            %   same refractive index with the coverglass. The PSFs are
            %   generated from emitters at z positions in the sample
            %   medium, defined in 'ZposMed'. The aberration phase is
            %   calculated at certain stage position, which is the first
            %   element in 'Zpos'.
            n=obj.PRstruct.RefractiveIndex;
            stagepos=obj.Zpos(1);
            % reference z position
            depth=stagepos*obj.nMed/n;
            N=numel(obj.Xpos);
            deltaH=depth*obj.nMed.*obj.Cos1-depth*n^2/obj.nMed.*obj.Cos3;
            % aberration phase from index mismatch
            IMMphase=exp(2*pi/obj.PRstruct.Lambda.*deltaH.*obj.NA_constrain.*1i);
            % z position of the emitters inside the sample medium
            zMed=obj.ZposMed;
            zMed(zMed<-depth)=-depth;
            
            R=obj.PSFsize;
            Ri=obj.Boxsize;
            psfs=zeros(Ri,Ri,N);
            for ii=1:N
                defocusMed=exp(2*pi/obj.PRstruct.Lambda*obj.nMed*zMed(ii).*obj.Cos1.*1i);
                shiftphase=-obj.k_r.*cos(obj.Phi).*obj.Xpos(ii).*obj.Pixelsize-obj.k_r.*sin(obj.Phi).*obj.Ypos(ii).*obj.Pixelsize;
                shiftphaseE=exp(-1i.*2.*pi.*shiftphase);
                pupil_complex=obj.PRstruct.Pupil.mag.*obj.PRstruct.Pupil.phase.*shiftphaseE.*defocusMed.*IMMphase;
                psfA=abs(fftshift(fft2(pupil_complex)));
                Fig2=psfA.^2;
                realsize0=floor(Ri/2);
                realsize1=ceil(Ri/2);
                startx=-realsize0+R/2+1;endx=realsize1+R/2;
                starty=-realsize0+R/2+1;endy=realsize1+R/2;
                psfs(:,:,ii)=Fig2(startx:endx,starty:endy)./R^2;
            end
            obj.IMMPSFs=psfs;
        end
        
        function scalePSF(obj)
            % scalePSF - generate OTF rescaled PSFs
            %   It operates 'PSFs' using the OTFrescale class. The OTF
            %   rescale acts as a 2D Gaussian filter, the resulting PSFs
            %   are smoother than the orignal PSFs.
            %
            %   see also OTFrescale
            otfobj=OTFrescale;
            otfobj.SigmaX=obj.PRstruct.SigmaX;
            otfobj.SigmaY=obj.PRstruct.SigmaY;
            otfobj.Pixelsize=obj.Pixelsize;
            otfobj.PSFs=obj.PSFs;
            otfobj.scaleRspace();
            obj.ScaledPSFs=otfobj.Modpsfs;
        end
    end
    
end

