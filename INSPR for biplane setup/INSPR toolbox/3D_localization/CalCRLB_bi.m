% Script for calculating CRLB in biplane setup
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, December 2019

classdef CalCRLB_bi < handle
    % CalCRLB class for calculating CRLB of simulated emitters, given a PSF model defined by PRstruct
    %   create object: obj = CalCRLB(PRstruct,PSFtype)
    %   inputs required: PRstruct - define necessary parameters for a PSF model
    %                   PSFtype - type of method to generate PSFs for CRLB
    %                   calculation, options are 'pupil' and 'zernike'
    %
    % CalCRLB Properties (Implemented object):
    %   PSFobj - 
    %
    % CalCRLB Properties (Input):
    %   Bg - 
    %   Boxsize - 
    %   Deltax - 
    %   Deltaz - 
    %   PRstruct - 
    %   Photon - 
    %   Pixelsize - 
    %   Xpos - 
    %   Ypos - 
    %   Zpos - 
    %
    % CalCRLB Properties (Output):
    %   CRLB - 
    %   X_STD - 
    %   Y_STD - 
    %   Z_STD - 
    %   Photon_STD - 
    %   Bg_STD - 
    %
    % CalCRLB Methods:
    %   prepInputparam - generate parameters of PSFs used for CRLB calculation
    %   calcrlb - calculate CRLB of simulated emitters, given a PSF model
    %   genfigs - generate plots of theoretical localization precision in x, y and z at z positions defined by 'Zpos'
    %
    %   see also PSF_pupil
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
        PSFobj1;% object of PSF_pupil or PSF_zernike class, used for generating PSFs 
        PSFobj2;
        Xpos;% x positions of simulated emitters, a vector of N elements, unit is pixel
        Ypos;% y positions of simulated emitters, a vector of N elements, unit is pixel
        Zpos;% z positions of simulated emitters, a vector of N elements, unit is micron
        Photon;% photon counts of simulated emitters, a vector of N elements
        Bg;% background photon counts of simulated emitters, a vector of N elements
        Pixelsize;% pixel size at sample plane, unit is micron
        Boxsize;% image size of simulated emitter
        PSFtype;% type of method to generate PSFs for CRLB
        Deltax;% increment in x and y directions for calculating first and second derivative of the objective function, unit is pixel
        Deltaz;% increment in z directions forcalculating first and second derivative of the objective function, unit is micron 
        PN = 4; % number of fitting parameters in CRLB calculation except of background
        FisherM;% Fisher information matrix
        Planedis = 0.4;
    end
    properties (SetAccess = private, GetAccess = public)
        Xin; % parameters of PSFs, a N*PN x PN matrix, N is the number of elements in Xpos. PN is the number of fitting parameters, including x, y, z, photon, bg  
    end
    % output parameters
    properties (SetAccess = private, GetAccess = public)
        CRLB; % CRLB of simulated emmiters, a N x PN matrix, N is the number of elements in Xpos. PN is the number of fitting parameters, including x, y, z, photon, bg  
        X_STD; % theoretical localization precision in X dimension, a vector of N elements, unit is pixel
        Y_STD; % theoretical localization precision in Y dimension, a vector of N elements, unit is pixel
        Z_STD; % theoretical localization precision in Z dimension, a vector of N elements, unit is micron
        Photon_STD; % theoretical localization precision in photon count, a vector of N elements
        Bg_STD; % theoretical localization precision in background count, a vector of N elements
    end
    
    methods
        function obj = CalCRLB_bi(PRstruct1,PRstruct2)
                    obj.PSFobj1 = PSF_zernike(PRstruct1);
                    obj.PSFobj2 = PSF_zernike(PRstruct2);
        end
        
        function prepInputparam(obj)
            % prepInputparam - generate parameters of PSFs use for CRLB calculation.
            %   the out put is Xin. it is a N*PN x PN matrix, N is the
            %   number of elements in Xpos. PN is the number of fitting
            %   parameters, including x, y, z, photon, bg
            N = numel(obj.Xpos);
            pN = obj.PN;
            PIx0 = zeros(1,pN);
            PIy0 = zeros(1,pN);
            PIz0 = zeros(1,pN);
            
            PIx0(2) = obj.Deltax;
            PIy0(3) = obj.Deltax;
            PIz0(4) = obj.Deltaz;
            
            x = [];
            y = [];
            z = [];
            x0 = cat(2,obj.Xpos,obj.Ypos,obj.Zpos);
            for t = 1:pN
                x = cat(2,x,x0(:,1)+PIx0(t));
                y = cat(2,y,x0(:,2)+PIy0(t));
                z = cat(2,z,x0(:,3)+PIz0(t));
            end
            obj.Xin = cat(2,reshape(x',N*pN,1),reshape(y',N*pN,1),reshape(z',N*pN,1));
            
        end
        
        function calcrlb(obj,tform,offset_int) 
            % calcrlb - calculate CRLB of simulated emitters, given a PSF model.
            %   It uses PSF_pupil or PSF_zernike class to generate PSFs
            %   from given parameters in 'Xin'
            %
            %   see also PSF_pupil
            psfL = [];
            for nn = 1:2
                eval(['psfobj=obj.PSFobj',num2str(nn),';'])
                psfobj.Xpos = obj.Xin(:,1);
                psfobj.Ypos = obj.Xin(:,2);
                psfobj.Zpos = obj.Xin(:,3)+obj.Planedis*(nn-1);
                psfobj.Boxsize = obj.Boxsize;
                psfobj.Pixelsize = obj.Pixelsize;
                psfobj.precomputeParam();
%                 psfobj.genZernike();   
%                 psfobj.genPupil();  
                psfobj.setPupil();
                psfobj.genZernikeMag();
                psfobj.genPSF_2();
                psfobj.scalePSF('normal');
                psfI = psfobj.ScaledPSFs;
                tmp = psfobj.Pupil.mag;
                normf = sum(sum(tmp.^2,1),2);
                psf = psfI./normf;  
                                      
                psfL = cat(4,psfL,psf);
            end
            % calculate Fisher Information matrix
            N = numel(obj.Xpos);
            pN = obj.PN;
            pN0 = 7;
            funFi = zeros(obj.Boxsize,obj.Boxsize,pN0);
            FisherM = zeros(pN0,pN0,2);
            xVar = zeros(N,pN0);
            funFi2 = zeros(obj.Boxsize,obj.Boxsize,pN0,2);
            psfIni2 = zeros(obj.Boxsize,obj.Boxsize,2);
            for s = 0:N-1
                for nn = 1:2
                    t = s+1;
                    psf = psfL(:,:,:,nn);
                    if nn == 2
                        tform_tmp = tform;
                        tform_tmp.T(3,1) = tform.T(3,1) + offset_int(s+1,1);
                        tform_tmp.T(3,2) = tform.T(3,2) + offset_int(s+1,2);
                        
                        edge_1 = mean(squeeze(psf(1,:,s*pN+1:(s+1)*pN)),1)';
                        edge_2 = mean(squeeze(psf(:,1,s*pN+1:(s+1)*pN)),1)';
                        edge_3 = mean(squeeze(psf(end,:,s*pN+1:(s+1)*pN)),1)';
                        edge_4 = mean(squeeze(psf(:,end,s*pN+1:(s+1)*pN)),1)';
                        edge_fill = min([edge_1 edge_2 edge_3 edge_4],[],2);
                    
                        psf(:,:,s*pN+1:(s+1)*pN) = imwarp(psf(:,:,s*pN+1:(s+1)*pN),tform_tmp,...
                        'cubic','OutputView',imref2d(size(psf(:,:,s*pN+1:(s+1)*pN))),'FillValues',edge_fill);                          
                    end
                    
                    vec_bg = zeros(obj.Boxsize,obj.Boxsize,2);
                    vec_bg(:,:,nn) = 1;
                    vec_I = zeros(obj.Boxsize,obj.Boxsize,2);
                    vec_I(:,:,nn) = psf(:,:,s*pN+1);
                    %x
                    funFi(:,:,1) = obj.Photon(t,nn).*(psf(:,:,s*pN+2)-psf(:,:,s*pN+1))./obj.Deltax; %FX chabge
                    %y
                    funFi(:,:,2) = obj.Photon(t,nn).*(psf(:,:,s*pN+3)-psf(:,:,s*pN+1))./obj.Deltax;
                    %z
                    funFi(:,:,3) = obj.Photon(t,nn).*(psf(:,:,s*pN+4)-psf(:,:,s*pN+1))./obj.Deltaz;
                    %I
                    funFi(:,:,4:5) = vec_I;
                    %bg
                    funFi(:,:,6:7) = vec_bg;
                    for j = 1:pN0
                        for k = 1:pN0
%                             psfIni = psf(:,:,s*pN+1).*obj.Photon(t,nn)+obj.Bg(t,nn); 
                            psfIni = max(psf(:,:,s*pN+1).*obj.Photon(t,nn)+obj.Bg(t,nn), 1e-4); 
                            FisherM(j,k,nn) = sum(sum(funFi(:,:,j).*funFi(:,:,k)./psfIni));
                        end
                    end
                    funFi2(:,:,:,nn) = funFi;
                    psfIni2(:,:,nn) = psfIni;
                end
                FisherMsum = squeeze(sum(FisherM,3));
                %FisherMsum(abs(FisherMsum)<1e-12)=1e-12;
                LowerBi = pinv(FisherMsum);
                xVar(t,:) = diag(LowerBi)';
            end
            obj.CRLB = abs(xVar);
            obj.X_STD = sqrt(obj.CRLB(:,1));
            obj.Y_STD = sqrt(obj.CRLB(:,2));
            obj.Z_STD = sqrt(obj.CRLB(:,3));
            obj.Photon_STD = sqrt(obj.CRLB(:,4:5));
            obj.Bg_STD = sqrt(obj.CRLB(:,6:7));
            obj.FisherM.foo = funFi2;
            obj.FisherM.psf = psfIni2;
        end
        
        function genfigs(obj)
            % genfigs - generate plots of theoretical localization
            % precision in x, y and z at z positions defined by Zpos
            figure('position',[100,200,500,300])
            plot(obj.Zpos.*1e3,obj.X_STD.*obj.Pixelsize.*1e3,'r.-')
            hold on
            plot(obj.Zpos.*1e3,obj.Y_STD.*obj.Pixelsize.*1e3,'b.-')
            plot(obj.Zpos.*1e3,obj.Z_STD.*1e3,'g.-')
            axis tight;
            xlabel('z positions (nm)')
            ylabel('precison from CRLB (nm)')
            legend('\sigmax','\sigmay','\sigmaz')
        end
    end
    
end
