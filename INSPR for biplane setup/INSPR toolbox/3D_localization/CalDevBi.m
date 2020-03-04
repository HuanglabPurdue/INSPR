% Script for calculating image derivatives in biplane setup
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, December 2019

classdef CalDevBi < handle
    
    properties
        PSFobj;% object of PSF_pupil or PSF_zernike class, used for generating PSFs
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
        PN = 7;
    end
    
    properties (SetAccess = private, GetAccess = public)
        Xin; % parameters of PSFs, a N*PN x PN matrix, N is the number of elements in Xpos. PN is the number of fitting parameters, including x, y, z, photon, bg  
    end
    % output parameters
    properties (SetAccess = private, GetAccess = public)
        Dev;
        Dev2;
        PSF0;
    end
    
    methods
        function obj = CalDevBi(PRstruct,PSFtype)
            switch PSFtype
                case 'pupil'
                    obj.PSFobj = PSF_pupil(PRstruct);
                case 'zernike'
                    obj.PSFobj = PSF_zernike(PRstruct);
                case 'IMM'
                    obj.PSFobj = PSF_zernike(PRstruct);
                case '4Pi'
                    obj.PSFobj = PSF_4pi(PRstruct);
                case 'interp'
                    obj.PSFobj = PSF_interp(PRstruct);
                case 'gauss'
                    obj.PSFobj = PSF_gauss(PRstruct);
                case 'DH'
                    obj.PSFobj = PSF_DH(PRstruct);
            end
            obj.PSFtype = PSFtype;
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
            
            PIx0([2,5]) = [obj.Deltax,2*obj.Deltax];
            PIy0([3,6]) = [obj.Deltax,2*obj.Deltax];
            PIz0([4,7]) = [obj.Deltaz,2*obj.Deltaz];
            
            x=[];
            y=[];
            z=[];

            x0=cat(2,obj.Xpos,obj.Ypos,obj.Zpos);
            for t=1:pN
                x = cat(2,x,x0(:,1)+PIx0(t));
                y = cat(2,y,x0(:,2)+PIy0(t));
                z = cat(2,z,x0(:,3)+PIz0(t));
            end
            obj.Xin = cat(2,reshape(x',N*pN,1),reshape(y',N*pN,1),reshape(z',N*pN,1));

        end
        
        function caldev(obj,tform,offset_int,flag)
            % caldev - calculate derivatives of simulated emitters, given a PSF model.
            %   It uses PSF_pupil or PSF_zernike class to generate PSFs
            %   from given parameters in 'Xin'
            %
            %   see also PSF_pupil
            obj.PSFobj.Xpos=obj.Xin(:,1);
            obj.PSFobj.Ypos=obj.Xin(:,2);
       
            if strcmp(obj.PSFtype,'IMM')
                obj.PSFobj.ZposMed = obj.Xin(:,3);
            else
                obj.PSFobj.Zpos=obj.Xin(:,3);
            end
            obj.PSFobj.Boxsize=obj.Boxsize;
            obj.PSFobj.Pixelsize=obj.Pixelsize;
            if ~sum(strcmp(obj.PSFtype,{'interp','gauss'}))
                obj.PSFobj.precomputeParam();
            end
            switch obj.PSFtype
                case 'zernike'
%                     obj.PSFobj.genPupil();
                    obj.PSFobj.genPSF_2();
                    obj.PSFobj.scalePSF('normal');
                case 'IMM'
                    obj.PSFobj.genPupil();
                    obj.PSFobj.genIMMPSF();
                    obj.PSFobj.scalePSF('IMM');
                case '4Pi'
                    obj.PSFobj.genPupil();
                    obj.PSFobj.genPupil_4pi('noIMMaber')
                    obj.PSFobj.genPSF(obj.PSFobj.Pupil4pi.p2);
                    obj.PSFobj.scalePSF();
                case 'pupil'
                    obj.PSFobj.genPSF();
                    obj.PSFobj.scalePSF();
                case 'interp'
                    obj.PSFobj.genSamplePSF();
                    obj.PSFobj.genPSF();
                case 'gauss'
                    obj.PSFobj.genPSF();
                case 'DH'
                    obj.PSFobj.genPupil();
                    obj.PSFobj.genPSF();
            end
            
            psfI = obj.PSFobj.ScaledPSFs;
            if ~sum(strcmp(obj.PSFtype,{'interp','gauss'}))
                tmp = obj.PSFobj.Pupil.mag;
                normf = sum(sum(tmp.^2,1),2);
            end
            switch obj.PSFtype
                case 'zernike'
                    psf = psfI./normf;
                case 'IMM'
                    psf = psfI./normf;
                case '4Pi'
                    psf = psfI./4;
                case 'pupil'
                    psf = psfI;
                case 'interp'
                    psf = psfI;
                case 'gauss'
                    psf = psfI;
                case 'DH'
                    psf = psfI;
            end
            
            if flag == 2
                for ii = 1 : size(offset_int,1)
                    tform_tmp = tform;
                    tform_tmp.T(3,1) = tform.T(3,1) + offset_int(ii,1);
                    tform_tmp.T(3,2) = tform.T(3,2) + offset_int(ii,2);
                    
                    
                    edge_1 = mean(squeeze(psf(1,:,(ii-1)*obj.PN+1:ii*obj.PN)),1)';
                    edge_2 = mean(squeeze(psf(:,1,(ii-1)*obj.PN+1:ii*obj.PN)),1)';
                    edge_3 = mean(squeeze(psf(end,:,(ii-1)*obj.PN+1:ii*obj.PN)),1)';
                    edge_4 = mean(squeeze(psf(:,end,(ii-1)*obj.PN+1:ii*obj.PN)),1)';
                    edge_fill = min([edge_1 edge_2 edge_3 edge_4],[],2);
                    
                    psf(:,:,(ii-1)*obj.PN+1:ii*obj.PN) = imwarp(psf(:,:,(ii-1)*obj.PN+1:ii*obj.PN),tform_tmp,...
                        'cubic','OutputView',imref2d(size(psf(:,:,(ii-1)*obj.PN+1:ii*obj.PN))),'FillValues',edge_fill);

                end
            end
%             % calculate numerical derivatives, 1st: dev, 2nd: dev2
            R = size(psf,1);
            N=numel(obj.Xpos);
            pN=obj.PN;
            dev = struct('x',zeros(R,R,N),'y',zeros(R,R,N),'z',zeros(R,R,N),'I',zeros(R,R,N),'bg',zeros(R,R,N));
            dev2 = struct('x',zeros(R,R,N),'y',zeros(R,R,N),'z',zeros(R,R,N),'I',zeros(R,R,N),'bg',zeros(R,R,N));
            for s=0:N-1
                t=s+1;
                %x
                dev.x(:,:,t) = obj.Photon(t).*(psf(:,:,s*pN+2)-psf(:,:,s*pN+1))./obj.Deltax;
                dev2.x(:,:,t) = obj.Photon(t).*(psf(:,:,s*pN+5)-2.*psf(:,:,s*pN+2)+psf(:,:,s*pN+1))./obj.Deltax./obj.Deltax;
                %y
                dev.y(:,:,t) = obj.Photon(t).*(psf(:,:,s*pN+3)-psf(:,:,s*pN+1))./obj.Deltax;
                dev2.y(:,:,t) = obj.Photon(t).*(psf(:,:,s*pN+6)-2.*psf(:,:,s*pN+3)+psf(:,:,s*pN+1))./obj.Deltax./obj.Deltax;
                %z
                dev.z(:,:,t) = obj.Photon(t).*(psf(:,:,s*pN+4)-psf(:,:,s*pN+1))./obj.Deltaz;
                dev2.z(:,:,t) = obj.Photon(t).*(psf(:,:,s*pN+7)-2.*psf(:,:,s*pN+4)+psf(:,:,s*pN+1))./obj.Deltaz./obj.Deltaz;
                %I
                dev.I(:,:,t) = psf(:,:,s*pN+1);
                dev2.I(:,:,t) = 0;
                %bg
                dev.bg(:,:,t) = 1;
                dev2.bg(:,:,t) = 0;
                
                
            end
            obj.Dev = dev;
            obj.Dev2 = dev2;
            obj.PSF0 = psf(:,:,[1:pN:N*pN]);
        end
        
        function [dL,dL2] = caldevLL(obj,data)
            namei = fields(obj.Dev);
            Np = numel(namei);
            N = size(data,3);
            dL = zeros(N,Np);
            dL2 = zeros(N,Np);
            parfor ii = 1:N
                datai = data(:,:,ii);
                psfI = obj.PSF0(:,:,ii).*obj.Photon(ii) + obj.Bg(ii);
                
                for nn = 1:Np
                    dev = obj.Dev.(namei{nn})(:,:,ii);
                    dev2 = obj.Dev2.(namei{nn})(:,:,ii);
                    dev(dev==0) = 1e-6;
                    dev2(dev2==0) = 1e-6;
                    tmp = (datai./psfI-1).*dev;
                    dL(ii,nn) = sum(tmp(:));
                    tmp = -1.*datai./psfI./psfI.*dev.*dev;% + (datai./psfI-1).*dev2; must change
                    dL2(ii,nn) = sum(tmp(:));
                end
                
            end
        end
        
        function [step] = getnextstep(obj,psfobj,data,I,bg,planedis,tform,offset_int)
            Nfit = size(I,1);
            Nplane = size(I,2);
            step = struct('x',zeros(Nfit,1),'y',zeros(Nfit,1),'z',zeros(Nfit,1),'I',zeros(Nfit,Nplane),'bg',zeros(Nfit,Nplane));
            namei = fields(step);
            dL_pos = zeros(Nfit,3);
            dL2_pos = zeros(Nfit,3);
            zpos = obj.Zpos;
            for ff = 1:Nplane
                obj.Photon = I(:,ff);
                obj.Bg = bg(:,ff);
                obj.Zpos = zpos+planedis(ff);
                obj.PSFobj = psfobj{ff};
                obj.prepInputparam();
                obj.caldev(tform,offset_int,ff); 
                [dL,dL2] = obj.caldevLL(data(:,:,:,ff));
                for nn = 1:3
                    dL_pos(:,nn) = dL_pos(:,nn)+dL(:,nn);
                    dL2_pos(:,nn) = dL2_pos(:,nn)+dL2(:,nn);
                end
                step.I(:,ff) = -1.*dL(:,4)./dL2(:,4);
                step.bg(:,ff) = -1.*dL(:,5)./dL2(:,5);
            end
            for nn = 1:3
                step.(namei{nn}) = -1.*dL_pos(:,nn)./dL2_pos(:,nn);
            end
        end
    end
    
end






