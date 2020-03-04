% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, August 2019
%
function genPupilfigs(obj,ImgType,datapath)
% genPuilfigs - generate figures of phase retrieval result, including various PSFs, pupil 
% functions, and plots of zernike coefficients. 
%
%   Input parameter: ImgType - type of image to be generated, select from
%   'PSF', 'pupil' and 'zernike'

switch ImgType
    case 'PSF'
        z = obj.Zpos;
        zind=[obj.Zindstart:obj.Zindstep:obj.Zindend];
        RC=64;
        L=length(zind);
        Modpsf=obj.PSFstruct.Modpsf(RC+1-RC/4:RC+1+RC/4,RC+1-RC/4:RC+1+RC/4,:);
        h1=[];
      
        figure('Color',[1,1,1],'Name','phase retrieved PSFs at sampled z positions','Resize','on','Units','normalized','Position',[0.3,0.3,0.43,0.1])

        for ii=1:L
            h1(ii)=subplot('position',[(ii-1)/(L+1),0.1,1/(L+1),0.8]);      %0.75
            image(double(squeeze(Modpsf(:,:,zind(ii)))),'CDataMapping','scaled','Parent',h1(ii))
            text(3,3,[num2str(z(zind(ii)),3),'\mum'],'color',[1,1,1]);

        end
        h1(ii+1)=subplot('position',[L/(L+1),0.1,1/(L+1),0.8]); %0.75
        image(double(permute(squeeze(Modpsf(17-10:17+10,17,:)),[2,1])),'CDataMapping','scaled','Parent',h1(ii+1))
        text(3,3,['x-z'],'color',[1,1,1]);
        colormap(jet)
        axis(h1,'equal')
        axis(h1,'off')

  
        saveas(gca,fullfile(datapath,'phase retrieved PSFs.tif'));   %edited by Fan Xu
        
    case 'pupil'
        figure('Color',[1,1,1],'Name',' phase retrieved and Zernike fitted pupil function','Resize','on','Units','normalized','Position',[0.3,0.3,0.22,0.42])
        h1=[];
        RC=64;
        Rsub=63;
        h1(1)=subplot('Position',[0,0.5,1/2,1/2]);
        image(double(obj.PRstruct.Pupil.mag(RC-Rsub:RC+Rsub,RC-Rsub:RC+Rsub)),'CDataMapping','scaled','Parent',h1(1))
        text(3,8,['PR pupil mag'],'color',[1,1,1],'FontSize',18);
        h1(2)=subplot('Position',[0,0,1/2,1/2]);
        image(double(obj.PRstruct.Fittedpupil.mag(RC-Rsub:RC+Rsub,RC-Rsub:RC+Rsub)),'CDataMapping','scaled','Parent',h1(2))
        text(3,8,['Zernike pupil mag'],'color',[1,1,1],'FontSize',18);
        h1(3)=subplot('Position',[0.5,0.5,1/2,1/2]);
        tmp=angle(obj.PRstruct.Pupil.phase);
        mag=obj.PRstruct.Pupil.mag;
        mag(mag>0)=1;
        PRphase=tmp.*mag;
        if obj.Enableunwrap==1
            PRphase=obj.PRstruct.Pupil.uwphase;
        end
        image(double(PRphase(RC-Rsub:RC+Rsub,RC-Rsub:RC+Rsub)),'CDataMapping','scaled','Parent',h1(3))
        text(3,8,['PR pupil phase'],'color',[0,0,0],'FontSize',18);
        h1(4)=subplot('Position',[0.5,0,1/2,1/2]);
        image(double(obj.PRstruct.Fittedpupil.phase(RC-Rsub:RC+Rsub,RC-Rsub:RC+Rsub)),'CDataMapping','scaled','Parent',h1(4))
        text(3,8,['Zernike pupil phase'],'color',[0,0,0],'FontSize',18);
        colormap(gray)
        axis(h1,'equal')
        axis(h1,'off')
     
        saveas(gca,fullfile(datapath,'phase retrieved and Zernike fitted pupil function.tif'));   %edited by Fan Xu
    case 'zernike'
        
        PlotZernikeC(obj.PRstruct.Zernike_phase,'phase');
        saveas(gca,fullfile(datapath,'phase_Zernike coefficinet.tif'));   %edited by Fan Xu
        
end
end

function PlotZernikeC(CN_phase,type)
nZ=length(CN_phase);
vec=linspace(max(CN_phase)-0.1,min(CN_phase)+0.1,8);
dinv=vec(1)-vec(2);
ftsz=12;
figure('position',[200,200,760,300])
plot(CN_phase,'o-')
text(nZ+3,vec(1),['x shift: ', num2str(CN_phase(2),'%.2f')],'fontsize',ftsz);
text(nZ+3,vec(2),['y shift: ', num2str(CN_phase(3),'%.2f')],'fontsize',ftsz);
text(nZ+3,vec(3),['z shift: ', num2str(CN_phase(4),'%.2f')],'fontsize',ftsz);
text(nZ+3,vec(4),['Astigmatism: ', num2str(CN_phase(5),'%.2f')],'fontsize',ftsz);
text(nZ+3,vec(5),['Astigmatism(45^o): ', num2str(CN_phase(6),'%.2f')],'fontsize',ftsz);
text(nZ+3,vec(6),['Coma(x): ', num2str(CN_phase(7),'%.2f')],'fontsize',ftsz);
text(nZ+3,vec(7),['Coma(y): ', num2str(CN_phase(8),'%.2f')],'fontsize',ftsz);
text(nZ+3,vec(8),['Spherical: ', num2str(CN_phase(9),'%.2f')],'fontsize',ftsz);
xlim([0,nZ+30])
ylim([min(CN_phase)-dinv/2,max(CN_phase)+dinv/2])
set(gca,'fontsize',ftsz)
xlabel('Zernike coefficient number','fontsize',ftsz)
ylabel('Value','fontsize',ftsz)
title(type)

end
