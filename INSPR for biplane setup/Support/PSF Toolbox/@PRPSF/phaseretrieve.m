function phaseretrieve(obj)
% phaseretrieve - generate pupil function based on a phase retrieval
% algorithm described in paper.
%z=[obj.Zstart:obj.Zstep:obj.Zend];
z = obj.Zpos;
zind=[obj.Zindstart:obj.Zindstep:obj.Zindend];
N=length(zind);
n=obj.PRstruct.RefractiveIndex;
Freq_max=obj.PRstruct.NA/obj.PRstruct.Lambda;
NA_constrain=obj.k_r<Freq_max;
k_z=sqrt((n/obj.PRstruct.Lambda)^2-obj.k_r.^2).*NA_constrain;
Fig=NA_constrain;
pupil_mag=Fig/sum(sum(Fig)); % initial pupil function, normalization

R=obj.PSFsize;
MpsfA=zeros(R,R,N);
RpsfA_phase=zeros(R,R,N);
Rpupil_mag=zeros(R,R,N);
Rpupil_phase=zeros(R,R,N);
pupil_phase=ones(R,R);
% weight = zeros(1,N); %Add by FX

for k=1:obj.IterationNum
    for j=1:N
        defocus_phase=2*pi*z(zind(j)).*k_z;
        pupil_complex=pupil_mag.*exp(defocus_phase.*1i).*pupil_phase;   %Apply defocus for each PSF section
        Fig1=abs(fftshift(fft2(pupil_complex))).^2;
        PSF0=Fig1./sum(sum(Fig1));
        Mpsfo=squeeze(obj.Mpsf_extend(:,:,zind(j)));
        
        % at iteration number greater than IterationNumK, add previous retrieved PSF information in measured PSF
        if k>obj.IterationNumK
            Mask=(Mpsfo==0);
            Mpsfo(Mask)=PSF0(Mask);
%             Mpsfo(Mask)=min(min(Mpsfo));

        end
        
        RpsfA=fft2(pupil_complex);  %FFT each complex section
        RpsfA_phase(:,:,j)=RpsfA./abs(RpsfA);
        Fig2=fftshift(sqrt(abs(Mpsfo)));
        MpsfA(:,:,j)=Fig2./sum(sum(Fig2));
        Rpupil=ifft2((MpsfA(:,:,j)).*RpsfA_phase(:,:,j));   %Replace calculated magnitude with measured PSF
%         Rpupil = wiener2(Rpupil,[5 5]); %FX
        
        Rpupil_mag(:,:,j)=abs(Rpupil);
        Rpupil_phase(:,:,j)=Rpupil./Rpupil_mag(:,:,j).*exp(-defocus_phase.*1i);
%         weight(1,j) =  exp(0*abs(z(j)));
%         Rpupil_phase(:,:,j)=Rpupil./Rpupil_mag(:,:,j).*exp(-defocus_phase.*1i) .* exp(0*abs(z(j))); %refocus for zero defocus

    end
%     A = fspecial('gaussian',5,1);
%     Rpupil_phase = imfilter(Rpupil_phase,A);
%     Rpupil_phase = wiener2(Rpupil_phase,[5 5]);
    % generate pupil phase
%     Fig5 = sum(Rpupil_phase,3) / sum(weight);
%     Fig5 = wiener2(Fig5,[5 5]);
%     A = fspecial('gaussian',5,1);
%     Fig5 = imfilter(Fig5,A);
    Fig5=mean(Rpupil_phase,3);
    pupil_phase=Fig5./abs(Fig5);
    
    %test filter
%     A = fspecial('gaussian',5,1);
%     tmp = imfilter(pupil_phase,A);
%     pupil_phase = tmp./abs(tmp);
%     if k > 5 && k < 15
%     if k <= 15
%         G = fspecial('Gaussian',9,2);
%         %meanPupil = conv2(abs(meanPupil),G,'same').*exp(1i.*angle(meanPupil));
% %         pupil_phase = imfilter(pupil_phase.*NA_constrain,G); %.*exp(1i.*angle(meanPupil));
%         pupil_phase = imfilter(pupil_phase,G); %.*exp(1i.*angle(meanPupil));
%     end

%     A = fspecial('gaussian',5,1);
%     Rpupil_mag = imfilter(Rpupil_mag,A);

    % generate pupil magnitude
    Fig3=mean(Rpupil_mag,3).*NA_constrain;  
%    Fig3 = wiener2(Fig3,[5 5]);
%     Fig3=NA_constrain;
    Fig4=abs(Fig5).*Fig3; % pupil magnitude before normalization
    Fig4=Fig4.^2;
    Fig4=Fig4./sum(sum(Fig4));
    pupil_mag=sqrt(Fig4); % pupil magnitude after normalization
    

%     if k<=15  
%        pupil_mag = Fig/sum(sum(Fig)); % added by FX, make magnitude flat
%     end
    
end
% generate phase retrieved PSF
psf=zeros(R,R,numel(z));

for j=1:numel(z)
    defocus_phase=2*pi*z(j).*k_z.*1i;
    pupil_complex=pupil_mag.*pupil_phase.*exp(defocus_phase);
    Fig2=abs(fftshift(fft2(pupil_complex))).^2;
    psf(:,:,j)=Fig2./R^2; % normalized PSF
end

% save pupil function and PSF in PRstruct
obj.PRstruct.Pupil.phase=pupil_phase;
obj.PRstruct.Pupil.mag=pupil_mag;
obj.PSFstruct.PRpsf=psf;
end