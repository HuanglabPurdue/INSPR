% Script for estimating pupil
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, August 2019
% 
function [probj,empupil] = PRPSF_aber_fromAveZ(ims_Zcal_ave_plane1,ims_Zcal_ave_plane2,index_record_Zplanes,empupil)


%% create object and set input properties of PRPSF class
disp('Generate pupil function and estimate aberration based on Phase retrieved method');
probj = PRPSF();
probj.CCDoffset = 0;
probj.Gain = 1;
probj.PRstruct.NA = empupil.NA; 
probj.PRstruct.Lambda = empupil.Lambda;
probj.PRstruct.RefractiveIndex = empupil.nMed;
probj.Pixelsize = empupil.Pixelsize;  
probj.PSFsize = 128;
probj.SubroiSize = empupil.imsz;    
probj.OTFratioSize = 60;
probj.ZernikeorderN = empupil.ZernikeorderN;  %Zernike coefficient
probj.Enableunwrap = 0; 
probj.IterationNum = 25;  
probj.IterationNumK = 5;  

%% prepare the PSFs and its positions
display('Prepare PR data...');
Zpos_plane1 = empupil.Z_pos - empupil.biplane_dist/2 + empupil.zshift;  %the Z position in two channels
Zpos_plane2 = empupil.Z_pos + empupil.biplane_dist/2 + empupil.zshift;

display(['Z-position index: ' num2str(index_record_Zplanes')]);

Zpos_plane1_sel = Zpos_plane1(index_record_Zplanes==1);    
Zpos_plane2_sel = Zpos_plane2(index_record_Zplanes==1);
ims_Zcal_ave_plane1_sel = ims_Zcal_ave_plane1(:,:,index_record_Zplanes==1);
ims_Zcal_ave_plane2_sel = ims_Zcal_ave_plane2(:,:,index_record_Zplanes==1);

probj.Zindstart = 1;
probj.Zindend = 2 * size(Zpos_plane1_sel,2);
probj.Zindstep = 1;

[Zpos_sort, index_sort] =  sort( cat(2, Zpos_plane1_sel, Zpos_plane2_sel) ); 
probj.Zpos = Zpos_sort;

%% generate PR result
temp_data = cat(3, ims_Zcal_ave_plane1_sel, ims_Zcal_ave_plane2_sel);
input_ZData = temp_data(:,:,index_sort);
probj.BeadData =  input_ZData;    
probj.prepdata();
probj.precomputeParam();


%% First PR
display('First time PR...');

probj.genMpsf();
probj.phaseretrieve();
probj.genZKresult();
probj.findOTFparam();

if empupil.XYshift_mode == 1
    
    [input_ZData, Zpos] = realign_Z_ave_img(probj.BeadData,probj.Zpos,probj);
    probj.Zpos = Zpos;
    probj.Zindend = size(input_ZData,3);
elseif empupil.XYshift_mode == 0
    [input_ZData, Zpos] = merge_Z_ave_img(probj.BeadData,probj.Zpos);

    probj.Zpos = Zpos;
    probj.Zindend = size(input_ZData,3);
end


display(['Shift mode in EMpupil is: ' num2str(empupil.Zshift_mode)]);

%% shift XYZ
if empupil.Zshift_mode

    for ii = 1 :3
        C4 = probj.PRstruct.Zernike_phase(4);
        est = fminsearch(@(x)probj.fitdefocus(x,C4),[0.2,1]);% calculate defocus
        zshift_tmp = est(1);

        CXY = probj.PRstruct.Zernike_phase([2,3]);
        xyshift = CXY*probj.PRstruct.Lambda/(2*pi*probj.Pixelsize*probj.PRstruct.NA);% calculate XY shift
        
        tmp_input = [];
        for i = 1 : size(input_ZData,3)
            tmp_input(:,:,i) = FourierShift2D(input_ZData(:,:,i), [-xyshift(2) -xyshift(1)]);
        end
        probj.BeadData = tmp_input;
        empupil.zshift = empupil.zshift + zshift_tmp;
        
        probj.Zpos = probj.Zpos + zshift_tmp;
        
        probj.prepdata();
        
        %% Second PR
        probj.genMpsf();
        probj.phaseretrieve();
        probj.genZKresult();
        probj.findOTFparam();
        
        if empupil.XYshift_mode == 1          
            [input_ZData, Zpos] = realign_Z_ave_img(probj.BeadData,probj.Zpos,probj);          
            probj.Zpos = Zpos;
            probj.Zindend = size(input_ZData,3);
        elseif empupil.XYshift_mode == 0
            [input_ZData, Zpos] = merge_Z_ave_img(probj.BeadData,probj.Zpos);
            probj.Zpos = Zpos;
            probj.Zindend = size(input_ZData,3);
        end
        
    end

end 

if empupil.Zshift_mode == 2
    empupil.zshift = 0;
end

display(['Z shift is: ' num2str(empupil.zshift)]);

display(['Updated Zernike(Z5-Z9): ' num2str(probj.PRstruct.Zernike_phase(5:9))]);

