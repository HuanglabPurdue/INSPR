%  Script for 3D reconstruction
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, July 2020
%       
%%
function srobj = analysis3D_fromPupil_ast(recon, setup_para)

global recon_stop  %control program stop
% input:
%   recon: reconstruction parameters from GUI,including data and pupil
%   model
%   setup: setup parameters

% output:
%   srobj: reconstruction results


%% setup parameters

setup.pixelsize = setup_para.Pixelsize * 1000;  %nm
setup.n_imm = setup_para.RefractiveIndex;
setup.n_sample = setup_para.nMed;
setup.workspace = setup_para.workspace;

srobj = SRscmos(setup.workspace); 

%% import data, biplane registration and localization

loc_x = [];
loc_y = [];
loc_z = [];
loc_t = [];
loc_photons = [];
loc_bg = [];
loc_ll = [];
loc_crlb = [];
loc_step = [];
for nn = 1:recon.dirN
    
    drawnow
    if recon_stop == 1
        return;
    end
    
    disp(['...Processing Data   ' int2str(nn)]);    
    if recon.isNewdata == 1
        load([recon.datapath, recon.datafile_name{nn}]);
    else
        ims = recon.ims;
    end
    

    %find subregion
    drawnow
    if recon_stop == 1
        return;
    end
    
    disp('Image segmentation');
    boxsz = 16;
    thresh = [recon.seg_thresh_low,recon.seg_thresh_high];      
    [subregion_ch1,subvar_ch1,frame_num,l,t] = crop_subregion_var_ast(ims,boxsz,thresh,setup_para,recon);
        
    %localzation
    drawnow
    if recon_stop == 1
        return;
    end
    
    disp('3D localization');
    
    step_num = mod(nn-1,recon.loopn) + 1;  
    
    if recon.isGPU
        [PM,crlbM,errM] = loc_ast_model(subregion_ch1,recon.probj_all{step_num},subvar_ch1,setup_para);
    else
        % message to GUI
        numcores = feature('numcores');
        msgbox(['Running in CPU version: the estimated time is ' num2str(round(length(t)/1000*numcores*0.3)) ' min in this cycle']);
        [PM,crlbM,errM] = loc_ast_model_CPU(subregion_ch1,recon.probj_all{step_num},setup_para);
    end
    loc_x = [loc_x; PM(:,2)+l];
    loc_y = [loc_y; PM(:,1)+t];
    
    offset = -(recon.probj_all{step_num}.Zpos(1)+ recon.probj_all{step_num}.Zpos(end))/2;
    loc_z = [loc_z; PM(:,3)+offset];
    loc_t = [loc_t; frame_num + (nn-1)*size(ims,3)];
    loc_photons = [loc_photons; PM(:,4)];
    loc_bg = [loc_bg; PM(:,5)];
    loc_ll = [loc_ll; errM(:,2)];
    loc_crlb = [loc_crlb; crlbM];
    
    loc_step = [loc_step; step_num * ones(size(frame_num,1),1)];

end

sobj = struct('loc_x',loc_x,'loc_y',loc_y,'loc_z',loc_z,'loc_t',loc_t,'loc_photons',loc_photons,'loc_bg',loc_bg,'loc_ll',loc_ll,'loc_crlb',loc_crlb, 'loc_step', loc_step);
save(fullfile(setup.workspace,'loc_backup'), 'sobj');

%% remove low confidence spots

if recon.isRej == 1
    disp('Localization rejection'); 
    
    llmask = sobj.loc_ll > recon.rej.llthreshold;
    intmask = sobj.loc_photons < recon.rej.min_photon;
    uncermask = sqrt(sobj.loc_crlb(:,3)) > recon.rej.loc_uncer_max;
    zmask= sobj.loc_z > recon.rej.zmask_high | sobj.loc_z < recon.rej.zmask_low;
    
    totmask=llmask | intmask | uncermask  | zmask;
    
    if recon.is_bg == 1
        totmask= uncermask  | zmask;
    end
    
    loc_x_keep = sobj.loc_x(~totmask);
    loc_y_keep = sobj.loc_y(~totmask);
    loc_z_keep = sobj.loc_z(~totmask);
    loc_t_keep = sobj.loc_t(~totmask);
    loc_step_keep = sobj.loc_step(~totmask);
    
    loc_ll_keep = sobj.loc_ll(~totmask);
    loc_bg_keep = sobj.loc_bg(~totmask);
    loc_photons_keep = sobj.loc_photons(~totmask);
    loc_crlb_keep = sobj.loc_crlb(~totmask,:);
else 
    loc_x_keep = sobj.loc_x;
    loc_y_keep = sobj.loc_y;
    loc_z_keep = sobj.loc_z;
    loc_t_keep = sobj.loc_t;
    loc_step_keep = sobj.loc_step;
    
    loc_ll_keep = sobj.loc_ll;
    loc_bg_keep = sobj.loc_bg;
    loc_photons_keep = sobj.loc_photons;
    loc_crlb_keep = sobj.loc_crlb;
end


%% 2D or 3D alignment
pixelsize = setup.pixelsize;    %nm

srobj.loc_x = loc_x_keep;% pixel
srobj.loc_y = loc_y_keep;% pixel
srobj.loc_z = loc_z_keep.*1e3 + recon.dc.z_offset;% nm, must be positive value
srobj.loc_t = loc_t_keep;
srobj.loc_step = loc_step_keep;

srobj.loc_ll = loc_ll_keep;
srobj.loc_bg = loc_bg_keep;
srobj.loc_photons = loc_photons_keep;
srobj.loc_crlb = loc_crlb_keep;

if recon.isDC == 1
    disp('Drift correction');

    srobj.frmpfile = recon.dc.frmpfile;
    srobj.loc_cycle = floor(loc_t_keep/srobj.frmpfile);
    srobj.Cam_pixelsz = pixelsize;
    
    srobj.Perform_DriftCorrection3D()
    
    save(fullfile(setup.workspace,'recon_3DDC_backup'),'srobj');

    drawnow
    if recon_stop == 1
        return;
    end

    n_imm = setup.n_imm;
    n_sample = setup.n_sample;
    srobj.step_ini = recon.dc.step_ini.*n_sample/n_imm;
    
    srobj.Perform_stackalignment()
end

save(fullfile(setup.workspace,'recon_final_backup'),'srobj');



