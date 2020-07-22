%  Script for 2D reconstruction
% (C) Copyright 2020                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, July 2020
%       
%%
function srobj = analysis2D_fromGaussian(recon, setup_para)


global recon_stop  %control program stop
% input:
%   recon: reconstruction parameters from GUI
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
    boxsz = 7;  % small subregion for 2D case
%     thresh = [recon.seg_thresh_low,recon.seg_thresh_high];  
    thresh = [recon.seg_thresh_low,recon.seg_thresh_low];      
    [subregion_ch1,subvar_ch1,frame_num,l,t] = crop_subregion_var_ast(ims,boxsz,thresh,setup_para);
        
    %localzation
    drawnow
    if recon_stop == 1
        return;
    end
    
    disp('2D localization');
        
    if recon.isGPU
        [sub_xM,sub_yM,sub_photonM,sub_bgM,crlbM,llM] = loc_gaussian_2D(subregion_ch1,subvar_ch1);
    else
        % message to GUI
        msgbox('Current 2D localization only supports GPU version, please select <Run GPU> checkbox!');
    end
    
    loc_x = [loc_x; sub_xM+l];
    loc_y = [loc_y; sub_yM+t];
    
    loc_t = [loc_t; frame_num + (nn-1)*size(ims,3)];
    loc_photons = [loc_photons; sub_photonM];
    loc_bg = [loc_bg; sub_bgM];
    loc_ll = [loc_ll; llM];
    loc_crlb = [loc_crlb; crlbM];
    
end

sobj = struct('loc_x',loc_x,'loc_y',loc_y,'loc_t',loc_t,'loc_photons',loc_photons,'loc_bg',loc_bg,'loc_ll',loc_ll,'loc_crlb',loc_crlb);
save(fullfile(setup.workspace,'loc_backup'), 'sobj');

%% remove low confidence spots

if recon.isRej == 1 && recon.is_bg == 0
    disp('Localization rejection'); 
    
    llmask = sobj.loc_ll > recon.rej.llthreshold;
    intmask = sobj.loc_photons < recon.rej.min_photon;
    
    totmask=llmask | intmask;
    
    
    loc_x_keep = sobj.loc_x(~totmask);
    loc_y_keep = sobj.loc_y(~totmask);
    loc_t_keep = sobj.loc_t(~totmask);
    
    loc_ll_keep = sobj.loc_ll(~totmask);
    loc_bg_keep = sobj.loc_bg(~totmask);
    loc_photons_keep = sobj.loc_photons(~totmask);
    loc_crlb_keep = sobj.loc_crlb(~totmask,:);
else 
    loc_x_keep = sobj.loc_x;
    loc_y_keep = sobj.loc_y;
    loc_t_keep = sobj.loc_t;
    
    loc_ll_keep = sobj.loc_ll;
    loc_bg_keep = sobj.loc_bg;
    loc_photons_keep = sobj.loc_photons;
    loc_crlb_keep = sobj.loc_crlb;
end


%% 2D alignment
pixelsize = setup.pixelsize;    %nm

srobj.loc_x = loc_x_keep;% pixel
srobj.loc_y = loc_y_keep;% pixel
srobj.loc_t = loc_t_keep;

srobj.loc_ll = loc_ll_keep;
srobj.loc_bg = loc_bg_keep;
srobj.loc_photons = loc_photons_keep;
srobj.loc_crlb = loc_crlb_keep;

if recon.isDC == 1
    disp('Drift correction');

    srobj.DC2D_fnum = 1000;
    srobj.DC2D_errthresh = 0.06;   % micron
    srobj.DC2D_interval = 10; % maximum segmented interval in analyzing drift     20
    srobj.Cam_pixelsz = pixelsize;
    srobj.loc_cycle = floor(loc_t_keep/srobj.DC2D_fnum);
    
    srobj.Perform_DriftCorrection2D()
    
    save(fullfile(setup.workspace,'recon_2DDC_backup'),'srobj');

    drawnow
    if recon_stop == 1
        return;
    end

end

save(fullfile(setup.workspace,'recon_final_backup'),'srobj');



