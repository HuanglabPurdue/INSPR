% Script for segmentation
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, July 2020
%%  
function [subregion_ch1,subregion_ch2,subvar_ch1,subvar_ch2,frame,l,t,offset_seg] = crop_subregion_without_transData(qd1,qd2,tform,boxsz,thresh,setup,recon)

if recon.isNewdata == 1
    if setup.is_sCMOS   %sCMOS case
        % sCMOS parameters
        offsetim_ch1 = repmat(setup.sCMOS_input.ccdoffset_ch1,[1 1 size(qd1,3)]);
        offsetim_ch2 = repmat(setup.sCMOS_input.ccdoffset_ch2,[1 1 size(qd2,3)]);
        
        varim_ch1 = repmat(setup.sCMOS_input.ccdvar_ch1,[1 1 size(qd1,3)]);
        varim_ch2 = repmat(setup.sCMOS_input.ccdvar_ch2,[1 1 size(qd2,3)]);
        
        gainim_ch1 = repmat(setup.sCMOS_input.gain_ch1,[1 1 size(qd1,3)]);
        gainim_ch2 = repmat(setup.sCMOS_input.gain_ch2,[1 1 size(qd2,3)]);
        
        qd1_in = (qd1 - offsetim_ch1) ./ gainim_ch1;
        qd2_in = (qd2 - offsetim_ch2) ./ gainim_ch2;
    else    %EMCCD case
        qd1_in = (qd1 - setup.offset) /setup.gain;
        qd2_in = (qd2 - setup.offset) /setup.gain;
    end
    
    qd1_in(qd1_in<=0) = 1e-6;
    qd2_in(qd2_in<=0) = 1e-6;
    
else
    qd1_in = qd1;
    qd2_in = qd2;
end


% background subtraction
if (recon.isNewdata == 1 && recon.is_bg == 1) || (recon.isNewdata == 0 && setup.is_bg == 0 && recon.is_bg == 1)
    filter_n = 101;
    bg_img_1 = medfilt1(qd1_in,filter_n,[],3);
    bg_img_2 = medfilt1(qd2_in,filter_n,[],3);
    
    subtract_img_1 = qd1_in - bg_img_1;
    subtract_img_2 = qd2_in - bg_img_2;
    subtract_img_1(subtract_img_1<=0) = 1e-6;
    subtract_img_2(subtract_img_2<=0) = 1e-6;
    
    qd1_in = subtract_img_1;
    qd2_in = subtract_img_2;
end

%transfer qd2 to find center peak
qd2_trans = imwarp(qd2_in,tform,'cubic','OutputView',imref2d(size(qd2_in)));

imsz_original = size(qd1_in);
rangemin=[5, 5];
rangemax=[imsz_original(1)-5, imsz_original(1)-5];

ims_ch1 = single(qd1_in(rangemin(1):rangemax(1),rangemin(2):rangemax(2),:));
ims_ch2 = single(qd2_trans(rangemin(1):rangemax(1),rangemin(2):rangemax(2),:));

ims_detect = ims_ch1 + ims_ch2;

%%  Set parameters
allcds = [];    
if setup.is_imgsz == 1
    thresh_dist = floor(boxsz / sqrt(2))-1;   %Two detected spots must be less than this distance
else 
    thresh_dist = 6;
end
%% Detect single molecules in both two channels 
%  Use uniform filter and maximum filter to obtain sub_region centers
display('Use uniform filter and maximum filter to obtain sub_region centers');

imsz = size(ims_detect, 1);
x=[];
y=[];
t=[];

% filter multiframe
sz = 3; 
[filteredim1] = unif_img(squeeze(ims_detect),sz);
sz = 9;
[filteredim2] = unif_img(squeeze(ims_detect),sz);
im_unif=filteredim1-filteredim2;

sz=4;
loc_max=(im_unif>=.999999999 * imdilate(im_unif, true(sz)));
im_max_L = loc_max & ((im_unif>thresh(1))&(im_unif<=thresh(2)));
im_max_H = loc_max & (im_unif>thresh(2));


centers_L = findcoord_seg(im_max_L);
centers_L(:,4) = 0;     %low threshold 
centers_H = findcoord_seg(im_max_H);
centers_H(:,4) = 1; 

allcds = cat(1,centers_L,centers_H);

%% Remove too close spots
allcds_keep = [];
fnum = size(ims_detect, 3);

for f = 0 : fnum - 1
    tmp_loc = allcds(allcds(:,3) == f, :);
    tmp_dist = pdist(tmp_loc);
    
    num_tmp = 1;
    tmp_remove = [];
    for ii = 1 : size(tmp_loc, 1) - 1
        for jj = ii + 1 : size(tmp_loc, 1)
            if (tmp_dist(num_tmp) < thresh_dist)
                tmp_remove = [tmp_remove ii jj];
            end
            num_tmp = num_tmp + 1;
        end
    end
    
    tmp_loc(tmp_remove,:) = [];
    allcds_keep = [allcds_keep; tmp_loc];
end


%%
boundmask=(allcds_keep(:,1)<=boxsz/2)|(allcds_keep(:,1)>=imsz-boxsz/2)|(allcds_keep(:,2)<=boxsz/2)|(allcds_keep(:,2)>=imsz-boxsz/2) | allcds_keep(:,4) == 0;
allcds_mask = allcds_keep(~boundmask,:);


%% Segmentation from two channels

% without channel2 transform
invertform = invert(tform);

allcds_mask_ch1 = [];

allcds_mask_ch1(:,1:2) = allcds_mask(:,1:2) + (rangemin - 1) + 1;  

[subregion_ch1 l1 t1]=cMakeSubregions(allcds_mask_ch1(:,2)-1,allcds_mask_ch1(:,1)-1,allcds_mask(:,3),boxsz,single(qd1_in)); %ims_ch1

allcds_affine = transformPointsForward(invertform,[allcds_mask_ch1(:,1) allcds_mask_ch1(:,2)]); %important (t1,l1)
allcds_affine_int = double(floor(allcds_affine));

allcds_affine_ch2 = allcds_affine_int - 1;
[subregion_ch2 l2 t2]=cMakeSubregions(allcds_affine_ch2(:,2),allcds_affine_ch2(:,1),...
    allcds_mask(:,3),boxsz, single(qd2_in)); %ims_ch2


%% Segmentation noise parameters
% sCMOS and EMCCD
if setup.is_sCMOS   %sCMOS case 
    % ch1
    [subvarim_ch1]=cMakeSubregions(allcds_mask_ch1(:,2)-1,allcds_mask_ch1(:,1)-1,...
        allcds_mask(:,3),boxsz,single(varim_ch1));
    [subgainim_ch1]=cMakeSubregions(allcds_mask_ch1(:,2)-1,allcds_mask_ch1(:,1)-1,...
        allcds_mask(:,3),boxsz,single(gainim_ch1));
    
    % ch2
    [subvarim_ch2]=cMakeSubregions(allcds_affine_ch2(:,2),allcds_affine_ch2(:,1),...
        allcds_mask(:,3),boxsz,single(varim_ch2));
    [subgainim_ch2]=cMakeSubregions(allcds_affine_ch2(:,2),allcds_affine_ch2(:,1),...
        allcds_mask(:,3),boxsz,single(gainim_ch2));
    
    % output var
    subvar_ch1 = subvarim_ch1 ./subgainim_ch1 ./subgainim_ch1;
    subvar_ch2 = subvarim_ch2 ./subgainim_ch2 ./subgainim_ch2;
else    %EMCCD case 
    subvar_ch1 = subregion_ch1 .* 0;
    subvar_ch2 = subregion_ch2 .* 0;
end
%% output

frame = allcds_mask(:,3);
l = l1 + rangemin(1) - 1;
t = t1 + rangemin(2) - 1;
offset_seg = allcds_affine(:,1:2) - allcds_affine_int(:,1:2);   %offset non-integer, output

end

