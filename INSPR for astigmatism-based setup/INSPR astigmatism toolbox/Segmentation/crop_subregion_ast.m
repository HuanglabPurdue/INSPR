%
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, July 2020
%
%% Script for segmentation from single molecule dataset
%  input: ims
%  output: subregion_ch1

%%
function [subregion_ch1,seg_display] = crop_subregion_ast(data1,boxsz,thresh,thresh_dist,setup)

% if setup.is_sCMOS   %sCMOS case 
%     % sCMOS parameters
%     offsetim_ch1 = repmat(setup.sCMOS_input.ccdoffset_ch1,[1 1 size(data1,3)]);    
%     gainim_ch1 = repmat(setup.sCMOS_input.gain_ch1,[1 1 size(data1,3)]);   
%     data1_in = (data1 - offsetim_ch1) ./ gainim_ch1;
% else    %EMCCD case
%     data1_in = (data1 - setup.offset) /setup.gain;
% end
% 
% data1_in(data1_in<=0) = 1e-6;
data1_in = data1;

imsz_original = size(data1_in);
if setup.is_imgsz == 1
    rangemin=[5, 5];
    rangemax=[imsz_original(1)-5, imsz_original(1)-5];
else 
    rangemin=[1, 1];
    rangemax=[imsz_original(1), imsz_original(1)];
end

ims_ch1 = single(data1_in(rangemin(1):rangemax(1),rangemin(2):rangemax(2),:));
ims_detect = ims_ch1;


%% Detect single molecules
%  Use uniform filter and maximum filter to obtain sub_region centers
display('Use uniform filter and maximum filter to obtain sub_region centers');
imsz = size(ims_detect, 1);
allcds = [];    

x=[];
y=[];
t=[];

% filter 
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
centers_H(:,4) = 1;   %high threshold

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

%% selected spots
boundmask=(allcds_keep(:,1)<=boxsz/2)|(allcds_keep(:,1)>=imsz-boxsz/2)|(allcds_keep(:,2)<=boxsz/2)|(allcds_keep(:,2)>=imsz-boxsz/2) | allcds_keep(:,4) == 0;
allcds_mask = allcds_keep(~boundmask,:);


%% Segmentation sub-region of single molecule
[subregion_ch1 l1 t1]=cMakeSubregions(allcds_mask(:,2),allcds_mask(:,1),allcds_mask(:,3),boxsz, ims_ch1); %ims_ch1


%% save display
seg_display.ims_ch1 = ims_ch1;
seg_display.allcds_mask = allcds_mask;
seg_display.t1 = t1;
seg_display.l1 = l1;


