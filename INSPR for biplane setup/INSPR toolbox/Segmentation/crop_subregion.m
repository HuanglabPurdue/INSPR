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
%  input: qd1 and qd2
%  output: subregion_ch1 and subregion_ch2

%%
function [subregion_ch1,subregion_ch2,seg_display] = crop_subregion(qd1,qd2,tform,boxsz,thresh,thresh_dist,setup)


% if setup.is_sCMOS   %sCMOS case 
%     % sCMOS parameters
%     offsetim_ch1 = repmat(setup.sCMOS_input.ccdoffset_ch1,[1 1 size(qd1,3)]);
%     offsetim_ch2 = repmat(setup.sCMOS_input.ccdoffset_ch2,[1 1 size(qd2,3)]);
%     
%     gainim_ch1 = repmat(setup.sCMOS_input.gain_ch1,[1 1 size(qd1,3)]);
%     gainim_ch2 = repmat(setup.sCMOS_input.gain_ch2,[1 1 size(qd2,3)]);
%    
%     qd1_in = (qd1 - offsetim_ch1) ./ gainim_ch1;
%     qd2_in = (qd2 - offsetim_ch2) ./ gainim_ch2;
% else    %EMCCD case
%     qd1_in = (qd1 - setup.offset) /setup.gain;
%     qd2_in = (qd2 - setup.offset) /setup.gain;
% end
% 
% qd1_in(qd1_in<=0) = 1e-6;
% qd2_in(qd2_in<=0) = 1e-6;

qd1_in = qd1;
qd2_in = qd2;

%transfer qd2 to find center peak
qd2_trans = imwarp(qd2_in,tform,'cubic','OutputView',imref2d(size(qd2_in)));


imsz_original = size(qd1_in);

if setup.is_imgsz == 1
    rangemin=[10, 10];
    rangemax=[imsz_original(1)-10, imsz_original(1)-10];
else 
    rangemin=[1, 1];
    rangemax=[imsz_original(1), imsz_original(1)];
end

ims_ch1 = single(qd1_in(rangemin(1):rangemax(1),rangemin(2):rangemax(2),:));
ims_ch2 = single(qd2_trans(rangemin(1):rangemax(1),rangemin(2):rangemax(2),:));


ims_detect = ims_ch1 + ims_ch2;



%% Detect single molecules in both two channels 
%  Use uniform filter and maximum filter to obtain sub_region centers
display('Use uniform filter and maximum filter to obtain sub_region centers');
imsz = size(ims_detect, 1);
allcds = [];    %select region in both channel

x=[];
y=[];
t=[];

% filter multiframe
sz = 3; 
[filteredim1] = unif_img(squeeze(ims_detect),sz);
sz = 9;
[filteredim2] = unif_img(squeeze(ims_detect),sz);
im_unif=filteredim1-filteredim2;


sz=3;
loc_max=(im_unif>=.999999999 * imdilate(im_unif, true(sz)));
im_max_L = loc_max & ((im_unif>thresh(1))&(im_unif<=thresh(2)));
im_max_H = loc_max & (im_unif>thresh(2));


centers_L = findcoord_seg(im_max_L);
centers_L(:,4) = 0;     %low threshold 
centers_H = findcoord_seg(im_max_H);
centers_H(:,4) = 1;     %high threshold

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


%% Segmentation from two channels
[subregion_ch1 l1 t1]=cMakeSubregions(allcds_mask(:,2),allcds_mask(:,1),allcds_mask(:,3),boxsz, ims_ch1); %ims_ch1
[subregion_ch2 l2 t2]=cMakeSubregions(allcds_mask(:,2),allcds_mask(:,1),allcds_mask(:,3),boxsz, ims_ch2); %ims_ch2


%% save display
seg_display.ims_ch1 = ims_ch1;
seg_display.ims_ch2 = ims_ch2;
seg_display.allcds_mask = allcds_mask;
seg_display.t1 = t1;
seg_display.l1 = l1;
