%% This function classifies single molecule to giving reference 
%  Z-position images, and average these Z-postion images in a certain
%  group.
%
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, October 2019
% 
function [ims_Zcal_ave_plane1, index_record_Zplanes] = classify_onePlane_par(subregion_ch1, ref_plane1, empupil) 


disp('Calculate the similarity between reference images and single molecules');
img_plane1 = single(subregion_ch1);    %normalization image

num_img = size(img_plane1,3);
num_ref = size(ref_plane1, 3);
imsz = empupil.imsz;  

similarity_in_Plane1 = zeros(num_ref, num_img);  %similarity

index_similarity = zeros(num_ref, num_img); %index
shift_row_plane1 = zeros(num_ref, num_img); %X,Y shift in two planes
shift_col_plane1 = zeros(num_ref, num_img);


parfor ii = 1 : num_img

    for jj = 1 : num_ref
        [shift1, shift2, tmpval1] = registration_in_each_channel(ref_plane1(:,:,jj), img_plane1(:,:,ii));
             
        if (abs(shift1) > 6 || abs(shift2) > 6)  %this value can be changed, now fixed
            continue;
        end
        
        %get X, Y shift in plane1
        shift_row_plane1(jj, ii) = shift1;
        shift_col_plane1(jj, ii) = shift2;
            
        if tmpval1 < empupil.min_similarity
            continue;
        end
     
        similarity_in_Plane1(jj, ii) = tmpval1;
    end    
end


%%
for ii = 1 : num_img
    % Sort the similarity
    [sort_similarity, index_sort] =  sort(similarity_in_Plane1(:,ii),'descend'); 
    if (sort_similarity(1) == 0)
        continue;
    end
    % Determine single molecule image belong which reference image
    index_similarity(index_sort(1), ii) = 1;
    for jj = 2 : num_ref
        if (sort_similarity(jj) >= sort_similarity(1)-0.0) && (abs(index_sort(jj)-index_sort(1)) == 1)    %now fixed
            index_similarity(index_sort(jj),ii) = 1;
        else
            break;
        end
    end
end


%% Updata average images
disp('Update average images');

ims_Zcal_ave_plane1 = zeros(imsz,imsz,num_ref);
index_record_Zplanes = zeros(num_ref,1);
for ii = 1 : num_ref
    index_selection = find(index_similarity(ii,:) == 1);
    sz_index = size(index_selection,2);
    if sz_index > empupil.bin_lowerBound    %Edited by FX, each group must have enough spots 
        ims_plane1_shift = zeros(imsz,imsz,sz_index);
        for jj = 1 : sz_index
            ims_plane1_shift(:,:,jj) = FourierShift2D(similarity_in_Plane1(ii, index_selection(jj)) .* img_plane1(:,:,index_selection(jj)), [shift_row_plane1(ii, index_selection(jj)) shift_col_plane1(ii, index_selection(jj))]);            
        end
        % average the images
        index_record_Zplanes(ii) = 1;
        ims_Zcal_ave_plane1(:,:,ii) = mean(ims_plane1_shift,3);       
    end
end


%% Remove too far away Reassembled PSF

tmp_record_keep = find(index_record_Zplanes == 1);
size_tmp = length(tmp_record_keep);

for ii = ceil(size_tmp/2) : -1 : 2
    if tmp_record_keep(ii) > tmp_record_keep(ii-1)+2
        index_record_Zplanes(tmp_record_keep(ii-1)) = 0;
        tmp_record_keep(ii-1) = tmp_record_keep(ii);
    end
end

for ii = ceil(size_tmp/2) : size_tmp-1
    if tmp_record_keep(ii) < tmp_record_keep(ii+1)-2
        index_record_Zplanes(tmp_record_keep(ii+1)) = 0;
        tmp_record_keep(ii+1) = tmp_record_keep(ii);
    end
end 
