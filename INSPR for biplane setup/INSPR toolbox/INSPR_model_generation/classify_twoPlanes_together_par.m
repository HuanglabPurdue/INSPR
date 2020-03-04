% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, August 2019
%       
%
%% This function classifies single molecule of two channels to giving reference 
%  Z-position images, and average these Z-postion images in a certain
%  group.

function [ims_Zcal_ave_plane1, ims_Zcal_ave_plane2, index_record_Zplanes] = classify_twoPlanes_together_par(subregion_ch1, subregion_ch2, ...
    ref_plane1, ref_plane2, empupil) 


disp('Calculate the similarity between reference images and single molecules');
img_plane1 = double(subregion_ch1);   
img_plane2 = double(subregion_ch2);

num_img = size(img_plane1,3);
num_ref = size(ref_plane1, 3);
imsz = empupil.imsz;  

similarity_in_twoPlanes = zeros(num_ref, num_img);

index_similarity = zeros(num_ref, num_img); %index
shift_row_biplane = zeros(num_ref, num_img);
shift_col_biplane = zeros(num_ref, num_img);


parfor ii = 1 : num_img

    for jj = 1 : num_ref
        
        [shift1, shift2, tmpval] = registration_in_biplane(ref_plane1(:,:,jj), ref_plane2(:,:,jj), img_plane1(:,:,ii), img_plane2(:,:,ii));
        
        if (tmpval < empupil.min_similarity || abs(shift1) > 6 || abs(shift2) > 6)
            continue;
        end
        %get X, Y shift
        shift_row_biplane(jj, ii) = shift1;
        shift_col_biplane(jj, ii) = shift2;
        similarity_in_twoPlanes(jj, ii) = tmpval;     
            
    end
   
end

%%
for ii = 1 : num_img
    % Sort the similarity
    [sort_similarity, index_sort] =  sort(similarity_in_twoPlanes(:,ii),'descend'); 
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


%% Updata average images in two planes, num_ave_image = num_ref_image
disp('Update average images in two planes');

ims_Zcal_ave_plane1 = zeros(imsz,imsz,num_ref);
ims_Zcal_ave_plane2 = zeros(imsz,imsz,num_ref);
index_record_Zplanes = zeros(num_ref,1);
for ii = 1 : num_ref
    index_selection = find(index_similarity(ii,:) == 1);
    sz_index = size(index_selection,2);
    if sz_index > empupil.bin_lowerBound    %each group must have enough spots 
        ims_plane1_shift = zeros(imsz,imsz,sz_index);
        ims_plane2_shift = zeros(imsz,imsz,sz_index);
        for jj = 1 : sz_index
            ims_plane1_shift(:,:,jj) = FourierShift2D(img_plane1(:,:,index_selection(jj)), [shift_row_biplane(ii, index_selection(jj)) shift_col_biplane(ii, index_selection(jj))]);
            ims_plane2_shift(:,:,jj) = FourierShift2D(img_plane2(:,:,index_selection(jj)), [shift_row_biplane(ii, index_selection(jj)) shift_col_biplane(ii, index_selection(jj))]);
            
        end
        
        % average the images
        index_record_Zplanes(ii) = 1;
        ims_Zcal_ave_plane1(:,:,ii) = mean(ims_plane1_shift,3);       
        ims_Zcal_ave_plane2(:,:,ii) = mean(ims_plane2_shift,3);
    end
end


%% Remove too far away Reassembled PSF

tmp_record_keep = find(index_record_Zplanes == 1);
size_tmp = length(tmp_record_keep);

for ii = ceil(size_tmp/2) : -1 : 2
    if tmp_record_keep(ii) > tmp_record_keep(ii-1)+3
        index_record_Zplanes(tmp_record_keep(ii-1)) = 0;
        tmp_record_keep(ii-1) = tmp_record_keep(ii);
    end
end

for ii = ceil(size_tmp/2) : size_tmp-1
    if tmp_record_keep(ii) < tmp_record_keep(ii+1)-3
        index_record_Zplanes(tmp_record_keep(ii+1)) = 0;
        tmp_record_keep(ii+1) = tmp_record_keep(ii);
    end
end 
