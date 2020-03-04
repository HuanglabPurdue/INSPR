


function allcds = findcoord_seg(im_max)

imsz = size(im_max, 1);

a = find(im_max);      %potential segmentation centers                               
z = floor(a/imsz/imsz);
pnum = mod(a,imsz*imsz);
y = mod(pnum,imsz);
x = floor(pnum/imsz);
allcds =[x y z];