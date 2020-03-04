
function unifim=unif_img(im,sz)
im = single(im);

kerim = ones(sz,sz) / (sz*sz);
unifim = imfilter(im, kerim, 'conv','replicate');