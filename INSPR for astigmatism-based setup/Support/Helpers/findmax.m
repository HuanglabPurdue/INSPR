function [out rowind colind]=findmax(m)
% [out rowind colind]=findmax(m)
[tmp inds]=max(m,[],1);
[out colind]=max(tmp,[],2);
rowind=inds(colind);