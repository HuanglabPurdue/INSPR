function [out colind rowind zind]=findmax3d(m)
[tmp inds]=max(m,[],1);
[tmp2 inds2]=max(tmp,[],2);
[out zind]=max(tmp2,[],3);
colind=inds2(1,1,zind);
rowind=inds(1,colind,zind);