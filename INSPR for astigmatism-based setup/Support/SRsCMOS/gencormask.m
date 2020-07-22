% Input:
%   t: vector of frame numbers, each single-molecule index in each image fame
%   nfs: number of frame per section
% Output:
%   cormask: section index. Single molecules of some sectiion will assign together

function [cormask] = gencormask(t,nfs)
fn = ceil(single(max(t))/nfs);
v = [0:fn];
cormask = zeros(size(t));
for ii = 1:numel(v)-1
    mask = (t >= v(ii)*nfs)&(t < v(ii+1)*nfs);
    cormask(mask) = ii;
end
cormask(cormask==0) = ii;
end