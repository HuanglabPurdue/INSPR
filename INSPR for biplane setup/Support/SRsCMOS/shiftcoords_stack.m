function [out]=shiftcoords_stack(xin,shiftx)

segnum=numel(xin);
xout{1}=xin{1};
for ii=2:1:segnum
    xout{ii}=xin{ii}-sum(shiftx(1:ii-1,1));
end
out=xout;