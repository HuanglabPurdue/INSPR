
function [out]=shiftcoords_LS(xout,shiftx,toutin,frmnum,reverseflag,currcy)

if nargin<=4
    segn=unique(currcy);

    for ii=2:1:numel(segn)
        maskt=(currcy==segn(ii));
        xout(maskt)=xout(maskt)-sum(shiftx(1:ii-1,1));
    end
    out=xout;
elseif reverseflag==0
        segn=unique(currcy);

    for ii=2:1:numel(segn)
        maskt=(currcy==segn(ii));
        xout(maskt)=xout(maskt)-sum(shiftx(1:ii-1,1));
    end
    out=xout;
else
    segn=unique(currcy);

    for ii=segn(1):1:segn(end)-1
        maskt=(currcy==ii);
        xout(maskt)=xout(maskt)+sum(shiftx(ii:end-1,1));
    end
    out=xout;
end
    