function [xout,yout,zout] = shiftcoords_all(x,y,z,cormask,shifts)

xout = x;
yout = y;
zout = z;

sz = size(shifts,1);

if ~isempty(x)
    for nn = 1:sz
    xout(cormask==nn) = x(cormask==nn)+shifts(nn,1);
    end
end

if ~isempty(y)
    for nn = 1:sz
    yout(cormask==nn) = y(cormask==nn)+shifts(nn,2);
    end
end

if ~isempty(z)
    for nn = 1:sz
    zout(cormask==nn) = z(cormask==nn)+shifts(nn,3);
    end
end

end
