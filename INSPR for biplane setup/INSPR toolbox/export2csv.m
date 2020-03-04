%
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, August 2019
%
% Export to CVS format
% Output X, Y, Z postions
%
function [flag]=export2csv(pathname,filename,colornum,xcell,ycell,zcell,tcell)
Mtot = [];
flag = 0;


input_data = fullfile(pathname,filename);

% 
c =['Frame-ID,x [nm],y [nm],z [nm]'];
fid = fopen(input_data,'wt');
fprintf(fid,'%s\n',c);
fclose(fid);

for cc=1:1:colornum
    Mtottmp = [];
    x = xcell{cc};
    y = ycell{cc};
    z = zcell{cc};
    col1 = tcell{cc};

    col2 = x;
    col3 = y;
    col4 = z;
  
    kk = 4;

    %% cat them together
    
    for ii = 1:1:kk
        eval(['Mtottmp=cat(2,Mtottmp,' ['col' num2str(ii)] ');']);
    end
    
    if cc == 1
        Mtot = Mtottmp;
    else
        Mtot = cat(1,Mtot,Mtottmp);
    end
end

try
    dlmwrite(input_data,Mtot,'-append')
    flag = 1;
catch
    error('Error in writing csv file, please check write permission in current folder');
end

