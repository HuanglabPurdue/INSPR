function strout = catstruct(str1,str2,dim)

fd = fields(str1);
for ii = 1:numel(fd)
    strout.(fd{ii}) = cat(dim,str1.(fd{ii}), str2.(fd{ii}));
end

end