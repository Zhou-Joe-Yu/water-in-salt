function a = readfile(name,skip)

fid = fopen(name,'rt'); 
indata = textscan(fid, '%f', 'HeaderLines',skip); 
fclose(fid); 
index = 1;
for i = 1: length(indata{1})
    if mod(i,2) == 1
        a(index,1) = indata{1}(i);
    end
    if mod(i,2) == 0
        a(index,2) = indata{1}(i);
        index = index + 1;
    end
end

end
