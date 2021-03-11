clear all
close all

fid=fopen('index_tfsi.ndx','w');

for i=1:2000
    
    fprintf(fid,'%s %d %s\n','[',i,']');
    
    for j=1:15
        fprintf(fid,'%d ',2000+15*(i-1)+j);
    end    
    
    fprintf(fid,'\n');

end

fclose(fid);

