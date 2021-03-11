%%build index 
clear all
close all

initial_index=2000;
ion_number=2000;
atom_in_ion=15;

fid=fopen('index_tfsi.ndx','w');

for i=1:ion_number
    
    fprintf(fid,'%s %d %s\n','[',i,']');
    
    for j=1:atom_in_ion
        fprintf(fid,'%d ',initial_index+atom_in_ion*(i-1)+j);
    end    
    
    fprintf(fid,'\n');

end

fclose(fid);

