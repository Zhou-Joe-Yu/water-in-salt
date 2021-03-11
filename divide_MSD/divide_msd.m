clear all
close all

load('index_li.mat');
frame=501;
ion_num=2000;
pos_li=dlmread('pos_li.xvg','',22+ion_num*3,1);
box_size=8.62149;
dt=20;%%ps

% %%exchange O and 1
% index_Li(find(index_Li==0))=2;
% index_Li(find(index_Li==1))=0;
% index_Li(find(index_Li==2))=1;

for m = 1:ion_num
    temp=index_Li(:,m);
    
    q = diff([0 temp' 0] == 1); 
    v1 = find(q == 1); %%%index of first Li
    v2 = find(q == -1) - find(q == 1); %%% length of each subset
    index_cell{m}=[v1' v2'];
end  

count=1;
for m = 1:ion_num
    
    for ii=1:length(index_cell{m}(:,1)) 
    
        if index_cell{m}(ii,2)>=2 %%only count the ion stays longer than x * dt
           temp_pos= pos_li(index_cell{m}(ii,1):(index_cell{m}(ii,1)+index_cell{m}(ii,2)-1),(3*(ion_num-1)+1:3*ion_num));
           temp_msd{count}= MSD(temp_pos,index_cell{m}(ii,2),box_size,dt);
           count=count+1;
        end
    
    end
    
end 

for jj=1:length(temp_msd)
    x(jj)=length(temp_msd{jj});
end

longest_length=max(x); %%find the longest time for Li stay in one condition

scaling=zeros(1,longest_length);
total_msd=zeros(1,longest_length);
for kk=1:length(temp_msd)
    
    total_msd=total_msd+[temp_msd{kk}(:,2)' zeros(1,longest_length-length(temp_msd{kk}(:,2)))];
    scaling=scaling+[ones(1,length(temp_msd{kk}(:,2))) zeros(1,longest_length-length(temp_msd{kk}(:,2)))];
     
end    

time=dt.*(1:longest_length);
scaling_msd=total_msd./scaling;

plot(time,scaling_msd);

% save time_water.dat time -ascii
% save msd_water.dat scaling_msd -ascii

save time_ionic.dat time -ascii
save msd_ionic.dat scaling_msd -ascii

