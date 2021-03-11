%%%this code is used to generate the index of Li, 1 represents Li is
%%%coordinated by TFSI, O represents Li is only coordinated by water
clear all
close all

atom_LI=2000;
atom_OS=atom_LI*4;
frame=501; %2ns~12ns
cutoff=0.28;
box_size=8.62149;

pos_li=dlmread('pos_li.xvg','',22+atom_LI*3,1);
pos_os=dlmread('pos_os.xvg','',22+atom_OS*3,1);

index_Li=zeros(frame,atom_LI);

for t=1:frame
    for i=1:atom_LI
        dist_temp1=zeros(1,atom_OS);
        for j=1:atom_OS
	        dx1=abs(pos_li(t,(i-1)*3+1)-pos_os(t,(j-1)*3+1));
            dy1=abs(pos_li(t,(i-1)*3+2)-pos_os(t,(j-1)*3+2));
            dz1=abs(pos_li(t,(i-1)*3+3)-pos_os(t,(j-1)*3+3));
            dist_temp1(j)=sqrt((min(dx1,box_size-dx1))^2+(min(dy1,box_size-dy1))^2+(min(dz1,box_size-dz1))^2);
%             if dist_temp1<=cutoff
%                 %sprintf('%s','yes')
%                 index_Li(t,i)=1;
%                 continue;
%             end
        end
        if min(dist_temp1)<=cutoff
           index_Li(t,i)=1;
%           continue;
        end  
    end
    sprintf('%d',t)
end

save index_li.mat index_Li -mat