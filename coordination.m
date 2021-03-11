clear all
close all

atom_LI=2000;
atom_OS=8000;
atom_OW=5600;
frame=251;
cutoff=0.28;
box_size=8.62149;

pos_li=dlmread('pos_li.xvg','',22+atom_LI*3,1);
pos_os=dlmread('pos_os.xvg','',22+atom_OS*3,1);
pos_ow=dlmread('pos_ow.xvg','',22+atom_OW*3,1);

count_OS=1;
count_OW=1;

Coord=zeros(10,10);
tfsi_num=zeros(1,10);
tfsi_num1=zeros(1,10);
tfsi_num2=zeros(1,10);
tfsi_num3=zeros(1,10);

for t=1:frame
    for i=1:atom_LI
        index_os=zeros(1,atom_OS);
        for j=1:atom_OS
            dx1=abs(pos_li(t,(i-1)*3+1)-pos_os(t,(j-1)*3+1));
            dy1=abs(pos_li(t,(i-1)*3+2)-pos_os(t,(j-1)*3+2));
            dz1=abs(pos_li(t,(i-1)*3+3)-pos_os(t,(j-1)*3+3));
            dist_temp1=sqrt((min(dx1,box_size-dx1))^2+(min(dy1,box_size-dy1))^2+(min(dz1,box_size-dz1))^2);
            if dist_temp1<=cutoff
                count_OS=count_OS+1;
                index_os(j)=1;
            end    
        end

        if sum(index_os)~=0
            tfsi_num(length(unique(ceil(find(index_os==1)/4))))=tfsi_num(length(unique(ceil(find(index_os==1)/4))))+1;
        end
        
        if sum(index_os)==2
            tfsi_num1(length(unique(ceil(find(index_os==1)/4))))=tfsi_num1(length(unique(ceil(find(index_os==1)/4))))+1;
        end  
        
        if sum(index_os)==3
            tfsi_num2(length(unique(ceil(find(index_os==1)/4))))=tfsi_num2(length(unique(ceil(find(index_os==1)/4))))+1;
        end  

        if sum(index_os)>1
            tfsi_num3(length(unique(ceil(find(index_os==1)/4))))=tfsi_num3(length(unique(ceil(find(index_os==1)/4))))+1;
        end  
        
        for k=1:atom_OW
            dx1=abs(pos_li(t,(i-1)*3+1)-pos_ow(t,(k-1)*3+1));
            dy1=abs(pos_li(t,(i-1)*3+2)-pos_ow(t,(k-1)*3+2));
            dz1=abs(pos_li(t,(i-1)*3+3)-pos_ow(t,(k-1)*3+3));
            dist_temp2=sqrt((min(dx1,box_size-dx1))^2+(min(dy1,box_size-dy1))^2+(min(dz1,box_size-dz1))^2);
            if dist_temp2<=cutoff
                count_OW=count_OW+1;
            end
        end
        
        Coord(count_OS,count_OW)=Coord(count_OS,count_OW)+1;
        count_OS=1;
        count_OW=1;
        
    end 
    sprintf('%d',t)
end    

contourf([0:9],[0:9],Coord,10);