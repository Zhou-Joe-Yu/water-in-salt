%%%This code is used to calculate the ionic cluster distribution. Start from one cation and search the neighboring oxygen atom in the counterion. 
%%%Then start from the all the four oxygen atoms from the counterion as the center and search the neighboring cation. Repeat this process.
clear all
close all

atom_num = 2000;
cut_off  = 0.3; %first valley rdf from Li to O in TFSI
frame    = 251;
box_size = 8.62149;

pos_li=dlmread('pos_li.xvg','',22+atom_num*3,1);
pos_os=dlmread('pos_os.xvg','',22+atom_num*4*3,1);

%frame = 2; %%%

for t=1:frame
    for size=1:atom_num*2
        cluster_statistic{t,size}=[];
    end
    dist=zeros(atom_num,atom_num*4);
    count=1;
    for n=1:atom_num
        temp_li=[pos_li(t,(n-1)*3+1), pos_li(t,(n-1)*3+2), pos_li(t,(n-1)*3+3)];
        for m=1:atom_num
            
            temp_os1=[pos_os(t,(m-1)*12+1), pos_os(t,(m-1)*12+2), pos_os(t,(m-1)*12+3)];
            temp_os2=[pos_os(t,(m-1)*12+1+3), pos_os(t,(m-1)*12+2+3), pos_os(t,(m-1)*12+3+3)];
            temp_os3=[pos_os(t,(m-1)*12+1+6), pos_os(t,(m-1)*12+2+6), pos_os(t,(m-1)*12+3+6)];
            temp_os4=[pos_os(t,(m-1)*12+1+9), pos_os(t,(m-1)*12+2+9), pos_os(t,(m-1)*12+3+9)];
            
            dx1=abs(temp_li(1)-temp_os1(1));dy1=abs(temp_li(2)-temp_os1(2));dz1=abs(temp_li(3)-temp_os1(3));
            dx2=abs(temp_li(1)-temp_os2(1));dy2=abs(temp_li(2)-temp_os2(2));dz2=abs(temp_li(3)-temp_os2(3));
            dx3=abs(temp_li(1)-temp_os3(1));dy3=abs(temp_li(2)-temp_os3(2));dz3=abs(temp_li(3)-temp_os3(3));
            dx4=abs(temp_li(1)-temp_os4(1));dy4=abs(temp_li(2)-temp_os4(2));dz4=abs(temp_li(3)-temp_os4(3));
            
            
            dist_temp1=sqrt((min(dx1,box_size-dx1))^2+(min(dy1,box_size-dy1))^2+(min(dz1,box_size-dz1))^2);
            dist_temp2=sqrt((min(dx2,box_size-dx2))^2+(min(dy2,box_size-dy2))^2+(min(dz2,box_size-dz2))^2);
            dist_temp3=sqrt((min(dx3,box_size-dx3))^2+(min(dy3,box_size-dy3))^2+(min(dz3,box_size-dz3))^2);
            dist_temp4=sqrt((min(dx4,box_size-dx4))^2+(min(dy4,box_size-dy4))^2+(min(dz4,box_size-dz4))^2);
            %%%get the matrix of a-b distance;
            dist(n,(m-1)*4+1:m*4)=[dist_temp1,dist_temp2,dist_temp3,dist_temp4];
            %%%get the index of ion pairs
            if find([dist_temp1,dist_temp2,dist_temp3,dist_temp4]<=cut_off)~=0
               PX(count)=n;
               PY(count)=m;
               count=count+1;
            end
        end
    end
    
    PP1=[PX' PY']; %% ion pair index A-B
    
    if isempty(PP1)==1
       %%%ssip
       cluster_statistic{t,1}=[atom_num atom_num];
       continue %% repeat to next time step
    end
    
    PP2=[PP1(:,2) PP1(:,1)];
    PP2=sortrows(PP2,1);
    
    unique_a=unique(PP1(:,1));%%get the index of unique A 
    %unique_b=unique(PP2(:,1));%%get the index of unique B
    
    unique_a1=unique_a;
    
    for i=1:length(unique_a) %% level 1 a index
        
        if sum(unique_a1)~=0 && unique_a(i)==unique_a1(1) 
        
        index_B1=PP1(find(PP1(:,1)==unique_a(i)),2); %% level 1 B index 
        
        for j=1:length(index_B1)
            
            temp_index_a2{j}=PP2(find(PP2(:,1)==index_B1(j)),2)';
        
        end
        
        index_A2=unique(cell2mat(temp_index_a2)); 
        
        index_A2_diff=setdiff(index_A2,unique_a(i)); %% level 2 a index
        
        index_A2_union=union(index_A2,unique_a(i));  %% level 1+2 a index
        
        %%% empty the temp_index_a2
        
        for j=1:length(index_B1)
            
            temp_index_a2{j}=[];
        
        end
        %%%
        
        if sum(index_A2_diff)==0
           
           cluster_statistic{t,1+length(index_B1)}=[cluster_statistic{t,1+length(index_B1)};[1 length(index_B1)]];
           
           %unique_a=setdiff(unique_a,unique_a(1));
           
           pp_n=1; %% only two levels
        
        else
            
           pp_n=0; %% need to find in next level
        
        end  
        
        while pp_n==0;
            
            for m=1:length(index_A2_diff)
                
                temp_index_b2{m}=PP1(find(PP1(:,1)==index_A2_diff(m)),2)';
                
            end    
                
            index_B2=unique(cell2mat(temp_index_b2)); 
            
            index_B2_diff=setdiff(index_B2,index_B1); %% level 2 b index
            
            index_B2_union=union(index_B2,index_B1);  %% leval 1+2 b index
            
            %%% empty the temp_index_B2
            for m=1:length(index_A2_diff)
                
                temp_index_b2{m}=[];
                
            end   
            %%%
            
            if sum(index_B2_diff)==0
                
               cluster_statistic{t,length(index_A2_union)+length(index_B1)}=[cluster_statistic{t,length(index_A2_union)+length(index_B1)};[length(index_A2_union) length(index_B1)]]; 
                
               %unique_a1=setdiff(unique_a,index_A2_union);
               
               break %% jump out of the while
            
            end
            
            for n=1:length(index_B2_diff) 
                
                temp_index_a3{n}=PP2(find(PP2(:,1)==index_B2_diff(n)),2)'; 
               
            end    
            
            index_A3=unique(cell2mat(temp_index_a3));  
            
            index_A3_diff=setdiff(index_A3,index_A2_union); %%level 3 a index
            
            index_A3_union=union(index_A3,index_A2_union);  %%level 1+2+3 a index
            
            %%% empty the temp_index_a3
            for n=1:length(index_B2_diff)
                
                temp_index_a3{n}=[]; 
                
            end   
            %%%
            
            if sum(index_A3_diff)==0
                
               cluster_statistic{t,length(index_A2_union)+length(index_B2_union)}=[cluster_statistic{t,length(index_A2_union)+length(index_B2_union)};[length(index_A2_union) length(index_B2_union)]]; 
                
               %unique_a1=setdiff(unique_a,index_A2_union);
               
               pp_n=1;
            
            else
               
               %cluster_statistic{t,length(index_A3_union)+length(index_B2_union)}=[cluster_statistic{t,length(index_A3_union)+length(index_B2_union)};[length(index_A3_union) length(index_B2_union)]];  
                
               %unique_a1=setdiff(unique_a,index_A3_union);
               
               pp_n=0;
            
            end
            
            index_A2_diff=index_A3_diff;
            
            index_A2_union=index_A3_union;
            
            index_B1=index_B2_union;
            
        end
        
        unique_a1=setdiff(unique_a1,index_A2_union);
        
        end
        
    end  
    
    fprintf('%d\n', t)
    %%% empty
%     PP1=[];
%     PP2=[];
%     unique_a=[];
%     unique_a1=[];
%     index_B1=[];
%     index_A2=[];
%     index_A2_diff=[];
%     index_A2_union=[];
%     index_B2=[];
%     index_B2_diff=[];
%     index_B2_union=[];
%     index_A3=[];
%     index_A3_diff=[];
%     index_A3_union=[];
    %%%
    
end    

for n=1:atom_num*2
    
    statistic{n}=[];
    
    for t=1:frame

        if isempty(cluster_statistic{t,n})~=1

           statistic{n}=[statistic{n};cluster_statistic{t,n}]; 

        end
    end    
end    

matrix=zeros(atom_num,atom_num);
for n=1:atom_num*2
    if sum(statistic{n})~=0
        for m=1:length(statistic{n}(:,1))
            x=statistic{n}(m,1);
            y=statistic{n}(m,2);
            matrix(x,y)=matrix(x,y)+1;
        end
    end
end

contourf([1:atom_num],[1:atom_num],matrix,10);
%%%x axis is tfsi
%%%y axis is Li

save matrix.dat matrix -ascii
    
cluster_large=0;

for n=3:atom_num*2
    if sum(statistic{n})~=0
        cluster_large=cluster_large+length(statistic{n}(:,1))*n;
    end
end     
    
