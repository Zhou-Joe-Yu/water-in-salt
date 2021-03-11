%%This code is used to calculate the structure factor from the radial distribution function among different atoms
clear all
close all

tot_atom=48800;
Li_num=2000;
TFSI_num=2000;
water_num=5600;
box=8.62149;
R=2.5;
delta_q=0.1;

neutral_factor=dlmread('neutral.dat','',1,0);
ion_factor=dlmread('ions.dat','',1,0);

for i=1:7

    for j=1:7

        P{i,j}=dlmread(['./rdf_results/rdf_',num2str(i,'%d'),num2str(j,'%d'),'.xvg'],'',24,0);
    
    end

end

b1=ion_factor(:,[1,3]); %Li
b2=neutral_factor(:,[1,10]);%F
b3=neutral_factor(:,[1,17]);%S
b4=neutral_factor(:,[1,9]);%O
b5=neutral_factor(:,[1,8]);%N
b6=neutral_factor(:,[1,7]);%C
b7=neutral_factor(:,[1,2]);%H

x1=Li_num/tot_atom;
x2=6*TFSI_num/tot_atom;
x3=2*TFSI_num/tot_atom;
x4=(4*TFSI_num+water_num)/tot_atom;
x5=TFSI_num/tot_atom;
x6=2*TFSI_num/tot_atom;
x7=2*water_num/tot_atom;

rho=tot_atom/box^3; %number/volume

b1_new=[[0:delta_q:20]',interp1(b1(1:56,1),b1(1:56,2),[0:delta_q/10:2])'];
b2_new=[[0:delta_q:20]',interp1(b2(1:56,1),b2(1:56,2),[0:delta_q/10:2])'];
b3_new=[[0:delta_q:20]',interp1(b3(1:56,1),b3(1:56,2),[0:delta_q/10:2])'];
b4_new=[[0:delta_q:20]',interp1(b4(1:56,1),b4(1:56,2),[0:delta_q/10:2])'];
b5_new=[[0:delta_q:20]',interp1(b5(1:56,1),b5(1:56,2),[0:delta_q/10:2])'];
b6_new=[[0:delta_q:20]',interp1(b6(1:56,1),b6(1:56,2),[0:delta_q/10:2])'];
b7_new=[[0:delta_q:20]',interp1(b7(1:56,1),b7(1:56,2),[0:delta_q/10:2])'];

denominator_temp1=x1.*b1_new(:,2)+x2.*b2_new(:,2)+x3.*b3_new(:,2)+x4.*b4_new(:,2)+x5.*b5_new(:,2)+x6.*b6_new(:,2)+x7.*b7_new(:,2);
denominator_temp2=denominator_temp1.*denominator_temp1;
denominator=[b1_new(:,1), denominator_temp2];

x_tot=[x1, x2, x3, x4, x5, x6, x7];
b_tot=[b1_new(:,2), b2_new(:,2), b3_new(:,2), b4_new(:,2), b5_new(:,2), b6_new(:,2), b7_new(:,2)];
S_tot=zeros(1,length(b1_new(:,1)));

for i=1:7

    for j=1:7        

        for q=0:delta_q:20
      
            index_q=round(q/delta_q)+1;
       	    
            for r=0.001:0.002:2.501
            
                index_r=floor(r/0.002)+1;
                
                temp(index_r)=4*pi*r^2*(P{i,j}(index_r,2)-1)*sin(q*r)/(q*r)*sin(pi*r/R)*R/(pi*r);
       	    
            end 
            
            S_test{i,j}(index_q)=rho*x_tot(i)*x_tot(j)*b_tot(index_q,i)*b_tot(index_q,j)*sum(temp.*0.002)/denominator(index_q,2);
        
        end
        
        S_tot=S_tot+S_test{i,j};
    
    end

end

xaxis=b1_new(:,1)';

figure(1)
plot(xaxis,S_tot,'-k','linewidth',2);hold on %all
plot(xaxis,S_test{2,6},':r','linewidth',2);hold on %F-C
plot(xaxis,S_test{2,2},':b','linewidth',2);hold on %F-F

legend('ALL','F-C','F-F','location','Southeast');

axis([0 20 -1.3 0.7]);

%xlabel('q (nm^{-1})');
set(gca,'ytick',[]);
set(gca,'xtick',[2 6 10 14 18]);
set(gca,'fontsize',20,'fontweight','bold');
set(gca,'linewidth',1);
set(get(gca,'xlabel'),'FontSize', 20, 'FontWeight', 'Bold','Fontname','Times New Roman');
set(get(gca,'ylabel'),'FontSize', 20, 'FontWeight', 'Bold','Fontname','Times New Roman');
title('MD 20m');
print('-dtiff','-r300','SAXS_20m.tif');

save x_01.dat xaxis -ascii
save y_01.dat S_tot -ascii

% figure(2) %ALL
% for i=1:7
%     for j=1:7
%         plot(xaxis,S_test{i,j});hold on
%     end
% end
