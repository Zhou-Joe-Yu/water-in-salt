clear all
close all

load('time_water.dat');
load('time_ionic.dat')
load('msd_ionic.dat');
load('msd_water.dat');
load('../li_x.dat');
load('../li_y.dat');

plot(time_ionic,msd_ionic,'-r','linewidth',2);hold on
%plot(time_water,msd_water,'-b','linewidth',2);
plot(li_x,li_y,'-g','linewidth',2);
legend('Aggregates','Water','Total');
axis([0 1000 0 1.0]);

set(gca,'fontsize',15,'fontweight','bold');
set(gca,'linewidth',1);
set(get(gca,'xlabel'),'FontSize', 15, 'FontWeight', 'Bold','Fontname','Times New Roman');
set(get(gca,'ylabel'),'FontSize', 15, 'FontWeight', 'Bold','Fontname','Times New Roman');
print('-dpng','-r300','divide_msd.png');




