clear all
close all

load('matrix.dat');

matrix=matrix./251;

v=log10(min(min(matrix(matrix~=0)))):0.5:log10(max(max(matrix)));
v_new=-2.5:0.5:2.5;
label=10.^v_new;

x=1:2000;
y=1:2000;
z=log10(matrix);
[h h]=contourf(x,y,z,v)

set(h,'LineStyle','none');


%colorbar;
colorbar('FontSize',15,'YTick',v(2:2:10),'YtickLabel',label(2:2:10));hold on

axis equal
axis([1 150 1 150])

set(gca,'xtick',[1 30 60 90 120 150]);
set(gca,'ytick',[1 30 60 90 120 150]);

set(gca,'fontsize',15,'fontweight','bold');
set(gca,'linewidth',1);
set(get(gca,'xlabel'),'FontSize', 15, 'FontWeight', 'Bold','Fontname','Times New Roman');
set(get(gca,'ylabel'),'FontSize', 15, 'FontWeight', 'Bold','Fontname','Times New Roman');


%print('-dtiff','-r300','cluster_analysis.tif');

[xx yy]=meshgrid(x,y);
xx_r=xx(:);
yy_r=yy(:);
m_r=matrix(:);
ft = fittype('poly1');
%ft = fittype({'x','0'});
cf = fit(xx_r,yy_r,ft,'weight',m_r);
test=plot(cf,'fit',0.95);
