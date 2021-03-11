clear all;
close all;

n = 5;

for i = 1 : n
    filename = ['msd_o' num2str(i) '.xvg'];
    %filename = ['msd_lat_C1' num2str(i) '.xvg'];
    EL{i} = readfile(filename,19);
end

for i = 1: n
    start = 0.3 * EL{i}(end,1);
    endp = 0.9 * EL{i}(end,1);
    D_EL(i) = GetDC(EL{i},start,endp)/6;
end

DC_EL = mean(D_EL); %%%unit is sigma2/tau
error_EL = std(D_EL);

figure;
for i = 1: n
    loglog(EL{i}(:,1),EL{i}(:,2),'b-'); hold on   
end

EL_yfirst=(EL{1}(:,2)+EL{2}(:,2)+EL{3}(:,2))/3;
EL_xfirst=EL{1}(:,1);
save poly_na_71_x.dat EL_xfirst -ascii
save poly_na_71_y.dat EL_yfirst -ascii

%xlabel('Time');
%ylabel('polymer com MSD');
%legend('polymer com');
%axis([1 3*10^5 10^-4 10]);
set(gca,'fontsize',15,'fontweight','bold');
set(gca,'linewidth',1);
set(get(gca,'xlabel'),'FontSize', 15, 'FontWeight', 'Bold','Fontname','Times New Roman');
set(get(gca,'ylabel'),'FontSize', 15, 'FontWeight', 'Bold','Fontname','Times New Roman');
%print('-dpng','-r300','msd_poly_com_71.png');

