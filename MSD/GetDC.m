function P = GetDC(a,start,endp)

judge = a(2,1) - a(1,1);
index_s = find(abs(a(:,1)-start) < judge/2);
index_e = find(abs(a(:,1)-endp) < judge/2);

b = polyfit(a(index_s:index_e,1),a(index_s:index_e,2),1);
P = b(1);

end
