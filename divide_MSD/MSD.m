function temp_msd = MSD(POS,LENGTH,box_size,dt)

distance=zeros(1,LENGTH-1);
count=zeros(1,LENGTH-1);
for i=1:LENGTH-1 %%initial
    for m=1:LENGTH-i %% delta t
        x_initial=POS(i,1);
        y_initial=POS(i,2);
        z_initial=POS(i,3);
        x_final=POS(i+m,1);
        y_final=POS(i+m,2);
        z_final=POS(i+m,3);

        %%PBC
        dx=abs(x_initial-x_final);
        dy=abs(y_initial-y_final);
        dz=abs(z_initial-z_final);
        distance(m)=distance(m)+((min(dx,box_size-dx))^2+(min(dy,box_size-dy))^2+(min(dz,box_size-dz))^2);
        count(m)=count(m)+1;
    end
end

distance = distance./count;
track=1:(LENGTH-1);
t=dt.*track;
temp_msd=[t' distance'];
