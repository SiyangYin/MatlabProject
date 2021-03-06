function [] = txupbt17GD(fluxrp,iy_,im_,id_,ih)
%[boxnmb,fluxrp]=eq25ea(0,eqrout(:,17),fluxrp);
%iy_=99; im_=07; id_=15; ih=21;
lon_=1.25:2.5:358.75;
lat_=-88.75:2.5:88.75;
[lon_,lat_]=meshgrid(lon_,lat_);
figure;
contourf(lon_,lat_,fluxrp');
colorbar;
%colorbar('Direction','Reverse');
hold on;
coast=load('coast');
idx=coast.long<0;
coast.long(idx)=coast.long(idx)+360;
plot(coast.long,coast.lat,'k');
axis([0 360 -90 90]);
grid on;
title({'Global distribution of LW upwelling for full sky at surface (w/m^2)',strcat('in',32,im_,32,id_,32,num2str(ih),'hour',32,iy_,'year')});
xlabel('Longitude');
ylabel('Latitude');