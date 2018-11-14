function [] = str16_17GD(txdwbt_16,txupbt_17,iy_,im_,id_,ih)
%STR16_17GD Summary of this function goes here
%   Detailed explanation goes here
lon_=1.25:2.5:358.75;
lat_=-88.75:2.5:88.75;
[lon_,lat_]=meshgrid(lon_,lat_);
figure;
contourf(lon_,lat_,(txdwbt_16-txupbt_17)');
colorbar;
%colorbar('Direction','Reverse');
hold on;
coast=load('coast');
idx=coast.long<0;
coast.long(idx)=coast.long(idx)+360;
plot(coast.long,coast.lat,'k');
axis([0 360 -90 90]);
grid on;
title({'Global distribution of Surface net thermal radiation for full sky (w/m^2)',strcat('in',32,im_,32,id_,32,num2str(ih),32,'hour',32,iy_,32,'year from ISCCP')});
xlabel('Longitude');
ylabel('Latitude');
end

