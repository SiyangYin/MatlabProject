function [] = sxdbcl22AM(sxdbcl_22)
figure;
t=83:109;
sxdbcl_22_1=sxdbcl_22(:,:,:,1);
sxdbcl_22_2=sxdbcl_22(:,:,:,2);
sxdbcl_22_3=sxdbcl_22(:,:,:,3);
sxdbcl_22_4=sxdbcl_22(:,:,:,4);
sxdbcl_22_5=sxdbcl_22(:,:,:,5);
sxdbcl_22_6=sxdbcl_22(:,:,:,6);
sxdbcl_22_7=sxdbcl_22(:,:,:,7);
sxdbcl_22_8=sxdbcl_22(:,:,:,8);
sxdbcl_22_9=sxdbcl_22(:,:,:,9);
sxdbcl_22_10=sxdbcl_22(:,:,:,10);
sxdbcl_22_11=sxdbcl_22(:,:,:,11);
sxdbcl_22_12=sxdbcl_22(:,:,:,12);
sxdbcl_22_13=sxdbcl_22(:,:,:,13);
sxdbcl_22_14=sxdbcl_22(:,:,:,14);
sxdbcl_22_15=sxdbcl_22(:,:,:,15);
sxdbcl_22_16=sxdbcl_22(:,:,:,16);
sxdbcl_22_17=sxdbcl_22(:,:,:,17);
sxdbcl_22_18=sxdbcl_22(:,:,:,18);
sxdbcl_22_19=sxdbcl_22(:,:,:,19);
sxdbcl_22_20=sxdbcl_22(:,:,:,20);
sxdbcl_22_21=sxdbcl_22(:,:,:,21);
sxdbcl_22_22=sxdbcl_22(:,:,:,22);
sxdbcl_22_23=sxdbcl_22(:,:,:,23);
sxdbcl_22_24=sxdbcl_22(:,:,:,24);
sxdbcl_22_25=sxdbcl_22(:,:,:,25);
sxdbcl_22_26=sxdbcl_22(:,:,:,26);
sxdbcl_22_27=sxdbcl_22(:,:,:,27);
sxdbcl_22_(1)=mean(sxdbcl_22_1(sxdbcl_22_1~=0));
sxdbcl_22_(2)=mean(sxdbcl_22_2(sxdbcl_22_2~=0));
sxdbcl_22_(3)=mean(sxdbcl_22_3(sxdbcl_22_3~=0));
sxdbcl_22_(4)=mean(sxdbcl_22_4(sxdbcl_22_4~=0));
sxdbcl_22_(5)=mean(sxdbcl_22_5(sxdbcl_22_5~=0));
sxdbcl_22_(6)=mean(sxdbcl_22_6(sxdbcl_22_6~=0));
sxdbcl_22_(7)=mean(sxdbcl_22_7(sxdbcl_22_7~=0));
sxdbcl_22_(8)=mean(sxdbcl_22_8(sxdbcl_22_8~=0));
sxdbcl_22_(9)=mean(sxdbcl_22_9(sxdbcl_22_9~=0));
sxdbcl_22_(10)=mean(sxdbcl_22_10(sxdbcl_22_10~=0));
sxdbcl_22_(11)=mean(sxdbcl_22_11(sxdbcl_22_11~=0));
sxdbcl_22_(12)=mean(sxdbcl_22_12(sxdbcl_22_12~=0));
sxdbcl_22_(13)=mean(sxdbcl_22_13(sxdbcl_22_13~=0));
sxdbcl_22_(14)=mean(sxdbcl_22_14(sxdbcl_22_14~=0));
sxdbcl_22_(15)=mean(sxdbcl_22_15(sxdbcl_22_15~=0));
sxdbcl_22_(16)=mean(sxdbcl_22_16(sxdbcl_22_16~=0));
sxdbcl_22_(17)=mean(sxdbcl_22_17(sxdbcl_22_17~=0));
sxdbcl_22_(18)=mean(sxdbcl_22_18(sxdbcl_22_18~=0));
sxdbcl_22_(19)=mean(sxdbcl_22_19(sxdbcl_22_19~=0));
sxdbcl_22_(20)=mean(sxdbcl_22_20(sxdbcl_22_20~=0));
sxdbcl_22_(21)=mean(sxdbcl_22_21(sxdbcl_22_21~=0));
sxdbcl_22_(22)=mean(sxdbcl_22_22(sxdbcl_22_22~=0));
sxdbcl_22_(23)=mean(sxdbcl_22_23(sxdbcl_22_23~=0));
sxdbcl_22_(24)=mean(sxdbcl_22_24(sxdbcl_22_24~=0));
sxdbcl_22_(25)=mean(sxdbcl_22_25(sxdbcl_22_25~=0));
sxdbcl_22_(26)=mean(sxdbcl_22_26(sxdbcl_22_26~=0));
sxdbcl_22_(27)=mean(sxdbcl_22_27(sxdbcl_22_27~=0));
plot(t,sxdbcl_22_,'k-*','LineWidth',2);
set(gca,'XTick',83:109);
axis([83 109 200 240]);
legend('SW downwelling for 100% overcast sky at surface (w/m^2)','Location','NorthWest');
title({'Annual mean time series of SW downwelling for 100% overcast sky at surface (w/m^2)','in Montreal(45°30′N，73°40′W) from 1983 to 2009'});
xlabel('Years from 1983 to 2009 (minus 1900)');
ylabel('SW downwelling for 100% overcast sky at surface (w/m^2)');