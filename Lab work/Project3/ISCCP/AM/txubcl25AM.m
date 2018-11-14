function [] = txubcl25AM(txubcl_25)
figure;
t=83:109;
txubcl_25_1=txubcl_25(:,:,:,1);
txubcl_25_2=txubcl_25(:,:,:,2);
txubcl_25_3=txubcl_25(:,:,:,3);
txubcl_25_4=txubcl_25(:,:,:,4);
txubcl_25_5=txubcl_25(:,:,:,5);
txubcl_25_6=txubcl_25(:,:,:,6);
txubcl_25_7=txubcl_25(:,:,:,7);
txubcl_25_8=txubcl_25(:,:,:,8);
txubcl_25_9=txubcl_25(:,:,:,9);
txubcl_25_10=txubcl_25(:,:,:,10);
txubcl_25_11=txubcl_25(:,:,:,11);
txubcl_25_12=txubcl_25(:,:,:,12);
txubcl_25_13=txubcl_25(:,:,:,13);
txubcl_25_14=txubcl_25(:,:,:,14);
txubcl_25_15=txubcl_25(:,:,:,15);
txubcl_25_16=txubcl_25(:,:,:,16);
txubcl_25_17=txubcl_25(:,:,:,17);
txubcl_25_18=txubcl_25(:,:,:,18);
txubcl_25_19=txubcl_25(:,:,:,19);
txubcl_25_20=txubcl_25(:,:,:,20);
txubcl_25_21=txubcl_25(:,:,:,21);
txubcl_25_22=txubcl_25(:,:,:,22);
txubcl_25_23=txubcl_25(:,:,:,23);
txubcl_25_24=txubcl_25(:,:,:,24);
txubcl_25_25=txubcl_25(:,:,:,25);
txubcl_25_26=txubcl_25(:,:,:,26);
txubcl_25_27=txubcl_25(:,:,:,27);
txubcl_25_(1)=mean(txubcl_25_1(txubcl_25_1~=0));
txubcl_25_(2)=mean(txubcl_25_2(txubcl_25_2~=0));
txubcl_25_(3)=mean(txubcl_25_3(txubcl_25_3~=0));
txubcl_25_(4)=mean(txubcl_25_4(txubcl_25_4~=0));
txubcl_25_(5)=mean(txubcl_25_5(txubcl_25_5~=0));
txubcl_25_(6)=mean(txubcl_25_6(txubcl_25_6~=0));
txubcl_25_(7)=mean(txubcl_25_7(txubcl_25_7~=0));
txubcl_25_(8)=mean(txubcl_25_8(txubcl_25_8~=0));
txubcl_25_(9)=mean(txubcl_25_9(txubcl_25_9~=0));
txubcl_25_(10)=mean(txubcl_25_10(txubcl_25_10~=0));
txubcl_25_(11)=mean(txubcl_25_11(txubcl_25_11~=0));
txubcl_25_(12)=mean(txubcl_25_12(txubcl_25_12~=0));
txubcl_25_(13)=mean(txubcl_25_13(txubcl_25_13~=0));
txubcl_25_(14)=mean(txubcl_25_14(txubcl_25_14~=0));
txubcl_25_(15)=mean(txubcl_25_15(txubcl_25_15~=0));
txubcl_25_(16)=mean(txubcl_25_16(txubcl_25_16~=0));
txubcl_25_(17)=mean(txubcl_25_17(txubcl_25_17~=0));
txubcl_25_(18)=mean(txubcl_25_18(txubcl_25_18~=0));
txubcl_25_(19)=mean(txubcl_25_19(txubcl_25_19~=0));
txubcl_25_(20)=mean(txubcl_25_20(txubcl_25_20~=0));
txubcl_25_(21)=mean(txubcl_25_21(txubcl_25_21~=0));
txubcl_25_(22)=mean(txubcl_25_22(txubcl_25_22~=0));
txubcl_25_(23)=mean(txubcl_25_23(txubcl_25_23~=0));
txubcl_25_(24)=mean(txubcl_25_24(txubcl_25_24~=0));
txubcl_25_(25)=mean(txubcl_25_25(txubcl_25_25~=0));
txubcl_25_(26)=mean(txubcl_25_26(txubcl_25_26~=0));
txubcl_25_(27)=mean(txubcl_25_27(txubcl_25_27~=0));
plot(t,txubcl_25_,'k-*','LineWidth',2);
set(gca,'XTick',83:109);
axis([83 109 320 360]);
legend('LW upwelling for 100% overcast sky at surface (w/m^2)','Location','NorthWest');
title({'Annual mean time series of LW upwelling for 100% overcast sky at surface (w/m^2)','in Montreal(45°30′N，73°40′W) from 1983 to 2009'});
xlabel('Years from 1983 to 2009 (minus 1900)');
ylabel('LW upwelling for 100% overcast sky at surface (w/m^2)');