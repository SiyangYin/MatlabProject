function [] = trubcr21AM(trubcr_21)
figure;
t=83:109;
trubcr_21_1=trubcr_21(:,:,:,1);
trubcr_21_2=trubcr_21(:,:,:,2);
trubcr_21_3=trubcr_21(:,:,:,3);
trubcr_21_4=trubcr_21(:,:,:,4);
trubcr_21_5=trubcr_21(:,:,:,5);
trubcr_21_6=trubcr_21(:,:,:,6);
trubcr_21_7=trubcr_21(:,:,:,7);
trubcr_21_8=trubcr_21(:,:,:,8);
trubcr_21_9=trubcr_21(:,:,:,9);
trubcr_21_10=trubcr_21(:,:,:,10);
trubcr_21_11=trubcr_21(:,:,:,11);
trubcr_21_12=trubcr_21(:,:,:,12);
trubcr_21_13=trubcr_21(:,:,:,13);
trubcr_21_14=trubcr_21(:,:,:,14);
trubcr_21_15=trubcr_21(:,:,:,15);
trubcr_21_16=trubcr_21(:,:,:,16);
trubcr_21_17=trubcr_21(:,:,:,17);
trubcr_21_18=trubcr_21(:,:,:,18);
trubcr_21_19=trubcr_21(:,:,:,19);
trubcr_21_20=trubcr_21(:,:,:,20);
trubcr_21_21=trubcr_21(:,:,:,21);
trubcr_21_22=trubcr_21(:,:,:,22);
trubcr_21_23=trubcr_21(:,:,:,23);
trubcr_21_24=trubcr_21(:,:,:,24);
trubcr_21_25=trubcr_21(:,:,:,25);
trubcr_21_26=trubcr_21(:,:,:,26);
trubcr_21_27=trubcr_21(:,:,:,27);
trubcr_21_(1)=mean(trubcr_21_1(trubcr_21_1~=0));
trubcr_21_(2)=mean(trubcr_21_2(trubcr_21_2~=0));
trubcr_21_(3)=mean(trubcr_21_3(trubcr_21_3~=0));
trubcr_21_(4)=mean(trubcr_21_4(trubcr_21_4~=0));
trubcr_21_(5)=mean(trubcr_21_5(trubcr_21_5~=0));
trubcr_21_(6)=mean(trubcr_21_6(trubcr_21_6~=0));
trubcr_21_(7)=mean(trubcr_21_7(trubcr_21_7~=0));
trubcr_21_(8)=mean(trubcr_21_8(trubcr_21_8~=0));
trubcr_21_(9)=mean(trubcr_21_9(trubcr_21_9~=0));
trubcr_21_(10)=mean(trubcr_21_10(trubcr_21_10~=0));
trubcr_21_(11)=mean(trubcr_21_11(trubcr_21_11~=0));
trubcr_21_(12)=mean(trubcr_21_12(trubcr_21_12~=0));
trubcr_21_(13)=mean(trubcr_21_13(trubcr_21_13~=0));
trubcr_21_(14)=mean(trubcr_21_14(trubcr_21_14~=0));
trubcr_21_(15)=mean(trubcr_21_15(trubcr_21_15~=0));
trubcr_21_(16)=mean(trubcr_21_16(trubcr_21_16~=0));
trubcr_21_(17)=mean(trubcr_21_17(trubcr_21_17~=0));
trubcr_21_(18)=mean(trubcr_21_18(trubcr_21_18~=0));
trubcr_21_(19)=mean(trubcr_21_19(trubcr_21_19~=0));
trubcr_21_(20)=mean(trubcr_21_20(trubcr_21_20~=0));
trubcr_21_(21)=mean(trubcr_21_21(trubcr_21_21~=0));
trubcr_21_(22)=mean(trubcr_21_22(trubcr_21_22~=0));
trubcr_21_(23)=mean(trubcr_21_23(trubcr_21_23~=0));
trubcr_21_(24)=mean(trubcr_21_24(trubcr_21_24~=0));
trubcr_21_(25)=mean(trubcr_21_25(trubcr_21_25~=0));
trubcr_21_(26)=mean(trubcr_21_26(trubcr_21_26~=0));
trubcr_21_(27)=mean(trubcr_21_27(trubcr_21_27~=0));
plot(t,trubcr_21_,'k-*','LineWidth',2);
set(gca,'XTick',83:109);
axis([83 109 320 360]);
legend('LW upwelling for clear sky at surface (w/m^2)','Location','NorthWest');
title({'Annual mean time series of LW upwelling for clear sky at surface (w/m^2)','in Montreal(45°30′N，73°40′W) from 1983 to 2009'});
xlabel('Years from 1983 to 2009 (minus 1900)');
ylabel('LW upwelling for clear sky at surface (w/m^2)');