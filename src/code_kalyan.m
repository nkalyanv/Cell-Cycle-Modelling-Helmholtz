%% G1Data

% 
% %figure 2e in paper
% %volume at budemergence vs volume at birth
% figure (1)
% hold on
% scatter(pulsedata.volumestartG1_daughters,pulsedata.volumeendG1_daughters, 'k','LineWidth',3)
% [volumebudemrg_mean,volumebudemrg_error,volumestartG1_binsmean]=KS_bindata_mean_20140916(pulsedata.volumestartG1_daughters,pulsedata.volumeendG1_daughters,10);
% errorbar(volumestartG1_binsmean,volumebudemrg_mean,volumebudemrg_error,'r','LineWidth',3);
% xlabel('volume at birth [fl]')
% ylabel('volume at bud emergence [fl]')
% xlim([10 80])
% ylim([10 80])
% [p,s]=polyfit(pulsedata.volumestartG1_daughters,pulsedata.volumeendG1_daughters,1)
% fit=polyval(p,pulsedata.volumestartG1_daughters)
% plot(pulsedata.volumestartG1_daughters,fit,'r')
% hold off
% 
% % %% growthrate for each cell between two frames 
% % figure(4)
% % hold on
% % 
% %z=pulsedata.volumeG1_daughters;
% z=pulsedata.totvolG1_daughters;
% growthrate_all=(diff(z,1,2))/3;
% volumeG1=pulsedata.volumeG1_daughters(:,[1:(endtimeanalysis-1)]);
% idxValid=~isnan(growthrate_all);
% growthrate_idx=growthrate_all(idxValid);
% volumeG1_idx=volumeG1(idxValid);
% [growthrate_all_mean,growthrate_all_error,volumeG1_binsmean]=KS_bindata_mean_20140506(volumeG1_idx,growthrate_idx,50);
% plot(volumeG1_idx,growthrate_idx,'b','Linestyle','none','marker','o')
% xlabel('volume [fl]')
% ylabel('growthrate G1 [a.u.]')
% hold off
% figure(2)
% hold on
% errorbar(volumeG1_binsmean,growthrate_all_mean,growthrate_all_error,'b','LineWidth',3)
% xlabel('volume [fl]')
% ylabel('growthrate G1 [a.u.]')
% % xlim([10 80])
% % ylim([-0.5 1.6])
% i=~isnan(growthrate_all_mean);
% [p1,s1]=polyfit(volumeG1_binsmean(i),growthrate_all_mean(i),0);
% fit1=polyval(p1,volumeG1_binsmean);
% plot(volumeG1_binsmean,fit1,'r')
% hold off
% 
% %% figure2c, rate of passing budemergence
% 
% 
A=pulsedata.volumeG1_daughters;
B=isnan(A);                        % all NaN set to 1 (pre budemergence =0, after budemergence=1)
for i=1:size(A,1)                  % for 1:rows
    k1(i)=find(~isnan(A(i,:)),1,'last');    %find last non-NaN element (tbudmergence)
  B(i,k1(i))=[1];                   %set last non-Nan element to 1
end
figure(3)
hold on
idxV=~isnan(A);
[B_mean,B_error,A_binsmean]=KS_bindata_mean_20140916(A(idxV),B(idxV),30);
plot(A,B,'b','Linestyle','none','marker','o')
%errorbar(A_binsmean,B_mean,B_error,'r','LineWidth',3)
%xlim([0 90])
%ylim([0 0.3]);
xlabel('volume [fl]')
ylabel('probability of passing budemergence')
idxV=~isnan(B_mean);
[p3,s3]=polyfit (A_binsmean(idxV),B_mean(idxV),1);
fit3=polyval(p3,A_binsmean(idxV));
%plot(A_binsmean(idxV),fit3,'b')
hold off
% 
% 
% 
% B_mean=1-e^(-kpre*t)
figure(4)
kpre=[];
kpre=-(log(1-B_mean(idxV)))/3;
hold on
plot (A_binsmean(idxV),kpre,'b')
[p4,s4]=polyfit (A_binsmean(idxV),kpre,1);
fit4=polyval(p4,A_binsmean(idxV));
plot(A_binsmean(idxV),fit4,'r')
xlabel('volume [fl]')
ylabel('kpre')
 xlim([0 70])
ylim([0 0.06]);

hold off
% 


% probability of division for whole cell (mother daughter)


A=pulsedata.totvolSG2M;
B=isnan(A);                        % all NaN set to 1 (pre division =0, after division=1)
for i=1:size(A,1)                  % for 1:rows
    k1(i)=find(~isnan(A(i,:)),1,'last');    %find last non-NaN element (tcytokinesis)
  B(i,k1(i))=[1];                   %set last non-Nan element to 1
end
figure(18)
hold on
idxV=~isnan(A);
[B_mean,B_error,pulsedata.volbud_binsmean]=KS_bindata_mean_20140916(pulsedata.volbud(idxV),B(idxV),40);
%plot(A,B,'b','Linestyle','none','marker','o')
errorbar(pulsedata.volbud_binsmean,B_mean,B_error,'r','LineWidth',3)
% xlim([0 120])
%ylim([0 0.3]);
xlabel('bud volume [fl]')
ylabel('probability of passing division')
idxV=~isnan(B_mean);
[p31,s31]=polyfit (pulsedata.volbud_binsmean(idxV),B_mean(idxV),1);
fit31=polyval(p31,pulsedata.volbud_binsmean(idxV));
plot(pulsedata.volbud_binsmean(idxV),fit31,'b')
hold off

%B_mean=1-e^(-kpost*t)

figure(19)
kpost=[];
kpost=-(log(1-B_mean(idxV)))/3;
hold on
plot (pulsedata.volbud_binsmean(idxV),kpost,'b')
x=pulsedata.volbud_binsmean(idxV);
[i]= find(x>10)    
[p41,s41]=polyfit (x(i),kpost(i),1);
fit41=polyval(p41,x(i));
plot(x(i),fit41,'r')
xlabel(' bud volume [fl]')
ylabel('kpost ')
%xlim([10 70])
%ylim([0 0.06]);
hold off
