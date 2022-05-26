%growthrate for each cell between two frames 
figure(4)
hold on

z=pulsedata.totvolG1_daughters;
growthrate_all=(diff(z,1,2))/3;
volumeG1=(z(:,1:(size(z,2)-1))+z(:,2:size(z,2)))/2;%pulsedata.totvolG1_daughters(:,[1:(endtimeanalysis-1)]);
idxValid=~isnan(growthrate_all);
growthrate_idx=growthrate_all(idxValid);
volumeG1_idx=volumeG1(idxValid);
[growthrate_all_mean,growthrate_all_error,volumeG1_binsmean]=KS_bindata_mean_20140506(volumeG1_idx,growthrate_idx,15);
plot(volumeG1,growthrate_all,'b','Linestyle','none','marker','o')
xlabel('volume [fl]')
ylabel('growthrate G1 [a.u.]')

hold off
figure(5)
hold on
errorbar(volumeG1_binsmean,growthrate_all_mean,growthrate_all_error,'r','LineWidth',3)
xlabel('volume [fl]')
ylabel('growthrate G1 [a.u.]')
% xlim([10 80])
% ylim([-0.5 1.6])
i=~isnan(growthrate_all_mean);
[p1,s1]=polyfit(volumeG1_binsmean(i),growthrate_all_mean(i),1);
fit1=polyval(p1,volumeG1_binsmean);
plot(volumeG1_binsmean,fit1,'b')
hold off

%% SG2M Data

%growthrate for mother+daughter 

figure(16)
hold on

z=pulsedata.totvolSG2M;
growthrate_tot=(diff(z,1,2))/3;
volumeSG2M=(z(:,1:(size(z,2)-1))+z(:,2:size(z,2)))/2;  % so that matrix has the same dimension
idxValid=~isnan(growthrate_tot);
growthrate_idx=growthrate_tot(idxValid);
volumetot_idx=volumeSG2M(idxValid);
[outliers]= find(volumetot_idx<100);
[growthrate_tot_mean,growthrate_tot_error,volumetot_binsmean]=KS_bindata_mean_20140506(volumetot_idx(outliers),growthrate_idx(outliers),25);
plot(volumetot_idx(outliers),growthrate_idx(outliers),'b','Linestyle','none','marker','o')
xlabel('volume [fl]')
ylabel('growthrate SG2M [a.u.]')
hold off

figure(17)
hold on
errorbar(volumetot_binsmean,growthrate_tot_mean,growthrate_tot_error,'b','LineWidth',3)
xlabel('volume [fL]')
ylabel('growth rate SG2M [a.u.]')
%xlim([0 150])
%ylim([-0.1 1])

[i]= find(volumetot_binsmean<70)    %fit one
[j]=~isnan(growthrate_tot_mean(i))
[p11,s11]=polyfit(volumetot_binsmean(j),growthrate_tot_mean(j),1);
fit11=polyval(p11,volumetot_binsmean(j));
plot(volumetot_binsmean(j),fit11,'r')

[ii]= find(volumetot_binsmean > 65) ;%fit two
vol=volumetot_binsmean(ii);
growth=growthrate_tot_mean(ii);
[jj]=~isnan(growth);
[p111,s111]=polyfit(vol(jj),growth(jj),0);
fit111=polyval(p111,vol(jj));
plot(vol(jj),fit111,'g')

i=~isnan(growthrate_tot_mean);  % combined fit
[p11,s11]=polyfit(volumetot_binsmean(i),growthrate_tot_mean(i),1);
fit11=polyval(p11,volumetot_binsmean(i));
plot(volumetot_binsmean(i),fit11,'r')
hold off