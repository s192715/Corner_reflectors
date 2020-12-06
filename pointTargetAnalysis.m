%% data load in and arrangement
btest=load('ProfilesNew.mat');
slcParameters = xml2struct('C:/Users/kisum/Documents/30350 - remote sensing/corner/data/s1a-20150809t083228-vv-iw2-burst5-deramped.xml');
data_azi=[btest.Profile_CR4_Azi btest.Profile_CR5_Azi btest.Profile_CR6_Azi btest.Profile_CR14_Azi];
data_range=[btest.Profile_CR4_Range' btest.Profile_CR5_Range' btest.Profile_CR6_Range' btest.Profile_CR14_Range'];


%% power conversion introduced on the for loop
%P = abs((data_azi(:,1)))^2;
%[powMax , idx] = max(P,[],[1,2],'linear');
%P = P/ powMax;
%[R_i_Max] = ind2sub(size(P), idx);
%% plotting and calculation of pslr and pulse width 
cridx=[4 5 6 14];
rangeidx=[0 1 2 3];
aziidx=[1 2 3 4];
oversamp=128;
figure(100)
clf
for i=1:4

    P = abs((data_azi(:,i))).^2;
    [powMax , idx] = max(P,[],[1,2],'linear');
    P = P/ powMax;
    [Az_i_Max] = ind2sub(size(P), idx);
P_h = P;%ZI(:,Az_i_Max(:,:,i)-(i-1)*size(ZI,1),i);
    P = abs(data_range(:,i)).^2;
    [powMax , idx] = max(P,[],[1,2],'linear');
    P = P/ powMax;
    [R_i_Max] = ind2sub(size(P), idx);
P_v = P;%ZI(R_i_Max(:,:,i),:,i);
% Intesity profiles for both range and azimuth

inter_ticks_azi=[-8:8];
inter_ticks_range=[-8:8];

subplot(4,2,rangeidx(i)+i)
plot(pow2db(abs(P_h)),'color','blue','LineWidth',1);
title(['Range intensity profile for corner reflector ID:',num2str(cridx(i))]);
grid
xlim([0 size(P_h,1)]);
ylim([-50 0]);
ylabel('Power[dB]');
xl=xlabel('Interpolated range pixels');
xl.Position(2)=xl.Position(2)-0.1;
%xticks([0:16]*oversamp);
%xticklabels(inter_ticks_range);
%xtickangle(-20);

subplot(4,2,aziidx(i)+i)
plot(pow2db(abs(P_v)),'color','red','LineWidth',1);
title(['Azimuth intensity profile for corner reflector ID:',num2str(cridx(i))]);
grid
xlim([0 size(P_v,1)]);
ylim([-50 0]);
ylabel('Power[dB]');
xl=xlabel('Interpolated azimuth pixels');
xl.Position(2)=xl.Position(2)-0.1;
%xticks([0:8*2]*oversamp);
%xticklabels(inter_ticks_azi);
%xtickangle(-20);

%calc -3db bw
pow_azi=(pow2db(abs(P_v)));
pow_range=(pow2db(abs(P_h)));
dbbw=-3;
[val_max_azi idx_max_azi]=max(pow_azi);
[val_max_range idx_max_range]=max(pow_range);

[val_dbbw_azi_low(i),idx_dbbw_azi_low(i)]=min(abs(pow_azi(1:idx_max_azi)-dbbw));
[val_dbbw_range_low(i),idx_dbbw_range_low(i)]=min(abs(pow_range(1:idx_max_range)-dbbw));
[val_dbbw_azi_high(i),idx_dbbw_azi_high(i)]=min(abs(pow_azi(idx_max_azi:end)-dbbw));
[val_dbbw_range_high(i),idx_dbbw_range_high(i)]=min(abs(pow_range(idx_max_range:end)-dbbw));

dbbw_azi_low(i)=pow_azi(idx_dbbw_azi_low(i));
dbbw_range_low(i)=pow_range(idx_dbbw_range_low(i));

dbbw_azi_high(i)=pow_azi(idx_dbbw_azi_high(i)+idx_max_azi);
dbbw_range_high(i)=pow_range(idx_dbbw_range_high(i)+idx_max_range);

bw_3db_azi_inter_px(i)=idx_max_azi+ idx_dbbw_azi_high(i)-idx_dbbw_azi_low(i);
bw_3db_range_inter_px(i)=idx_max_range+ idx_dbbw_range_high(i)-idx_dbbw_range_low(i);

% pslr calculations
[peaks_azi,peaks_idx_azi]=findpeaks(pow_azi);
[peaks_azi_sorted,peaks_azi_sorted_idx]=sort(peaks_azi);
plsr_azi(i)=peaks_azi(peaks_azi_sorted_idx(end-1));

[peaks_range,peaks_idx_range]=findpeaks(pow_range);
[peaks_range_sorted,peaks_range_sorted_idx]=sort(peaks_range);
plsr_range(i)=peaks_range(peaks_range_sorted_idx(end-1));


%plotting pslr and -3dbbw
p1=[peaks_idx_range(peaks_range_sorted_idx(end-1)) peaks_range(peaks_range_sorted_idx(end-1))];
p2= [peaks_idx_range(peaks_range_sorted_idx(end-1)) 0];
subplot(4,2,rangeidx(i)+i)
line([peaks_idx_range(peaks_range_sorted_idx(end-1)) peaks_idx_range(peaks_range_sorted_idx(end-1))], [peaks_range(peaks_range_sorted_idx(end-1)) 0],'color','green','LineWidth',2,'linestyle','-.')
%text(peaks_idx_range(peaks_range_sorted_idx(end-1)),peaks_range(peaks_range_sorted_idx(end-1))/2,'\leftarrow PSLR')
%text(p1(1),p1(2), sprintf('(%.0f,%.0f)',p1))
%text(p2(1),p2(2), sprintf('(%.0f,%.0f)',p2))
p1=[idx_dbbw_range_low(i) -3];
p2=[idx_max_range+idx_dbbw_range_high(i) -3];
line([p1(1) p2(1)], [-3 -3],'color',[.5 .5 0],'LineWidth',2)
%text(p2(1),p2(2),'  \leftarrow 3 dB pulsewidth')
legend('intensity[dB]',['PSLR \approx ',num2str(round(abs(peaks_range(peaks_range_sorted_idx(end-1))))),' dB'],'-3 dB pulsewidth','location','NE')


p1=[peaks_idx_azi(peaks_azi_sorted_idx(end-1)) peaks_azi(peaks_azi_sorted_idx(end-1))];
p2= [ peaks_idx_azi(peaks_azi_sorted_idx(end-1)) 0];
subplot(4,2,aziidx(i)+i)
%text(peaks_idx_azi(peaks_azi_sorted_idx(end-1)),peaks_azi(peaks_azi_sorted_idx(end-1))/2,'\leftarrow PSLR')
line([peaks_idx_azi(peaks_azi_sorted_idx(end-1)) peaks_idx_azi(peaks_azi_sorted_idx(end-1))], [peaks_azi(peaks_azi_sorted_idx(end-1)) 0],'color','green','LineWidth',2,'linestyle','-.')
%text(p1(1),p1(2), sprintf('(%.0f,%.0f)',p1))
%text(p2(1),p2(2), sprintf('(%.0f,%.0f)',p2))
p1=[idx_dbbw_azi_low(i) -3];
p2=[ idx_max_azi+idx_dbbw_azi_high(i) -3];
line([p1(1) p2(1)], [-3 -3],'color',[.5 .5 0],'LineWidth',2)
%text(p2(1),p2(2),'  \leftarrow 3 dB pulsewidth')
legend('intensity[dB]',['PSLR \approx ',num2str(round(abs(peaks_azi(peaks_azi_sorted_idx(end-1))))),' dB'],'-3 dB pulsewidth','location','NE')

end
%%
% conversion from interpolated pixels to meters
bw_3db_azi_px=bw_3db_azi_inter_px/oversamp;
bw_3db_range_px=bw_3db_range_inter_px/oversamp;
prf=14.0660469;%from bjorn, just for the presentation
%prf=vg/str2num(slcParameters.burstParameterFile.slcParameters.pulseRepetitionFrequency.Text);
bw_3db_azi_m=bw_3db_azi_px*prf

rsf=physconst('lightspeed')/(2*str2num(slcParameters.burstParameterFile.slcParameters.rangeSamplingFrequency.Text));
bw_3db_range_m=bw_3db_range_px*rsf