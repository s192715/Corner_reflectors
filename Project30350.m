%Reading in corner reflector geographic coordinates
%Format: IDnumber, Lat[deg], Lon[deg], Height[m].
CRgeo= readmatrix('C:/Users/kisum/Documents/30350 - remote sensing/corner/Data/CRs_ITRF2008.txt');

%% Calculating corner reflectors in ECEF coordinates

%WGS84 parameters:
f = 1/298.257223563;
a = 6378137;
b = a*(1-f);
e = sqrt((a^2 - b^2)/a^2);
e2 = sqrt((a^2 - b^2)/b^2);

N = zeros(1,length(CRgeo(:,1))); %allocating for radius of curvature
CR_ECEF = zeros(size(CRgeo)); %allocating for Corner Reflector ECEF coordinates
CR_ECEF(:,1) = CRgeo(:,1); %Adding the station numbers
for i = 1:length(CRgeo(:,1))%number of stations
    N(i) = a/sqrt(1-(e^2)*sind(CRgeo(i,2))^2); 
    CR_ECEF(i,2:4) = [(N(i)+CRgeo(i,4))*cosd(CRgeo(i,2))*cosd(CRgeo(i,3)) (N(i)+CRgeo(i,4))*cosd(CRgeo(i,2))*sind(CRgeo(i,3)) ((b^2/a^2)*N(i)+CRgeo(i,4))*sind(CRgeo(i,2))];
end

%We now have the corner reflectors in ECEF coordinates:
%CR_ECEF format: IDnumber, X[m], Y[m], Z[m].


%% reading radar parameters into structure:
slcParameters = xml2struct('C:/Users/kisum/Documents/30350 - remote sensing/corner/data/s1a-20150809t083228-vv-iw2-burst5-deramped.xml');

%Reading binary data matrix into complex matrix:
numOfColumns = str2num(slcParameters.burstParameterFile.rasterFileDescription.numOfColumns.Text);
slc = readBinFile('C:/Users/kisum/Documents/30350 - remote sensing/corner/data/s1a-20150809t083228-vv-iw2-burst5-deramped.slc', numOfColumns, 2);
%So the slc variable is our single look complex image

%% Computing the azimuth time and slant range from each Corner Reflector to all satellite state vectors:

%allocating for slant range and azimuth time, one line is for each CR and column denotes state vectors
Range = zeros(4,length(slcParameters.burstParameterFile.orbitStateVectorList.orbitStateVector));
DopplerEq = zeros(4,length(slcParameters.burstParameterFile.orbitStateVectorList.orbitStateVector));

for i = 1:length(CRgeo(:,1)) %number of Corner Reflectors (They are 4)
    for j = 1:length(Range(1,:)) %number of state vectors (They are 22)
        Range(i,j) = sqrt((CR_ECEF(i,2)-str2num(slcParameters.burstParameterFile.orbitStateVectorList.orbitStateVector{j}.x.Text))^2 + (CR_ECEF(i,3)-str2num(slcParameters.burstParameterFile.orbitStateVectorList.orbitStateVector{j}.y.Text))^2 + (CR_ECEF(i,4)-str2num(slcParameters.burstParameterFile.orbitStateVectorList.orbitStateVector{j}.z.Text))^2);
        DopplerEq(i,j) = dot([str2num(slcParameters.burstParameterFile.orbitStateVectorList.orbitStateVector{j}.vx.Text) str2num(slcParameters.burstParameterFile.orbitStateVectorList.orbitStateVector{j}.vy.Text) str2num(slcParameters.burstParameterFile.orbitStateVectorList.orbitStateVector{j}.vz.Text)], ...
            [[CR_ECEF(i,2) CR_ECEF(i,3) CR_ECEF(i,4)] - [str2num(slcParameters.burstParameterFile.orbitStateVectorList.orbitStateVector{j}.x.Text) str2num(slcParameters.burstParameterFile.orbitStateVectorList.orbitStateVector{j}.y.Text) str2num(slcParameters.burstParameterFile.orbitStateVectorList.orbitStateVector{j}.z.Text)]]);
    end
end

%% Finding the State Vector Time of each state vector (This is just to
%include the times on the x-axis in the plots)
SVT = num2cell(zeros(1,length(Range(1,:)))); %allocating, string of State Vector Times
for k = 1:length(Range(1,:))
    SVT{1,k} = string(extractBetween(slcParameters.burstParameterFile.orbitStateVectorList.orbitStateVector{k}.time.Text,14,20));
end
%%

%plotting slant range for each station
for i = 1:length(CRgeo(:,1)) %number of stations
    figure(i); plot(1:22,Range(i,:)/1000,'b* ') %diving by 1000 to have units in km
    xlabel('Azimuth Time')
    xticks([1:22])
    xticklabels({SVT{1},SVT{2},SVT{3},SVT{4},SVT{5},SVT{6},SVT{7},SVT{8},SVT{9},SVT{10},SVT{11},SVT{12},SVT{13},SVT{14},SVT{15} ...
    ,SVT{16},SVT{17},SVT{18},SVT{19},SVT{20},SVT{21},SVT{22}})
    xtickangle(90)
    ylabel('Slant range [km]')
    title(['Corner Reflector ID: ', num2str(CRgeo(i,1)),', 9th of August 2015.'])
end
%%
%plotting Doppler equation for each station
for i = 1:length(CRgeo(:,1)) %number of stations
    figure(i+4); plot(1:22,DopplerEq(i,:),'b* ')
    xlabel('Azimuth Time')
    xticks([1:22])
    xticklabels({SVT{1},SVT{2},SVT{3},SVT{4},SVT{5},SVT{6},SVT{7},SVT{8},SVT{9},SVT{10},SVT{11},SVT{12},SVT{13},SVT{14},SVT{15} ...
    ,SVT{16},SVT{17},SVT{18},SVT{19},SVT{20},SVT{21},SVT{22}})
    xtickangle(90)
    ylabel('Doppler equation $V_S \cdot (X_P - X_S)$ [m$^2$/s]','interpreter','latex')
    title(['Corner Reflector ID: ', num2str(CRgeo(i,1)),', 9th of August 2015.'])
end

%% testing interpolation in the dopler domain 
% this iis done by the function testPolyfit.m
n_max=5;%the highest order of polynoial for the testing
res_max=5;%the allowed threshold for the interpolation
i_true=11;%the index for wich is used for the testing
[residual, degree]=testPolyfit(n_max,i_true,res_max,DopplerEq);
%This result in a linear fit 
%%
%finding zero doppler with linear interpolation
steps=20e-7;%step size for the interpolation

for m = 1:size(DopplerEq,1)
zero=0;
[val(m),idx(m)]=min(abs(DopplerEq(m,:)-zero));
minVal(m)=DopplerEq(m,idx(m));
i=[idx(m)-1 , idx(m)+1];
fit  = polyfit(i, DopplerEq(m,i),1);%the fitting
i=[idx(m)-1:steps:idx(m)+1];
fit_val=polyval(fit,i);
[val_dopp(m),idx_dopp(m)]=min(abs(fit_val));

figure()
clf
hold on 
plot([i(1) i(end)],DopplerEq(m,[10 12]))
plot(i(idx_dopp(m)),val_dopp(m),'x')
xticks([i(1) i(end)])
xticklabels([SVT{10},SVT{12}])
xtickangle(90)
grid on
hold off
end

%%
%figure(9)
%clf
%hold on 
%plot([i(1) i(end)],DopplerEq(4,[10 12]))
%plot(i(idx_dopp),val_dopp,'+')
%xticks([i(1) i(end)])
%xticklabels([SVT{10},SVT{12}])
%xtickangle(90)
%hold off
%% calculate the time at found index
for m = 1:size(DopplerEq,1)
time(1)=datetime(SVT{1,idx(m)-1},'InputFormat','HH:mm:ss');
time(2)=datetime(SVT{1,idx(m)+1},'InputFormat','HH:mm:ss');
time_azi(m)=seconds(time(2)-time(1))/length(i)*idx_dopp(m);%whole interpolated time interval, divided by number of samples, multiplyed by the index for zero doppler
time_z_doppler(m)=time(1)+seconds(time_azi(m));
time_z_doppler_string=datestr(datenum(time_z_doppler(m)),'HH:mm:ss');%21 minutes are missing,dunno why
end
%%
figure(10)
clf
hold on 
plot([i(1) i(idx_dopp(1)) i(end)],[DopplerEq(4,10) val_dopp(1) DopplerEq(4,12)])
plot(i(idx_dopp(1)),val_dopp(1),'+')
%plot([i(1) i(idx_dopp) i(end)],[0 0 0])
xticks([i(1) i(idx_dopp(1)) i(end)])
xticklabels([SVT{10},time_z_doppler_string,SVT{12}])
xtickangle(90)
grid on
hold off

%%
%Read in Bjørns coordinates
corners=load('SampleLineDoppEq.mat');
%%
%impulse response analysis
data = slc;
P = abs(data).^2;
[powMax,idx]= max(P(:));
P_norm=P/powMax;

%[rowMax , colMax] = ind2sub(size(P), idx);
%fprintf('rowMax , colMax = %d, %d\n', rowMax , colMax);
%data patch extraction for interpolation
halfPatchWindowSize = 8;
overSampFactor = 128;
rowStart =round( corners.azi_pix_loc_1_corr - halfPatchWindowSize);
rowEnd = round(corners.azi_pix_loc_1_corr + halfPatchWindowSize);
colStart = round(corners.range_pix_loc_1_corr - halfPatchWindowSize);
colEnd = round(corners.range_pix_loc_1_corr + halfPatchWindowSize);
dataPatch = zeros(halfPatchWindowSize*2+1,halfPatchWindowSize*2+1,length(rowStart));
for i=1:length(rowStart)
dataPatch(:,:,i) = P_norm(rowStart(i) : rowEnd(i) , colStart(i) : colEnd(i));    
end

%% power conversions
ZI = interpft(interpft(dataPatch , size(dataPatch ,1)*overSampFactor ,1),size(dataPatch ,1)*overSampFactor ,2);
P = ZI.*conj(ZI);


[powMax , idx] = max(P,[],[1,2],'linear');


for i=1:size(P,3)
    P(:,:,i) = P(:,:,i)/ powMax(i);
end
[R_i_Max , Az_i_Max] = ind2sub(size(P), idx);





%% plot of interpolated patch before and after
figure()
imagesc([],[],abs(dataPatch(:,:,1).^0.3).^2)
xlabel("Azimuth (samples)");
ylabel("Range (samples)");
colormap('Gray');
colorbar;

figure()
imagesc([],[],abs(ZI(:,:,1).^0.3).^2)
xlabel("Azimuth (samples)");
ylabel("Range (samples)");
colormap('Gray');
colorbar;




%% 
cr_loc_azi=round(corners.azi_pix_loc_1_corr);
cr_loc_range=round(corners.range_pix_loc_1_corr);


%% plotting intesity profiles and calculating pslr and pulse width

cridx=[4 5 6 14];
rangeidx=[0 1 2 3];
aziidx=[1 2 3 4];
figure(100)
clf
for i=1:size(P,3)

P_h = ZI(:,Az_i_Max(:,:,i)-(i-1)*size(ZI,1),i);
P_v = ZI(R_i_Max(:,:,i),:,i);
% Intesity profiles for both range and azimuth

inter_ticks_azi=[-halfPatchWindowSize:halfPatchWindowSize]+cr_loc_azi(i);
inter_ticks_range=[-halfPatchWindowSize:halfPatchWindowSize]+cr_loc_range(i);

subplot(4,2,rangeidx(i)+i)
plot(pow2db(abs(P_h(:,1)/max(P_h(:,1))).^2),'color','blue','LineWidth',1);
title(['Range intensity profile for corner reflector #',num2str(cridx(i))]);
grid
%xlim([-30 994]);
ylim([-50 5]);
ylabel('Power[dB]');
xl=xlabel('Range pixels');
xl.Position(2)=xl.Position(2)-0.1;
xticks([0:halfPatchWindowSize*2]*overSampFactor);
xticklabels(inter_ticks_range);
xtickangle(-20);

subplot(4,2,aziidx(i)+i)
plot(pow2db(abs(P_v(1,:)/max(P_v(1,:))).^2),'color','red','LineWidth',1);
title(['Azimuth intensity profile for corner reflector #',num2str(cridx(i))]);
grid
%xlim([-18 1006]);
ylim([-50 5]);
ylabel('Power[dB]');
xl=xlabel('Azimuth pixels');
xl.Position(2)=xl.Position(2)-0.1;
xticks([0:halfPatchWindowSize*2]*overSampFactor);
xticklabels(inter_ticks_azi);
xtickangle(-20);

%calc -3db bw
pow_azi=(pow2db(abs(P_v/max(P_v)).^2));
pow_range=(pow2db(abs(P_h/max(P_h)).^2));
dbbw=-3;

[val_dbbw_azi_low(i),idx_dbbw_azi_low(i)]=min(abs(pow_azi(1:Az_i_Max(i)-(i-1)*size(ZI,1))-dbbw));
[val_dbbw_range_low(i),idx_dbbw_range_low(i)]=min(abs(pow_range(1:R_i_Max(i))-dbbw));
[val_dbbw_azi_high(i),idx_dbbw_azi_high(i)]=min(abs(pow_azi(Az_i_Max(i)-(i-1)*size(ZI,1):end)-dbbw));
[val_dbbw_range_high(i),idx_dbbw_range_high(i)]=min(abs(pow_range(R_i_Max(i):end)-dbbw));

dbbw_azi_low(i)=pow_azi(idx_dbbw_azi_low(i));
dbbw_range_low(i)=pow_range(idx_dbbw_range_low(i));

dbbw_azi_high(i)=pow_azi(idx_dbbw_azi_high(i)+Az_i_Max(i)-(i-1)*size(ZI,1)-1);
dbbw_range_high(i)=pow_range(idx_dbbw_range_high(i)+R_i_Max(i)-1);

bw_3db_azi_inter_px(i)=(idx_dbbw_azi_high(i)+Az_i_Max(i)-(i-1)*size(ZI,1))-idx_dbbw_azi_low(i);
bw_3db_range_inter_px(i)=(idx_dbbw_range_high(i)+R_i_Max(i))-idx_dbbw_range_low(i);

% pslr calculations
[peaks_azi,peaks_idx_azi]=findpeaks(pow_azi);
[peaks_azi_sorted,peaks_azi_sorted_idx]=sort(peaks_azi);
plsr_azi(i)=peaks_azi(peaks_azi_sorted_idx(end-1));

[peaks_range,peaks_idx_range]=findpeaks(pow_range);
[peaks_range_sorted,peaks_range_sorted_idx]=sort(peaks_range);
plsr_range(i)=peaks_range(peaks_range_sorted_idx(end-1));


%plotting pslr and -3dbbw
p1=[peaks_idx_range(peaks_range_sorted_idx(end-1)) peaks_range(peaks_range_sorted_idx(end-1))];
p2= [ peaks_idx_range(peaks_range_sorted_idx(end-1)) 0];
subplot(4,2,rangeidx(i)+i)
line([peaks_idx_range(peaks_range_sorted_idx(end-1)) peaks_idx_range(peaks_range_sorted_idx(end-1))], [peaks_range(peaks_range_sorted_idx(end-1)) 0],'color','green','LineWidth',2,'linestyle','-.')
%text(peaks_idx_range(peaks_range_sorted_idx(end-1)),peaks_range(peaks_range_sorted_idx(end-1))/2,'\leftarrow PSLR')
%text(p1(1),p1(2), sprintf('(%.0f,%.0f)',p1))
%text(p2(1),p2(2), sprintf('(%.0f,%.0f)',p2))
p1=[idx_dbbw_range_low(i) -3];
p2=[idx_dbbw_range_high(i)+R_i_Max(i) -3];
line([p1(1) p2(1)], [-3 -3],'color',[.5 .5 0],'LineWidth',2)
%text(p2(1),p2(2),'  \leftarrow 3 dB bandwidth')
legend('intensity[dB]',['PSLR = ',num2str(round(abs(peaks_range(peaks_range_sorted_idx(end-1))),2))],'-3 dB bandwidth','location','NE')


p1=[peaks_idx_azi(peaks_azi_sorted_idx(end-1)) peaks_azi(peaks_azi_sorted_idx(end-1))];
p2= [ peaks_idx_azi(peaks_azi_sorted_idx(end-1)) 0];
subplot(4,2,aziidx(i)+i)
%text(peaks_idx_azi(peaks_azi_sorted_idx(end-1)),peaks_azi(peaks_azi_sorted_idx(end-1))/2,'\leftarrow PSLR')
line([peaks_idx_azi(peaks_azi_sorted_idx(end-1)) peaks_idx_azi(peaks_azi_sorted_idx(end-1))], [peaks_azi(peaks_azi_sorted_idx(end-1)) 0],'color','green','LineWidth',2,'linestyle','-.')
%text(p1(1),p1(2), sprintf('(%.0f,%.0f)',p1))
%text(p2(1),p2(2), sprintf('(%.0f,%.0f)',p2))
p1=[idx_dbbw_azi_low(i) -3];
p2=[idx_dbbw_azi_high(i)+Az_i_Max(i)-(i-1)*size(ZI,1) -3];
line([p1(1) p2(1)], [-3 -3],'color',[.5 .5 0],'LineWidth',2)
%text(p2(1),p2(2),'  \leftarrow 3 dB bandwidth')
legend('intensity[dB]',['PSLR = ',num2str(round(abs(peaks_azi(peaks_azi_sorted_idx(end-1))),2))],'-3 dB bandwidth','location','NE')

end

%%
% conversion from interpolated pixels to meters
bw_3db_azi_px=bw_3db_azi_inter_px/overSampFactor;
bw_3db_range_px=bw_3db_range_inter_px/overSampFactor;
prf=14.0660469;%from bjorn, just for the presentation
%prf=vg/str2num(slcParameters.burstParameterFile.slcParameters.pulseRepetitionFrequency.Text);
bw_3db_azi_m=bw_3db_azi_px*prf

rsf=physconst('lightspeed')/(2*str2num(slcParameters.burstParameterFile.slcParameters.rangeSamplingFrequency.Text));
bw_3db_range_m=bw_3db_range_px*rsf

%sub pixel accuracy



