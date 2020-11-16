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
plot([i(1) i(idx_dopp) i(end)],[DopplerEq(4,10) val_dopp DopplerEq(4,12)])
plot(i(idx_dopp),val_dopp,'+')
%plot([i(1) i(idx_dopp) i(end)],[0 0 0])
xticks([i(1) i(idx_dopp) i(end)])
xticklabels([SVT{10},time_z_doppler_string,SVT{12}])
xtickangle(90)
hold off



