%Reading in corner reflector geographic coordinates
%Format: IDnumber, Lat[deg], Lon[deg], Height[m].
CRgeo= readmatrix('C:\Users\Björn Bergsson\Desktop\DTUcourses\30350 RemoteSensing\ProjectCornerReflectors\Data\CRs_ITRF2008');

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


%reading radar parameters into structure:
slcParameters = xml2struct('C:\Users\Björn Bergsson\Desktop\DTUcourses\30350 RemoteSensing\ProjectCornerReflectors\Data\s1a-20150809t083228-vv-iw2-burst5-deramped.xml');

%Reading binary data matrix into complex matrix:
numOfColumns = str2num(slcParameters.burstParameterFile.rasterFileDescription.numOfColumns.Text);
slc = readBinFile('C:\Users\Björn Bergsson\Desktop\DTUcourses\30350 RemoteSensing\ProjectCornerReflectors\Data\s1a-20150809t083228-vv-iw2-burst5-deramped.slc', numOfColumns, 2);
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

%Finding the State Vector Time of each state vector (This is just to
%include the times on the x-axis in the plots)
SVT = num2cell(zeros(1,length(Range(1,:)))); %allocating, string of State Vector Times
for k = 1:length(Range(1,:))
    SVT{1,k} = string(extractBetween(slcParameters.burstParameterFile.orbitStateVectorList.orbitStateVector{k}.time.Text,14,20));
end


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

