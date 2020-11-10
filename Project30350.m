%Corner reflectot porject on DTU 30350 Remote Sensing
%Team members: Björn Bergsson, Dániel Levente Bényei, Jesper Asbjørn Juul Kisum


%Reading in corner reflector geographic coordinates
%Format: IDnumber, Lat[deg], Lon[deg], Height[m].
%Bjorn
%matrixPath='C:\Users\Björn Bergsson\Desktop\DTUcourses\30350 RemoteSensing\ProjectCornerReflectors\Data\CRs_ITRF2008';
%Daniel
matrixPath='Q:\DOKSI\EGYETEM\DTU\3_semester\30350_Remote_sensing\project\data\CRs_ITRF2008';
%Jesper

CRgeo= readmatrix(matrixPath);

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
%Björn
%xmlPath='C:\Users\Björn Bergsson\Desktop\DTUcourses\30350 RemoteSensing\ProjectCornerReflectors\Data\s1a-20150809t083228-vv-iw2-burst5-deramped.xml'
%Daniel
xmlPath='Q:\DOKSI\EGYETEM\DTU\3_semester\30350_Remote_sensing\project\data\s1a-20150809t083228-vv-iw2-burst5-deramped.xml';
%Jesper

slcParameters = xml2struct(xmlPath);

%Reading binary data matrix into complex matrix:
numOfColumns = str2num(slcParameters.burstParameterFile.rasterFileDescription.numOfColumns.Text);

%Björn
%binFilePath='C:\Users\Björn Bergsson\Desktop\DTUcourses\30350 RemoteSensing\ProjectCornerReflectors\Data\s1a-20150809t083228-vv-iw2-burst5-deramped.slc';
%Daniel
binFilePath='Q:\DOKSI\EGYETEM\DTU\3_semester\30350_Remote_sensing\project\data\s1a-20150809t083228-vv-iw2-burst5-deramped.slc';
%Jesper
slc = readBinFile(binFilePath, numOfColumns, 2);
%So the slc variable is our single look complex image

%% Computing the Doppler equation and slant range from each Corner Reflector to all satellite state vectors:

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


%plotting slant range for each station, interpolates vectors 10,11 and 12 and finds the minimum value,
%returning the slant range closest approach and zero-doppler azimuth time for each Corner Reflector.
R_0 = zeros(1,length(CRgeo(:,1))); %allocating for slant range closest approach
Eta_0num = zeros(1,length(CRgeo(:,1))); %allocating for zero-doppler azimuth time, displays seconds from last state vector.
Eta_0str = num2cell(zeros(1,length(CRgeo(:,1)))); %displays absolute time
%in these vectors, first element is the value for CR nr. 1 (ID:4), 2nd element is the
%value for CR nr. 2 (ID:5) and so on.
for i = 1:length(CRgeo(:,1)) %number of stations
    figure(i); plot(1:length(Range(1,:)),Range(i,:)/1000,'b* ') %diving by 1000 to have units in km
    hold on;
    xlabel('Azimuth Time')
    xticks([1:length(Range(1,:))])
    xticklabels({SVT{1},SVT{2},SVT{3},SVT{4},SVT{5},SVT{6},SVT{7},SVT{8},SVT{9},SVT{10},SVT{11},SVT{12},SVT{13},SVT{14},SVT{15} ...
    ,SVT{16},SVT{17},SVT{18},SVT{19},SVT{20},SVT{21},SVT{22}})
    xtickangle(90)
    ylabel('Slant range [km]')
   
    pfit = polyfit(10:12,Range(i,10:12),2); %fitting vectors 10:12
    pval = polyval(pfit,10:10^(-7):12); %evaluating with resolution of a microsecond, using 10^(-7) because.
    plot(10:10^(-7):12,pval/1000,'k-')                                              %10 secs between vectors
    R_0(i) = min(pval)/1000; %units [km]
    Test_vec = find(pval==min(pval))/10^6; %if minimum is in more than one location,
    Eta_0num(i) = 33 + Test_vec(1);         %take the first minimum (could happen because the resolution is so fine).
    Eta_0str{1,i} = ['8:32:',  num2str(Eta_0num(i))];
    
    title({['Corner Reflector ID: ' num2str(CRgeo(i,1)) ', 9th of August 2015.'] [' $R_0 =$ ' num2str(R_0(i)) ' km and $\eta_0 =$ ' char(Eta_0str(1,i)) '.']},'interpreter','latex')
end
R_0_meters = R_0*1000; %to also have the slant ranges in [m]

%plotting Doppler equation for each station
for i = 1:length(CRgeo(:,1)) %number of stations
    figure(i+4); plot(1:length(Range(1,:)),DopplerEq(i,:),'b* ')
    xlabel('Azimuth Time')
    xticks([1:length(Range(1,:))])
    xticklabels({SVT{1},SVT{2},SVT{3},SVT{4},SVT{5},SVT{6},SVT{7},SVT{8},SVT{9},SVT{10},SVT{11},SVT{12},SVT{13},SVT{14},SVT{15} ...
    ,SVT{16},SVT{17},SVT{18},SVT{19},SVT{20},SVT{21},SVT{22}})
    xtickangle(90)
    ylabel('Doppler equation $V_S \cdot (X_P - X_S)$ [m$^2$/s]','interpreter','latex')
    title(['Corner Reflector ID: ', num2str(CRgeo(i,1)),', 9th of August 2015.'])
end

%-----------------------------------------------------------------------------------------------------
%Testing of fitting polynomial can removed to Test/slantRange.m ? What
%parameters do we need to pass? @Björn pls look at it later


%% Test section, fitting polynomial to slant ranges, experimenting with number of state vectors to fit.
%fitting vectors 10:12 gives best fit

figure(9); plot(1:length(Range(1,:)),Range(1,:)/1000,'b* ') %diving by 1000 to have units in km
hold on;
xlabel('Azimuth Time')
xticks([1:length(Range(1,:))])
xticklabels({SVT{1},SVT{2},SVT{3},SVT{4},SVT{5},SVT{6},SVT{7},SVT{8},SVT{9},SVT{10},SVT{11},SVT{12},SVT{13},SVT{14},SVT{15} ...
    ,SVT{16},SVT{17},SVT{18},SVT{19},SVT{20},SVT{21},SVT{22}})
xtickangle(90)
ylabel('Slant range [km]')
title(['Corner Reflector ID: ', num2str(CRgeo(1,1)),', 9th of August 2015.'])

[pfit S1] = polyfit(1:length(Range(1,:)),Range(1,:),2); %fitting all state vectors
[pval delta1] = polyval(pfit,1:10^(-6):length(Range(1,:)),S1); %we evaluate the polynomial at 10^(-6) spacing
plot(1:10^(-6):length(Range(1,:)),pval/1000,'r-')

[pfit2 S2] = polyfit(3:20,Range(1,3:20),2); %fitting vectors 3:20
[pval2 delta2] = polyval(pfit2,3:10^(-6):20,S2);
plot(3:10^(-6):20,pval2/1000,'g-')
[pfit3 S3] = polyfit(8:15,Range(1,8:15),2); %fitting vectors 8:15
[pval3 delta3] = polyval(pfit3,8:10^(-6):15,S3);
plot(8:10^(-6):15,pval3/1000,'b-')
[pfit4 S4] = polyfit(10:12,Range(1,10:12),2); %fitting vectors 10:12
[pval4 delta4] = polyval(pfit4,10:10^(-6):12,S4);
plot(10:10^(-6):12,pval4/1000,'k-')

legend('Range','Fit22','Fit3:20','Fit8:15','Fit10:12')

Norm_Of_Residuals = [S1.normr S2.normr S3.normr S4.normr];

%% Finding the radar coordinates of CR's (i.e. pixel numbers in azimuth and range)

%extracting radar parameters
f_s = str2num(slcParameters.burstParameterFile.slcParameters.rangeSamplingFrequency.Text); %range sampling frequency
PRF = str2num(slcParameters.burstParameterFile.slcParameters.pulseRepetitionFrequency.Text); %PulseRepetitionFreq
midIncAngle = str2num(slcParameters.burstParameterFile.slcParameters.sceneCentreIncidenceAngle.Text);
near_range = str2num(slcParameters.burstParameterFile.slcParameters.nearRange.Text);
ncols = str2num(slcParameters.burstParameterFile.rasterFileDescription.numOfColumns.Text);
nlines = str2num(slcParameters.burstParameterFile.rasterFileDescription.numOfRows.Text);
TimeLengthOfImage = nlines/PRF; % units: [s]

%calculating pixel spacings
dr = physconst('lightspeed')/(2*f_s);
%dx = calcVg()/

%calculating far and mid slangt ranges:
%to calculate far slant range, (R_far - R_near)/dr = nsample
far_range = near_range + dr*ncols;
mid_range = near_range + dr*(ncols/2);




%finding slant range pixel location for each CR:
range_pix_loc = (R_0_meters - near_range)/dr; %first element corresponds to first CR (ID: 4) and so on...

%finding azimuth pixel location for each CR:
AziStartTimeSeconds = str2num(string(extractBetween(slcParameters.burstParameterFile.slcParameters.azimuthStartTime.Text,19,28)));
Delta_AziTime = Eta_0num - AziStartTimeSeconds;
azi_pix_loc = Delta_AziTime*PRF;

%showing the amplitude image and the location of the CR's in radar coordinates
figure(10); imagesc(abs(slc),[0 1]);
colorbar
colormap(gray)
hold on;
plot(range_pix_loc,azi_pix_loc,'*y ')
xlabel('Range pixels')
ylabel('Azimuth pixels')

%labeling CR ID numbers to plot
text(range_pix_loc + 375,azi_pix_loc - 25,num2cell(CRgeo(:,1)'),'Fontsize', 15,'Color','yellow')

