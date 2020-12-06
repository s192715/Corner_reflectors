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

%% TEST SECTION ONLY, fitting slant ranges (method 2) by leaving out StateVector11
for i = 1:length(CRgeo(:,1)) %number of stations
    figure(i); plot(1:length(Range(1,:)),Range(i,:)/1000,'b* ') %diving by 1000 to have units in km
    hold on;
    xlabel('Azimuth Time')
    xticks([1:length(Range(1,:))])
    xticklabels({SVT{1},SVT{2},SVT{3},SVT{4},SVT{5},SVT{6},SVT{7},SVT{8},SVT{9},SVT{10},SVT{11},SVT{12},SVT{13},SVT{14},SVT{15} ...
    ,SVT{16},SVT{17},SVT{18},SVT{19},SVT{20},SVT{21},SVT{22}})
    xtickangle(90)
    ylabel('Slant range [km]')
   
    [n_x res error] = FitTest2_30350(Range(i,:),20,11,10^(-3)); 
    
    i_true = 11;
if rem(n_x,2) == 1 %if n is odd
        i_x = [i_true-round((n_x+1)/2):1:i_true+round((n_x+1)/2)];
        i_x= i_x(i_x~=i_true);
else %if n is even
        i_x = [i_true-(n_x/2):1:i_true+((n_x+2)/2)];
        i_x = i_x(i_x~=i_true);
end

    
    X_axis = 1:length(Range(1,:));
    [pfit_x S_x] = polyfit(X_axis(i_x),Range(i,i_x)/1000,n_x); 
    [pval_x] = polyval(pfit_x,1:10^(-4):length(Range(i,:)),S_x); 
    
    %pfit = polyfit(10:12,Range(i,10:12),2); %fitting vectors 10:12
    %pval = polyval(pfit,10:10^(-7):12); %evaluating with resolution of a microsecond, using 10^(-7) because.
    plot(1:10^(-4):length(Range(i,:)),pval_x,'k-')                                              %10 secs between vectors
    R_0(i) = min(pval_x); %units [km]
    Test_vec = find(pval_x==min(pval_x))/10^3; %if minimum is in more than one location,
    Eta_0num(i) = Test_vec(1);         %take the first minimum (could happen because the resolution is so fine).
    Eta_0num = rem(Eta_0num,57); %azi time in seconds after 8:32
    Eta_0str{1,i} = ['8:32:',  num2str(Eta_0num(i))];
    %first state vector starts at 8:31, but here we start at 8:32 and take
    %the remainder of Eta_0num after dividing by 60
    
    title({['Corner Reflector ID: ' num2str(CRgeo(i,1)) ', 9th of August 2015.'] [' $R_0 =$ ' num2str(R_0(i)) ' km and $\eta_0 =$ ' char(Eta_0str(1,i)) '.']},'interpreter','latex')

end
R_0_meters = R_0*1000; %to also have the slant ranges in [m]
%Eta_0num is now in seconds from first state vector

%difference in range between this method and the above is about 0.3 m (0.6
%m average)
%difference in range pix is about 0.1, or 0.25 in average.
%but we're getting unstable polynomials when using this method so we're
%going with the other.

%difference in azi-time (Eta_0num) is about millisecond (we put a millisecond
%threshold), average -0.001621. 
%Difference in azi_pix_loc about 2/3 of a pix, mean -0.79.

%Can experiment more by putting lower threshold
%was the point maybe all along to just fit parabola to the curve?
%[we at least get unstable polynomials when fitting higher orders]

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

legend('Range','Fit1-22','Fit3-20','Fit8-15','Fit10-12')

Norm_Of_Residuals = [S1.normr S2.normr S3.normr S4.normr];
title({['Corner Reflector ID: ', num2str(CRgeo(1,1)),', 9th of August 2015.'],['|r|_{1-22} = ' , num2str(round(S1.normr/1000,1)), ' km, |r|_{3-20} = ' , num2str(round(S2.normr/1000,1)), ' km.'],['|r|_{8-15} = ' , num2str(round(S3.normr/1000,1)), ' km, |r|_{10-12} = ' , num2str(round(S4.normr/1000,1)), ' km.']})

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
range_pix_loc = 1 + (R_0_meters - near_range)/dr; %first element corresponds to first CR (ID: 4) and so on...

%finding azimuth pixel location for each CR:
AziStartTimeSeconds = str2num(string(extractBetween(slcParameters.burstParameterFile.slcParameters.azimuthStartTime.Text,19,28)));
Delta_AziTime = Eta_0num - AziStartTimeSeconds;
azi_pix_loc = 1 + Delta_AziTime*PRF;

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


%% Working on method 1, interpolating satellite position coordinates

X_sat = zeros(1,length(Range(1,:))); %allocating, X positions for all satellite state vectors.
Y_sat = zeros(1,length(Range(1,:)));
Z_sat = zeros(1,length(Range(1,:)));

for i = 1:length(Range(1,:)) %getting position coordinates for all satellite state vectors
    X_sat(i) = str2num(slcParameters.burstParameterFile.orbitStateVectorList.orbitStateVector{i}.x.Text);
    Y_sat(i) = str2num(slcParameters.burstParameterFile.orbitStateVectorList.orbitStateVector{i}.y.Text);
    Z_sat(i) = str2num(slcParameters.burstParameterFile.orbitStateVectorList.orbitStateVector{i}.z.Text);
end

%Using FitTest30350 function, finding what degree polynomial best fits the data
%here we fitted all the state vectors
% [n_x NormRes_x] = FitTest30350(1:length(X_sat),X_sat,21,10^(-3));
% [n_y NormRes_y] = FitTest30350(1:length(Y_sat),Y_sat,21,10^(-3));
% [n_z NormRes_z] = FitTest30350(1:length(Z_sat),Z_sat,21,10^(-3));

%% Improved method of fitting the data, leaves out vector 11 and predicts it value, to get a direct error estimation.

%Using FitTest2_30350 function, finding what degree polynomial best fits the data
[n_x NormRes_x error_x] = FitTest2_30350(X_sat/1000,20,11,10^(-8)); %divided by 1000 to have units in km (does not matter if in meters or km)
[n_y NormRes_y error_y] = FitTest2_30350(Y_sat/1000,20,11,10^(-8));
[n_z NormRes_z error_z] = FitTest2_30350(Z_sat/1000,20,11,10^(-8));
%we implement threshold, otherwise polynomial degree becomes too high and
%might result in unstable fit.

%finding indexes i_x, i_y, i_z
i_true = 11;
if rem(n_x,2) == 1 %if n is odd
        i_x = [i_true-round((n_x+1)/2):1:i_true+round((n_x+1)/2)];
        i_x= i_x(i_x~=i_true);
else %if n is even
        i_x = [i_true-(n_x/2):1:i_true+((n_x+2)/2)];
        i_x = i_x(i_x~=i_true);
end

if rem(n_y,2) == 1 %if n is odd
        i_y = [i_true-round((n_y+1)/2):1:i_true+round((n_y+1)/2)];
        i_y= i_y(i_y~=i_true);
else %if n is even
        i_y = [i_true-(n_y/2):1:i_true+((n_y+2)/2)];
        i_y = i_y(i_y~=i_true);
end

if rem(n_z,2) == 1 %if n is odd
        i_z = [i_true-round((n_z+1)/2):1:i_true+round((n_z+1)/2)];
        i_z= i_z(i_z~=i_true);
else %if n is even
        i_z = [i_true-(n_z/2):1:i_true+((n_z+2)/2)];
        i_z = i_z(i_z~=i_true);
end

%polynomial evaluated at a millisecond interval
X_sat_xaxis = 1:length(X_sat);
[pfit_x S_x] = polyfit(X_sat_xaxis(i_x),X_sat(i_x)/1000,n_x); 
[pval_x] = polyval(pfit_x,1:10^(-4):length(X_sat),S_x); 

Y_sat_xaxis = 1:length(Y_sat);
[pfit_y S_y] = polyfit(Y_sat_xaxis(i_y),Y_sat(i_y)/1000,n_y); 
[pval_y] = polyval(pfit_y,1:10^(-4):length(Y_sat),S_y); 

Z_sat_xaxis = 1:length(Z_sat);
[pfit_z S_z] = polyfit(Z_sat_xaxis(i_z),Z_sat(i_z)/1000,n_z); 
[pval_z] = polyval(pfit_z,1:10^(-4):length(Z_sat),S_z); 

%%
%plotting the position coordinates and the fitted polynomial (polynomial
%might be generated by fitting n+1 state vectors with n=5 degree
%polynomial, but it is evaluated over the full x-axis (all azimuth time).
%And it predicts value at vector 11 really well.
figure(11); plot(1:length(X_sat),X_sat/1000,'b*')
hold on;
plot(1:10^(-4):length(X_sat),pval_x,'r-')
xlabel('Azimuth Time')
xticks([1:length(Range(1,:))])
xticklabels({SVT{1},SVT{2},SVT{3},SVT{4},SVT{5},SVT{6},SVT{7},SVT{8},SVT{9},SVT{10},SVT{11},SVT{12},SVT{13},SVT{14},SVT{15} ...
    ,SVT{16},SVT{17},SVT{18},SVT{19},SVT{20},SVT{21},SVT{22}})
xtickangle(90)
ylabel('Satellite X position [km]')

figure(12); plot(1:length(Y_sat),Y_sat/1000,'b*')
hold on;
plot(1:10^(-4):length(Y_sat),pval_y,'r-')
xlabel('Azimuth Time')
xticks([1:length(Range(1,:))])
xticklabels({SVT{1},SVT{2},SVT{3},SVT{4},SVT{5},SVT{6},SVT{7},SVT{8},SVT{9},SVT{10},SVT{11},SVT{12},SVT{13},SVT{14},SVT{15} ...
    ,SVT{16},SVT{17},SVT{18},SVT{19},SVT{20},SVT{21},SVT{22}})
xtickangle(90)
ylabel('Satellite Y position [km]')

figure(13); plot(1:length(Z_sat),Z_sat/1000,'b*')
hold on;
plot(1:10^(-4):length(Y_sat),pval_z,'r-')
xlabel('Azimuth Time')
xticks([1:length(Range(1,:))])
xticklabels({SVT{1},SVT{2},SVT{3},SVT{4},SVT{5},SVT{6},SVT{7},SVT{8},SVT{9},SVT{10},SVT{11},SVT{12},SVT{13},SVT{14},SVT{15} ...
    ,SVT{16},SVT{17},SVT{18},SVT{19},SVT{20},SVT{21},SVT{22}})
xtickangle(90)
ylabel('Satellite Z position [km]')



%then when we have the zero doppler time, we find the value of the
%polynomials at that time (use pval_x pval_y and pval_z).

Time_BB = Eta_0num-33; %zero-doppler times obtained from method 2

AziTimesMethod1 = load('matlab.mat');
Time_Jesp = AziTimesMethod1.time_z_doppler.Second -33; %times obtained from method 1 in seconds from last minute
%we do -33 because we want seconds since state vector 10

%finding the value at x-axis corresponding to the azimuth zero-doppler time
time_axis_nr = round((90 + Time_Jesp)*1000); %number of seconds since first vector times 1000 to get milliseconds
%rounded to nearest millisecond
%first element corresponds to CR ID:4 and so on...

%Satellite ECEF coordinates at zero-doppler time, first element corresponds to CR ID:4 and so on...
X_eta0 = pval_x(time_axis_nr); %values in km
Y_eta0 = pval_y(time_axis_nr);
Z_eta0 = pval_z(time_axis_nr);

%Different representation:
%ECEF coordinates X Y Z of the satellite at zero-doopler time for each CR
Sat_CR4 = [X_eta0(1) Y_eta0(1) Z_eta0(1)]*1000; %multiply by 1000 to get values in [m]
Sat_CR5 = [X_eta0(2) Y_eta0(2) Z_eta0(2)]*1000;
Sat_CR6 = [X_eta0(3) Y_eta0(3) Z_eta0(3)]*1000;
Sat_CR14 = [X_eta0(4) Y_eta0(4) Z_eta0(4)]*1000; 

%calculating the range of closest approach:
R0_CR4 = norm(CR_ECEF(1,2:4)-Sat_CR4);
R0_CR5 = norm(CR_ECEF(2,2:4)-Sat_CR5);
R0_CR6 = norm(CR_ECEF(3,2:4)-Sat_CR6);
R0_CR14 = norm(CR_ECEF(4,2:4)-Sat_CR14);

%string of R0 distances, first element corresponds to R0 for CR ID4 and so on...
R0_method1 = [R0_CR4 R0_CR5 R0_CR6 R0_CR14];

%comparison between method1 and method2
DeltaR0 = R0_method1 - R_0_meters; %~m difference, all R_0_meters value (obtained by method 2), seem to be a bit larger than from method 1.
DeltaEta0 = AziTimesMethod1.time_z_doppler.Second - Eta_0num; %~2 ms difference

%calculating pixel locations for method 1
range_pix_loc_1 = 1 + (R0_method1 - near_range)/dr; %first element corresponds to first CR (ID: 4) and so on...
Delta_AziTime_1 =  (Time_Jesp + 33) - AziStartTimeSeconds;
azi_pix_loc_1 = 1 + Delta_AziTime_1*PRF;


%showing the amplitude image and the location of the CR's in radar coordinates
figure(14); imagesc(abs(slc),[0 1]);
colorbar
colormap(gray)
hold on;
plot(range_pix_loc,azi_pix_loc,'*y ') %method number 2
plot(range_pix_loc_1,azi_pix_loc_1,'*r ') %method number 1, number 1 seems to be more accurate.
xlabel('Range pixels')
ylabel('Azimuth pixels')

%labeling CR ID numbers to plot
text(range_pix_loc + 375,azi_pix_loc - 25,num2cell(CRgeo(:,1)'),'Fontsize', 15,'Color','yellow')
%text(range_pix_loc_1 + 375,azi_pix_loc_1 - 25,num2cell(CRgeo(:,1)'),'Fontsize', 15,'Color','red')

%Comparing with John's function llh2sali
slcParameters2 = ParseMatSlcPar('C:\Users\Björn Bergsson\Desktop\DTUcourses\30350 RemoteSensing\ProjectCornerReflectors\Data\s1a-20150809t083228-vv-iw2-burst5-deramped.xml',1);
sali = llh2sali(CRgeo(:,2:4),slcParameters2);
plot(sali(:,1)',sali(:,2)','*b ') 
legend('RangeEq', 'DopplerEq', 'llh2sali')
% plot(range_pix_loc_FromProfiles,azi_pix_loc_FromProfiles,'*m ') 
% legend('RangeEq', 'DopplerEq', 'llh2sali','InterpolatedProfiles')


%% Calculating azimuth bistatic delay for each CR
% c = physconst('lightspeed'); %[m/s]
% Delta_t = (2*R0_method1)/c;
% 
% %to find Azimuth Offset due to bistatic delay we need the satellite velocity
% %run lines 170-298 with lines 178-180 changed to vx, vy and vz.
% %Then vectors Sat_CR4,...,Sat_CR14 contain the satellite velocity vx, vy, vz at zero doppler time.
% 
% %magnitude of the satellite velocity at zero doppler time, first element for CR4 and so on...
% Vmag = [norm(Sat_CR4) norm(Sat_CR5) norm(Sat_CR6) norm(Sat_CR14)]; units m/s
Vmag = [7587.79913914198, 7587.74023321573, 7587.78923080371, 7587.79300155905]; 

% AziOffset = Vmag.*Delta_t/2;
% %Results, ~22 m for all CR:
%AziOffset = [22.1795947065058, 21.9925382558491, 21.8735187482024, 21.6575479310834];

%ESTIMATING THE AziOffset IN A MORE THOROUGH WAY (Based on 2020 paper)
rank = 8;
Delta_PRI = 1/slcParameters2.PRF; %Pulse Repetition Interval, obtained from Pulse Repetition Frequency.
tau_0 = rank*Delta_PRI;

Delta_swst = 1.620099906725725e-04;
tau=@(k) tau_0 + Delta_swst + k/slcParameters2.fs;

%<startTime>2015-08-09T08:32:28.049414</startTime>, time in seconds from midnight: 30748.049414
%r_PRI0 = 28.049414;

%t_stop_go = @(i) r_PRI0 + i*Delta_PRI;

%t_zero_Doppler_1 = t_stop_go(azi_pix_loc_1) + tau(range_pix_loc_1)./2 - tau_0;

%%%% THIS IS THE CORRECT CORRECTION, eq. 21.
%the first term in eq. 21 had already been applied to the data, so we just
%apply the latter 3 terms to our geolocation to match it with the data.

t_zero_Doppler_1 = tau(slcParameters2.nsa/2)/2 + tau(range_pix_loc_1)./2 - tau_0;
%t_zero_Doppler_1 = tau(range_pix_loc_1)./2 - tau_0;
AziOffset_1 = Vmag.*t_zero_Doppler_1; %for method 1, DopplerEq

t_zero_Doppler_2 = tau(slcParameters2.nsa/2)/2 + tau(range_pix_loc)./2 - tau_0;
%t_zero_Doppler_2 = tau(range_pix_loc)./2 - tau_0;
AziOffset_2 = Vmag.*t_zero_Doppler_2; %for method 2, fitting range curve

%Important thing to realize here is that we were calculating the real zero_doppler time (eq. 21),
%but the time the data is registered at is the received times (later time than real zero-doppler),
%so to correct the data one has to find the difference between received and zero_doppler_time and
%subtract the (received time - zero_doppler_time)*Vmag from the data to correct for the Azimuth Offset.
%The data we have has already had a bulk correction i.e. been corrected for t_IPF (first term in
%eq.21), so to match the geolocation with the data we need to add the last
%three terms to the geolocation.


%satellite ground velocity
Vg = CalcVg(slcParameters2, slcParameters2.nsa/2, slcParameters2.nli/2);
%Azimuth pixel resolution
dx = Vg/slcParameters2.PRF;

%estimating range tropospheric delay. From VMF3:
%https://vmf.geo.tuwien.ac.at/
%https://vmf.geo.tuwien.ac.at/trop_products/GRID/1x1/VMF3/VMF3_EI/2015/VMF3_20150809.H06

%lat lon ah aw zhd zwd
%-27.5  151.5  0.00125715  0.00047432  2.2131  0.0580
zhd = 2.2131; %zenith delay hydrostatic
zwd = 0.0580; %zenith delay wet

incAngle = str2double(slcParameters.burstParameterFile.slcParameters.sceneCentreIncidenceAngle.Text);
Range_delay = (zhd+zwd)/cosd(incAngle); %units m

%corrected pixel locations
range_pix_loc_1_corr = range_pix_loc_1 + Range_delay/dr; %corrected for method  nr 1 (DopplerEq)
azi_pix_loc_1_corr = azi_pix_loc_1 - AziOffset_1/dx;
range_pix_loc_2_corr = range_pix_loc + Range_delay/dr; %corrected for method nr 2 (fitting range curve)
azi_pix_loc_2_corr = azi_pix_loc - AziOffset_2/dx;
%%Bulk correction is applied with respect to the mid range, and our CR's
%%are all behind mid-range, so bulk correction actually sets the received
%%time to before zero-doppler time. So we subtract the small correction from the geolocation to match with data. 





%plot with corrected positions
figure(15); imagesc(abs(slc),[0 1]);
colorbar
colormap(gray)
hold on;
plot(range_pix_loc_2_corr,azi_pix_loc_2_corr,'*y ') %method number 2, Fitting Range curve
plot(range_pix_loc_1_corr,azi_pix_loc_1_corr,'*r ') %method number 1, DopplerEq
plot(sali(:,1)',sali(:,2)','*b ')  %John's function
xlabel('Range pixels')
ylabel('Azimuth pixels')
%labeling CR ID numbers to plot
text(range_pix_loc + 375,azi_pix_loc - 25,num2cell(CRgeo(:,1)'),'Fontsize', 15,'Color','yellow')



%% TEST SECTION ONLY, used to create function Response30350, constructing how to interpolate the azimuth and range profiles.

%starting with CR4 azimuth profile, _azi in variable names denotes that
%we're working on the azimuth profile
oversample_azi = 32;
profile_pixels_azi = 1295:1311; %choose these numbers from the amplitude image, take a reasonable interval.
profile_length_azi = length(profile_pixels_azi);

ProfileFirst_azi = 11260; %interpolating azimuth profiles at every range pixel 11260:11280,
ProfileLast_azi = 11280;    %finding at which range pixel the maximum is.

MaximumAbs_azi = zeros(1,length(ProfileFirst_azi:ProfileLast_azi)); %allocating for maximum values
k = 1;
for i = [ProfileFirst_azi:1:ProfileLast_azi] %this loop finds at which range pixel the azimuth profile has the greatest peak magnitude
    Test_profile = interpft(slc(profile_pixels_azi,i),oversample_azi*profile_length_azi); %oversampling factor 16
    MaximumAbs_azi(k) = max(abs(Test_profile));
    k = k+1;
end
At_pixel_azi = ProfileFirst_azi + find(MaximumAbs_azi == max(MaximumAbs_azi)) -1; %at what range pixel the profile is taken, i.e. where maximum value of all azimuth profiles is found

Profile_CR4_Azi = interpft(slc(profile_pixels_azi,At_pixel_azi),oversample_azi*profile_length_azi);

figure(15); plot(1:oversample_azi*profile_length_azi,abs(Profile_CR4_Azi),'b* ')
xlabel('Azimuth pixel')
xticks([1:oversample_azi:oversample_azi*profile_length_azi])
xticklabels({profile_pixels_azi})
xtickangle(90)
ylabel('Magnitude')
title({['Corner reflector ID: 4. Azimuth profile at range pixel nr. ' num2str(At_pixel_azi) '.'] ['Peak magnitude: ' num2str(max(abs(Profile_CR4_Azi))) '.']})


%Now doing CR4 range profile, _r in variable names denotes that
%we're working on the range profile
oversample_r = 32;
profile_pixels_r = 11260:11280;
profile_length_r = length(profile_pixels_r);

ProfileFirst_r = 1295; %selecting ranges of azimuth pixels to try range profiles.
ProfileLast_r = 1311; 

MaximumAbs_r = zeros(1,length(ProfileFirst_r:ProfileLast_r)); %allocating for maximum values
k = 1;
for i = [ProfileFirst_r:1:ProfileLast_r] %this loop finds at which range pixel the azimuth profile has the greatest peak magnitude
    Test_profile = interpft(slc(i,profile_pixels_r),oversample_r*profile_length_r); %oversampling factor 16
    MaximumAbs_r(k) = max(abs(Test_profile));
    k = k+1;
end
At_pixel_r = ProfileFirst_r + find(MaximumAbs_r == max(MaximumAbs_r)) -1; %at what range pixel the profile is taken, i.e. where maximum value of all azimuth profiles is found

Profile_CR4_Range = interpft(slc(At_pixel_r,profile_pixels_r),oversample_r*profile_length_r);

%for dB: 10*log10(abs(Profile_CR4_Range)/max(abs(Profile_CR4_Range))
figure(16); plot(1:oversample_r*profile_length_r,abs(Profile_CR4_Range),'b* ')
xlabel('Range pixel')
xticks([1:oversample_r:oversample_r*profile_length_r])
xticklabels({profile_pixels_r})
xtickangle(90)
ylabel('Magnitude')
title({['Corner reflector ID: 4. Range profile at azimuth pixel nr. ' num2str(At_pixel_r) '.'] ['Peak magnitude: ' num2str(max(abs(Profile_CR4_Range))) '.']})

figure(17); plot3(ones(1,oversample_azi*profile_length_azi)*At_pixel_azi, profile_pixels_azi(1):(1/oversample_azi):profile_pixels_azi(end) +1-(1/oversample_azi), abs(Profile_CR4_Azi),'b- ')
hold on;
xlabel('Range pixel')
ylabel('Azimuth Pixel')
zlabel('Magnitude')
grid on;
plot3(profile_pixels_r(1):(1/oversample_r):profile_pixels_r(end) +1-(1/oversample_r),ones(1,oversample_r*profile_length_r)*At_pixel_r,abs(Profile_CR4_Range),'r- ')

%in dB:
figure(100); plot3(ones(1,oversample_azi*profile_length_azi)*At_pixel_azi, profile_pixels_azi(1):(1/oversample_azi):profile_pixels_azi(end) +1-(1/oversample_azi), 10*log10(abs(Profile_CR4_Azi)/max(abs(Profile_CR4_Azi))),'b- ')
hold on;
xlabel('Range pixel')
ylabel('Azimuth Pixel')
zlabel('dB')
zlim([-20 0])
grid on;
plot3(profile_pixels_r(1):(1/oversample_r):profile_pixels_r(end) +1-(1/oversample_r),ones(1,oversample_r*profile_length_r)*At_pixel_r,10*log10(abs(Profile_CR4_Range)/max(abs(Profile_CR4_Range))),'r- ')
%view(0,0) to see range profile
%view(90,0) to see azimuth profile



%% Checking impulse response of corner reflectors, finding azimuth and range profiles for the Corner Reflectors.

Oversample = 128; %oversampling factor
% [Profile_CR4_Azi Profile_CR4_Range peak_pixel_numbers_CR4] = Response30350(slc,18,Oversample,1295:1311,11260:11280);
% title('CR4')
% [Profile_CR5_Azi Profile_CR5_Range peak_pixel_numbers_CR5] = Response30350(slc,19,Oversample,314:330,8090:8110);
% title('CR5')
% [Profile_CR6_Azi Profile_CR6_Range peak_pixel_numbers_CR6] = Response30350(slc,20,Oversample,1132:1142,6070:6090);
% title('CR6')
% [Profile_CR14_Azi Profile_CR14_Range peak_pixel_numbers_CR14] = Response30350(slc,21,Oversample,1190:1210,2410:2422);
% title('CR14')

[Profile_CR4_Azi Profile_CR4_Range peak_pixel_numbers_CR4] = Response30350_2(slc, 18, Oversample, range_pix_loc_2_corr(1), azi_pix_loc_2_corr(1));
[CR4_Range_Resolution CR4_Range_PSLR] = Calc3dBRes(Profile_CR4_Range,Oversample,dr);
[CR4_Azi_Resolution CR4_Azi_PSLR] = Calc3dBRes(Profile_CR4_Azi,Oversample,dx);

[Profile_CR5_Azi Profile_CR5_Range peak_pixel_numbers_CR5] = Response30350_2(slc, 20, Oversample, range_pix_loc_2_corr(2), azi_pix_loc_2_corr(2));
[CR5_Range_Resolution CR5_Range_PSLR] = Calc3dBRes(Profile_CR5_Range,Oversample,dr);
[CR5_Azi_Resolution CR5_Azi_PSLR] = Calc3dBRes(Profile_CR5_Azi,Oversample,dx);

[Profile_CR6_Azi Profile_CR6_Range peak_pixel_numbers_CR6] = Response30350_2(slc, 22, Oversample, range_pix_loc_2_corr(3), azi_pix_loc_2_corr(3));
[CR6_Range_Resolution CR6_Range_PSLR]= Calc3dBRes(Profile_CR6_Range,Oversample,dr);
[CR6_Azi_Resolution CR6_Azi_PSLR] = Calc3dBRes(Profile_CR6_Azi,Oversample,dx);

[Profile_CR14_Azi Profile_CR14_Range peak_pixel_numbers_CR14] = Response30350_2(slc, 24, Oversample, range_pix_loc_2_corr(4), azi_pix_loc_2_corr(4));
[CR14_Range_Resolution CR14_Range_PSLR] = Calc3dBRes(Profile_CR14_Range,Oversample,dr);
[CR14_Azi_Resolution CR14_Azi_PSLR] = Calc3dBRes(Profile_CR14_Azi,Oversample,dx);

CR_Range_Resolutions = [CR4_Range_Resolution, CR5_Range_Resolution, CR6_Range_Resolution, CR14_Range_Resolution];
CR_Azi_Resolutions = [CR4_Azi_Resolution, CR5_Azi_Resolution, CR6_Azi_Resolution, CR14_Azi_Resolution];
Azi_PSLR = [CR4_Azi_PSLR CR5_Azi_PSLR CR6_Azi_PSLR CR14_Azi_PSLR];
Range_PSLR = [CR4_Range_PSLR CR5_Range_PSLR CR6_Range_PSLR CR14_Range_PSLR];

c = physconst('lightspeed');
%The theoretical resolutions
Theo_Range_Res = c/(2*slcParameters2.rangeBandwidth);
Theo_Azi_Res = Vg/slcParameters2.Bp;

%difference between resolution and theoretical resolution:
ResDifferenceRange = CR_Range_Resolutions - Theo_Range_Res; %Outcome: calculated smaller than theoretical, hmm shouldn't be)
ResDifferenceAzi = CR_Azi_Resolutions - Theo_Azi_Res; %Outcome: calculated quite bigger (over 4m) than theoretical.


%comparing the peak sub-pixel location obtained from the integrated
%profiles with the geographical sub-pixel location
azi_pix_loc_FromProfiles = [peak_pixel_numbers_CR4(1), peak_pixel_numbers_CR5(1), peak_pixel_numbers_CR6(1), peak_pixel_numbers_CR14(1)];
range_pix_loc_FromProfiles = [peak_pixel_numbers_CR4(2), peak_pixel_numbers_CR5(2), peak_pixel_numbers_CR6(2), peak_pixel_numbers_CR14(2)];

%pixel difference multiplied with pixel spacing to get geolocation accuracy
%in meters.
GeolocationAccuracyAzi = (azi_pix_loc_FromProfiles - azi_pix_loc_1_corr )*dx;
GeolocationAccuracyRange = (range_pix_loc_1_corr - range_pix_loc_FromProfiles)*dr;
GeolocationAccuracy = sqrt(GeolocationAccuracyAzi.^2 + GeolocationAccuracyRange.^2);

% GeolocationAccuracyAzi_mean = mean(abs(GeolocationAccuracyAzi));
% GeolocationAccuracyAzi_std = std(abs(GeolocationAccuracyAzi));
% GeolocationAccuracyRange_mean = mean(abs(GeolocationAccuracyRange));
% GeolocationAccuracyRange_std = std(abs(GeolocationAccuracyRange));

GeolocationAccuracyAzi_mean = mean(GeolocationAccuracyAzi);
GeolocationAccuracyAzi_std = std(GeolocationAccuracyAzi);
GeolocationAccuracyRange_mean = mean(GeolocationAccuracyRange);
GeolocationAccuracyRange_std = std(GeolocationAccuracyRange);

GeolocationAccuracy_mean = mean(GeolocationAccuracy);
GeolocationAccuracy_std = std(GeolocationAccuracy);

%plot with corrected positions
figure(22); imagesc(abs(slc),[0 1]);
colorbar
colormap(gray)
hold on;
plot(range_pix_loc_2_corr,azi_pix_loc_2_corr,'*y ') %method number 2, Fitting Range curve
plot(range_pix_loc_1_corr,azi_pix_loc_1_corr,'*r ') %method number 1, DopplerEq
plot(sali(:,1)',sali(:,2)','*b ')  %John's function
plot(range_pix_loc_FromProfiles,azi_pix_loc_FromProfiles,'*m ') %subpixel accuracy from interpolated profiles
xlabel('Range pixels')
ylabel('Azimuth pixels')
legend('RangeEq (corrected)', 'DopplerEq (corrected)', 'llh2sali', 'InterpolatedProfiles')
%labeling CR ID numbers to plot
text(range_pix_loc + 375,azi_pix_loc - 25,num2cell(CRgeo(:,1)'),'Fontsize', 15,'Color','yellow')


figure(23); plot(0,0,'*m ')
hold on;
plot(GeolocationAccuracyRange*100,GeolocationAccuracyAzi*100, 'r* '); %times 100 to get units in cm
grid on;
xlim([-100 100])
ylim([-100 100])
xlabel('cm')
ylabel('cm')
title({['Geolocation accuracy of Sentinel-1A IW'], ['\mu = ' num2str(round(GeolocationAccuracy_mean,3)*100) ' cm, \sigma = ' num2str(round(GeolocationAccuracy_std,3)*100) ' cm.'], ['\mu_{azi} = ' num2str(round(GeolocationAccuracyAzi_mean,3)*100) ' cm, \sigma_{azi} = ' num2str(round(GeolocationAccuracyAzi_std,3)*100) ' cm.'], ['\mu_{range} = ' num2str(round(GeolocationAccuracyRange_mean,3)*100) ' cm, \sigma_{range} = ' num2str(round(GeolocationAccuracyRange_std,3)*100) ' cm.']})
%title({['Geolocation accuracy of Corner Reflectors'], ['\mu = ' num2str(round(GeolocationAccuracy_mean,3)*100) ' cm, \sigma = ' num2str(round(GeolocationAccuracy_std,3)*100) ' cm.']})
text(GeolocationAccuracyRange*100 + 2, GeolocationAccuracyAzi*100, num2cell(CRgeo(:,1)'),'Fontsize', 15,'Color','black')


text(GeolocationAccuracyRange*100 + 2, GeolocationAccuracyAzi*100, num2cell(CRgeo(:,1)'),'Fontsize', 15,'Color','black')

