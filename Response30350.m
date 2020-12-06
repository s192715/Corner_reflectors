function [Profile_azi, Profile_range, peak_pixel_number] = Response30350(slc, FigureNumber, oversample, profile_pixels_azi, profile_pixels_r)
%This function generates the best azimuth and range profiles over a given corner
%reflector by calculating many side-by-side range and azimuth profiles and selecting
%only the two with the highest peak magnitude. Also returns a 3-D figure
%showing these profiles. Furthermore, it returns the peak pixel numbers in
%azimuth and range.

%Takes in slc image and figure number (says what number the figure should,
%not really important), an oversampling factor,
%and azimuth and range profile pixels (range of pixels the profiles
%covers).

%The function calculates an azimuth profile at all the range pixel
%specified in the range profile, and a range profile at all the azimuth
%pixels specified in the azimuth profile. It then selects the best azimuth
%and range profiles by comparing the maximum magnitude at each profile.
%The profile that has the highest maximum magnitude is considered the one that
%best extends over the peak response of the Corner Reflector.


%Added: function also returns the sub-pixel numbers of the profile peak.


%doing azi profile
oversample_azi = oversample;
profile_length_azi = length(profile_pixels_azi);
ProfileFirst_azi = profile_pixels_r(1);
ProfileLast_azi = profile_pixels_r(end);
MaximumAbs_azi = zeros(1,length(ProfileFirst_azi:ProfileLast_azi)); %allocating for maximum values
k = 1;
for i = [ProfileFirst_azi:1:ProfileLast_azi] %this loop finds at which range pixel the azimuth profile has the greatest peak magnitude
    Test_profile = interpft(slc(profile_pixels_azi,i),oversample_azi*profile_length_azi); %oversampling factor 16
    MaximumAbs_azi(k) = max(abs(Test_profile));
    k = k+1;
end
At_pixel_azi = ProfileFirst_azi + find(MaximumAbs_azi == max(MaximumAbs_azi)) -1;
Profile_azi = interpft(slc(profile_pixels_azi,At_pixel_azi),oversample_azi*profile_length_azi);


%Now doing range profile
oversample_r = oversample;
ProfileFirst_r = profile_pixels_azi(1);
ProfileLast_r = profile_pixels_azi(end);
profile_length_r = length(profile_pixels_r);
MaximumAbs_r = zeros(1,length(ProfileFirst_r:ProfileLast_r));
k = 1;
for i = [ProfileFirst_r:1:ProfileLast_r] %this loop finds at which range pixel the azimuth profile has the greatest peak magnitude
    Test_profile = interpft(slc(i,profile_pixels_r),oversample_r*profile_length_r); %oversampling factor 16
    MaximumAbs_r(k) = max(abs(Test_profile));
    k = k+1;
end
At_pixel_r = ProfileFirst_r + find(MaximumAbs_r == max(MaximumAbs_r)) -1; %at what range pixel the profile is taken, i.e. where maximum value of all azimuth profiles is found
Profile_range = interpft(slc(At_pixel_r,profile_pixels_r),oversample_r*profile_length_r);


%PlottingResults, in magnitude:
% figure(FigureNumber); plot3(ones(1,oversample_azi*profile_length_azi)*At_pixel_azi, profile_pixels_azi(1):(1/oversample_azi):profile_pixels_azi(end) +1-(1/oversample_azi), abs(Profile_azi),'b- ')
% hold on;
% xlabel('Range pixel')
% ylabel('Azimuth Pixel')
% zlabel('Magnitude')
% grid on;
% plot3(profile_pixels_r(1):(1/oversample_r):profile_pixels_r(end) +1-(1/oversample_r),ones(1,oversample_r*profile_length_r)*At_pixel_r,abs(Profile_range),'r- ')


%in dB:
figure(FigureNumber); plot3(ones(1,oversample_azi*profile_length_azi)*At_pixel_azi, profile_pixels_azi(1):(1/oversample_azi):profile_pixels_azi(end) +1-(1/oversample_azi),10*log10(abs(Profile_azi).^2/max(abs(Profile_azi)).^2),'b- ')
hold on;
xlabel('Range pixel')
ylabel('Azimuth Pixel')
zlabel('Power [dB]')
zlim([-50 0])
grid on;
plot3(profile_pixels_r(1):(1/oversample_r):profile_pixels_r(end) +1-(1/oversample_r),ones(1,oversample_r*profile_length_r)*At_pixel_r,10*log10(abs(Profile_range).^2/max(abs(Profile_range)).^2),'r- ')




[peak_azi index_azi] = findpeaks(10*log10(abs(Profile_azi).^2/max(abs(Profile_azi)).^2));
peak_azipix = index_azi(find(peak_azi== max(peak_azi)))/oversample + profile_pixels_azi(1);
[peak_range index_range] = findpeaks(10*log10(abs(Profile_range).^2/max(abs(Profile_range)).^2));
peak_rangepix = index_range(find(peak_range== max(peak_range)))/oversample + profile_pixels_r(1);


%format [azi range]
peak_pixel_number = [peak_azipix, peak_rangepix];

end

