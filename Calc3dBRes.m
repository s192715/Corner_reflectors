function [Res PSLR] = Calc3dBRes(Profile,Oversample,PixSpacing)
%calculates resolution from a azimuth/range dB profile, also takes in interpolation oversampling factor
%and pixel spacing.

%Finds the element in the profile that has the closest to -3dB value, then
%changes it greatly and does it again to find the element closest to -3dB  on the other side.
%We than know how many elements the -3dB peak spans, we divide by the
%Oversampling factor to get how many pixels the peak spans, and then
%multiply with pixel spacing to get the resolution in meters.


a = 10*log10(abs(Profile).^2/max(abs(Profile)).^2);

%finds the two pixel values that are closest to -3dB
ind = interp1(a,1:length(a),-3,'nearest');
b = a;
b(ind) = 100;
ind2 = interp1(b,1:length(b),-3,'nearest');


[pks loc] = findpeaks(abs(Profile).^2/max(abs(Profile)).^2);
pks = sort(pks);
PSLR = 10*log10(pks(end-1)/pks(end));





%Number of pixels between -3dB beamwidth multiplied with pixel spacing to get resolution in m:
Res = (abs(ind-ind2)/Oversample)*PixSpacing; %Oversample is the oversampling factor used in the interpolation.