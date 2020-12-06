%Daniel's solution for all problems from the beginning. Reading of data
%copied from Bjorn

ITRFCoordinatesFile='Q:\DOKSI\EGYETEM\DTU\3_semester\30350_Remote_sensing\project\data\CRs_ITRF2008';

%Have data in memory
ReflecotrsITRFCoordinates= readmatrix(ITRFCoordinatesFile);
% Calculating corner reflectors in ECEF coordinates

%WGS84 parameters:
f = 1/298.257223563;
a = 6378137;
b = a*(1-f);
e = sqrt((a^2 - b^2)/a^2);
e2 = sqrt((a^2 - b^2)/b^2);


