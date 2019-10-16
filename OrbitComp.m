
clc
clear all

%Initial Conditions & Observations
mu = 398600.4354            %[km^3/s^2]
Re = 6378.1366              %[km]
omegaE = 7.2921159e-5       %[rad/s]
lat = 40%32.248814             %site Latitude [deg]
lst = -110%110.987419            %site Longitude [deg]
alt = 2000%757                   %site altitude [m]
ra = [.939913 45.025748 67.886655];%[30.859159090717 14.564451739639 0.829762748762];          %right ascension [deg]
dec = [18.667717 35.664741 36.996583];%[79.318796817875 78.120651560859 75.903618501209];        %declination [deg]
JD = [2454872.241766892 2454872.241940503 2454872.242114115];   %julian dates

%Change in time between observations
T1 = -8*60%(JD(1,1)-JD(1,2))*60                               %delta t btwn OBS 1&2 [sec]
T3 = 4*60%(JD(1,3)-JD(1,2))*60                               %delta t btwn OBS 2&3 [sec]

%Time differences
a1 = T3/(T3-T1)
a3 = -T1/(T3-T1)

%Change in times
a1u = (T3*((T3-T1)^2-T3^2))/(6*(T3-T1))
a2u = -(T1*((T3-T1)^2-T1^2))/(6*(T3-T1))

%Line Of Sight (LOS) unit vector to satellite observations
for i=1:3
    L(i,1) = cosd(dec(1,i))*cosd(ra(1,i));
    L(i,2) = cosd(dec(1,i))*sind(ra(1,i));
    L(i,3) = sind(dec(1,i));
end
L=L'

%Inverse of the unit vector
L_inv = inv(L)

%Convert Site Latitude, Longitude, and Altitude to ECEF (SEZ) frame
%lla = [40 -110 2000];%[32.248814 -110.987419 757];
%rsiteECEF = lla2ecef(lla);                  %ECEF Vector
%rsiteECEF = rsiteECEF'
rsiteSEZ = [Re*cosd(lat)*cosd(lst) Re*cosd(lat)*sind(lst) sind(lat)]'
%not sure if this should be IJK or SEZ

%Calculate M (mean anomaly?)
M = L_inv*rsiteIJK
%not sure if M uses SEZ site or IJK site

%Slant range position vector
rhoSEZ = [-L(1,1)*cosd(dec(1,1))*cosd(ra(1,1)); 
L(1,1)*cosd(dec(1,1))*sind(ra(1,1));
L(1,1)*sind(dec(1,1))]

%Transformation matrix from SEZ to IJK


%TN DCM
gamma = 110
phi = 40

NT = [-sind(gamma) cosd(gamma) 0;
    -cosd(gamma)*sind(phi) -sind(gamma)*sind(phi) cosd(phi);
    cosd(gamma)*cosd(phi) sind(gamma)*cosd(phi) sind(phi)]'

