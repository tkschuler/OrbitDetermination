%--------------------------------------------------------------------------
%
%  main.m
%
%  this script takes 3 satellite observations from a site on earth and 
%  and determines the orbit of the satellite using Gauss
% 
%  TODO:
%   -Iterate for Gauss to find better solution
%   -Add TOF functionality
%   -Do the same method for Laplace and see if we get the same solutions
%
% Last modified:   10/15/2019   T. Schuler
% 
% -------------------------------------------------------------------------


format long;
format compact;
clear;
clc;

%% Observations

lat = 40        %degrees
lon = -110      %degrees
alt = 2         %km

%Julian Date of observations
JD(1,1) = juliandate(datetime('2012-08-20 11:40:28'));
JD(2,1) = juliandate(datetime('2012-08-20 11:48:28'));
JD(3,1) = juliandate(datetime('2012-08-20 11:52:28'));

% right ascension angle for each observation
ra(1,1) = 0.939913;
ra(2,1) = 45.025748;
ra(3,1) = 67.886655;

% Angle of Declanation for each observation
dec(1,1) = 18.667717;
dec(2,1) = 35.664741;
dec(3,1) = 36.996583;

%Convert Latitude, Longitude, Altitude to r_site ECEF vector
lla = [lat lon alt];
rho = lla2ecef(lla); 
rho = rho'/1000 

%Find Local Siderial Time for each Observation
for i = 1:3
    [GST, LST] = siderial_time(JD(i,1),lon)
    lst(i,1) = LST; %deg
end

TOF = 5 %time of flight

%% Orbit Determination

[r0, v0, oe0, rf, vf, oef] = OrbitDetermination(lat,lst,rho,ra,dec,JD,TOF)
