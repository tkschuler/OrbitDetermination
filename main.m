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
% Last modified:   10/21/2019   T. Schuler
% 
% -------------------------------------------------------------------------

format long;
format compact;
clear;
clc;

%% Observations

%{
% Example from Rosengren's Slides:

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

%Find Local Siderial Time for each Observation
for i = 1:3
    [GST, LST] = siderial_time(JD(i,1),lon)
    lst(i,1) = LST; %deg
end

%}

%% I-data

lat = 32.248814        	%degrees
lon = 110.987419      	%degrees
alt = 0.757         	%km

JD(1,1) = 2454872.241766892;
JD(2,1) = 2454872.241940503;
JD(3,1) = 2454872.242114115;

ra(1,1) = 30.859159090717;
ra(2,1) = 14.564451739639;
ra(3,1) = 0.829762748762;

dec(1,1) = 79.318796817875;
dec(2,1) = 78.120651560859;
dec(3,1) = 75.903618501209;

lst(1,1) = 295.996368384251;
lst(2,1) = 296.059039499724;
lst(3,1) = 296.121710615405;

tf = 2454873.2055555555; %Final Julian Date to Propagate to
TOF = tf - JD(2,1); 


%% C-Data
%{
JD(1,1) = 2454871.514010361;
JD(2,1) = 2454871.514183972;
JD(3,1) = 2454871.514357583;

ra(1,1) = 5.931355414284;
ra(2,1) = 6.369337583606;
ra(3,1) = 6.814572192903;

dec(1,1) = -26.399712354399;
dec(2,1) = -23.712111094605;
dec(3,1) = -20.736087662850;

lst(1,1) = 33.286705754698;
lst(2,1) = 33.349376869963;
lst(3,1) = 33.412047985644;

%}

%% Convert Latitude, Longitude, Altitude to r_site ECEF vector
lla = [lat lon alt];
rho = lla2ecef(lla); 
rho = rho'/1000 

%% Orbit Determination
[r0, v0, oe0, rf, vf, oef] = GaussAngles(lat,lst,rho,ra,dec,JD,TOF)

