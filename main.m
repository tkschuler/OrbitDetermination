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
clear
clc
close all;

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

JD_I(1,1) = 2454872.241766892;
JD_I(2,1) = 2454872.241940503;
JD_I(3,1) = 2454872.242114115;

ra_I(1,1) = 30.859159090717;
ra_I(2,1) = 14.564451739639;
ra_I(3,1) = 0.829762748762;

dec_I(1,1) = 79.318796817875;
dec_I(2,1) = 78.120651560859;
dec_I(3,1) = 75.903618501209;

lst_I(1,1) = 295.996368384251;
lst_I(2,1) = 296.059039499724;
lst_I(3,1) = 296.121710615405;


%% C-Data

JD_C(1,1) = 2454871.514010361;
JD_C(2,1) = 2454871.514183972;
JD_C(3,1) = 2454871.514357583;

ra_C(1,1) = 5.931355414284;
ra_C(2,1) = 6.369337583606;
ra_C(3,1) = 6.814572192903;

dec_C(1,1) = -26.399712354399;
dec_C(2,1) = -23.712111094605;
dec_C(3,1) = -20.736087662850;

lst_C(1,1) = 33.286705754698;
lst_C(2,1) = 33.349376869963;
lst_C(3,1) = 33.412047985644;


%% Convert Latitude, Longitude, Altitude to r_site ECEF vector

Re = 6378.1366          %km
lat = 32.2227       	%degrees
lon = -110.0101     	%degrees
alt = 0.757         	%km

lla = [lat lon alt];
%rho = lla2ecef(lla);
rho = [(Re+alt)*cosd(lat)*cosd(lon) (Re+alt)*cosd(lat)*sind(lon) (Re+alt)*sind(lon)]
rho = rho'/1000 


%% Orbit Determination & Visualization
JD_Prop = 2454873.2055555555; %Final Julian Date to Propagate to
[r0, v0, oe0, rf, vf, oef] = GaussAngles(lat,lst_I,alt,rho,ra_I,dec_I,JD_I,JD_Prop)

OrbitViz(r0,v0)

[r0, v0, oe0, rf, vf, oef] = GaussAngles(lat,lst_C,alt,rho,ra_C,dec_C,JD_C,JD_Prop)

OrbitViz(r0,v0)


%% Gauss Distribution

%Do you guys have a better way to do this?

%Add Normal Distribution to computer possible errors
figure(2)

[r0, v0, oe0, rf, vf, oef] = GaussAngles(lat,lst_I,alt,rho,ra_I,dec_I,JD_I,JD_Prop)

t0 = JD_I(2,1)
tf = JD_Prop

r_min = -2.0;
r_max = 2.0;
n = 100;
r_error = r_min+rand(1,n)*(r_max-r_min)

v_min = -2.0;
v_max = 2.0;
n = 100;
v_error = v_min+rand(1,n)*(v_max-v_min)
    
for i = 1:100
    [rf, vf, oef] = OrbitPropagation(r0*r_error(i),v0*v_error(i),t0,tf);
    
    fprintf('I Error point calculation: %i.\n', i);
    
    plot3(rf(1),rf(2),rf(3),'*','Color','r','MarkerSize',4);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    hold on
    grid on
end

%% C- Group

[r0, v0, oe0, rf, vf, oef] = GaussAngles(lat,lst_C,alt,rho,ra_C,dec_C,JD_C,JD_Prop)
    
t0 = JD_C(2,1)
tf = JD_Prop

r_min = -2.0;
r_max = 2.0;
n = 100;
r_error = r_min+rand(1,n)*(r_max-r_min)

v_min = -2.0;
v_max = 2.0;
n = 100;
v_error = v_min+rand(1,n)*(v_max-v_min)
    
for i = 1:100
    [rf, vf, oef] = OrbitPropagation(r0*r_error(i),v0*v_error(i),t0,tf);
    
    fprintf('C Error point calculation: %i.\n', i);
    
    plot3(rf(1),rf(2),rf(3),'*','Color','b','MarkerSize',4);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    hold on
    grid on
end



