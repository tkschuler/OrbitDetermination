%--------------------------------------------------------------------------
%
%  main.m
%
%  this script takes 3 satellite observations from a site on earth and 
%  and determines the orbit of the satellite using Gauss
%
% Last modified:   10/28/2019   J. Hammes
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

Re = 6378.1366;         %km
lat = 32.2227;       	%degrees
lon = -110.0101;     	%degrees
alt = 0.757;         	%km

%lla = [lat lon alt];
%rho = lla2ecef(lla);
rho = [(Re+alt)*cosd(lat)*cosd(lon) (Re+alt)*cosd(lat)*sind(lon) (Re+alt)*sind(lon)];
rho = rho'/1000; 


%% Orbit Determination & Visualization
JD_Prop = 2454873.2055555555; %Final Julian Date to Propagate to
[r0, v0, oe0, rf, vf, oef] = GaussAngles(lat,lst_I,alt,rho,ra_I,dec_I,JD_I,JD_Prop)
rfI = rf;
dt = datetime(JD_I(2,1),'convertfrom','juliandate');
[Y,M,D] = ymd(dt);
[h,m,s] = hms(dt);
utc1 = [Y M D h m s];
o1_I = eci2lla(r0',utc1);

dt2 = datetime(JD_Prop,'convertfrom','juliandate');
[Y,M,D] = ymd(dt2);
[h,m,s] = hms(dt2);
utc2 = [Y M D h m s];
o2_I = eci2lla(rf',utc2);


OrbitViz(r0,v0)

[r0, v0, oe0, rf, vf, oef] = GaussAngles(lat,lst_C,alt,rho,ra_C,dec_C,JD_C,JD_Prop)
rfC = rf;
dt = datetime(JD_C(2,1),'convertfrom','juliandate');
[Y,M,D] = ymd(dt);
[h,m,s] = hms(dt);
utc1 = [Y M D h m s];
o1_C = eci2lla(r0',utc1);

dt2 = datetime(JD_Prop,'convertfrom','juliandate');
[Y,M,D] = ymd(dt2);
[h,m,s] = hms(dt2);
utc2 = [Y M D h m s];
o2_C = eci2lla(rf',utc2);

OrbitViz(r0,v0)

figure(2)
geoplot([lat],[lon],'r*')
hold on
geoplot([o1_I(1) o2_I(1)],[o1_I(2) o2_I(2)],'g-*')
hold on
geoplot([o1_C(1) o2_C(1)],[o1_C(2) o2_C(2)],'b-*')
hold off
NominalDistance = norm(rfI-rfC)

%% Gauss Distribution

RI = zeros(3,100);
RC = zeros(3,100);

%% ODE Propagation
%Add Normal Distribution to computer possible errors
figure(3)

[r0, v0, oe0, rf, vf, oef] = GaussAngles(lat,lst_I,alt,rho,ra_I,dec_I,JD_I,JD_Prop)

t0 = JD_I(2,1)*24*60*60;
tf = JD_Prop*24*60*60

%tf= juliandate(datetime('2009-02-10 17:15:59'))*24*60*60

fprintf('Calculating I Error Distribution...\n');
    
for i = 1:100
    
    p = 3; %3 sigma

    r_min = -2.0;
    r_max = 2.0;
    n=3; %size of array to generate

    r_error = r_min + (r_max - r_min)*sum(rand(n,p),1)/p;
    
    v_min = -0.1;
    v_max = 0.1;
    n = 3;
    v_error = v_min+(v_max-v_min)*sum(rand(n,p),1)/p;
    
    
    [rf, vf, oef] = KeplerPropagation(r0'+r_error,v0'+v_error,t0,tf);
    %[rf, vf, oef] = OrbitPropagation(r0+r_error',v0+v_error',t0,tf);
    RI(:,i)=rf;
    
    plot3(rf(1),rf(2),rf(3),'*','Color','r','MarkerSize',4);
    hold on
    %grid on
    fprintf('Calculating I distribution point... %d \n', i);
end


%% C- Group

[r0, v0, oe0, rf, vf, oef] = GaussAngles(lat,lst_C,alt,rho,ra_C,dec_C,JD_C,JD_Prop)
    
t0 = JD_C(2,1)*24*60*60;
tf = JD_Prop*24*60*60;

fprintf('Calculating C Error Distribution... \n');

for i = 1:100
    
    p = 3; %3 sigma

    r_min = -2.0;
    r_max = 2.0;
    n=3; %size of array to generate

    r_error = r_min + (r_max - r_min)*sum(rand(n,p),1)/p;
    
    v_min = -0.1;
    v_max = 0.1;
    n = 3;
    v_error = v_min+(v_max-v_min)*sum(rand(n,p),1)/p;
    
    [rf, vf, oef] = KeplerPropagation(r0'+r_error,v0'+v_error,t0,tf);
    %[rf, vf, oef] = OrbitPropagation(r0+r_error',v0+v_error',t0,tf);
    RC(:,i)=rf;
    
    plot3(rf(1),rf(2),rf(3),'*','Color','b','MarkerSize',4);
    hold on
    grid on
    fprintf('Calculating C distribution point... %d \n', i);
end

%Why is the legend not doing the right color?
legend({'I-Orbit','C-Orbit'})
xlabel('I distance (km)')
ylabel('J distance (km)')
zlabel('K distance (km)')
title('3 Sigma Distribution of Two Overlapping Orbits')

%Better for Visualization
hold off

%Matrix contains all the distances

AM= zeros(100,100);
for i=1:100
    for c=1:100
        AM(c,i)=norm(RI(:,i)-RC(:,c));
    end
end
minAM = min(AM,[],[1 2]);
maxAM = max(AM,[],[1 2]);
app_range = (maxAM-minAM)/50;                        %gives the intervals
app_num = zeros(1,51);
for n=0:50
    interval(n+1,1) = minAM+n*app_range;
    for i=1:100
        for c=1:100
            if minAM+n*app_range <= AM(i,c) && AM(i,c)< minAM+(n+1)*app_range
                app_num(1,n+1) = app_num(1,n+1)+1;
            end
        end
    end
end
figure(4)
bar(interval,app_num);
xlabel('Approach distance (km)')
ylabel('Number of approaches')
title('Distribution of approach distance for 100 samples per orbit each')
