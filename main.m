%--------------------------------------------------------------------------
%
%  main.m
%
%  this script takes 3 satellite observations from a site on Earth and 
%  and determines the orbit of the satellite using the Gauss Method
%
% Last modified:   10/28/2019   T. Schuler
% 
% -------------------------------------------------------------------------

format long;
format compact;
clear;
clc;
close all;

%% Observations
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

Re = 6378.1366;                     %km
lat = 32.2227;                      %degrees
lon = -110.0101;                    %degrees
alt = 0.757;                        %km
JD_Prop = 2454873.2055555555;       %Final Julian Date to Propagate to

%Determine r_site ecef
rho = [(Re+alt)*cosd(lat)*cosd(lon) (Re+alt)*cosd(lat)*sind(lon) (Re+alt)*sind(lon)];
rho = rho'/1000; 

%% Orbit Determination & Visualization

cprintf('*red','I Group Orbit Determination and Propagation\n\n');
[r0_I, v0_I, oe0_I, rf_I, vf_I, oef_I] = GaussAngles(lat,lst_I,alt,rho,ra_I,dec_I,JD_I,JD_Prop);
OrbitViz(r0_I,v0_I)
hold on

cprintf('*blue','C Group Orbit Determination and Propagation\n\n');
[r0_C, v0_C, oe0_C, rf_C, vf_C, oef_C] = GaussAngles(lat,lst_C,alt,rho,ra_C,dec_C,JD_C,JD_Prop);
OrbitViz(r0_C,v0_C)
hold off;

nominal_distance = norm(rf_I-rf_C);
fprintf('The distance between the satellites at the final epoch is %g [km]\n',nominal_distance)

%% Gauss 3-Sigma Orbit Propagation (n=100)
%Comment out the line to do either Kepler propagation or ode45 propagation
figure(2)
tf = JD_Prop*24*60*60;

cprintf('*green', 'Gauss 3-Sigma Orbit Propagation (n=100) \n');
cprintf('*black','Calculating I Error Distribution...\n');
t0 = JD_I(2,1)*24*60*60;
    
for i = 1:100
    %3-Sigma Normal distribution
    p = 3; %3-sigma
    r_min = -2.0;
    r_max = 2.0;
    n=3; %size of array to generate
    r_error = r_min + (r_max - r_min)*sum(rand(n,p),1)/p;
    
    v_min = -0.1;
    v_max = 0.1;
    n = 3;
    v_error = v_min+(v_max-v_min)*sum(rand(n,p),1)/p;

    %Kepler Propagation
    [rf, vf, oef] = KeplerPropagation(r0_I'+r_error,v0_I'+v_error,t0,tf);
    
    %ode45 Propagation
    %[rf, vf, oef] = ode45_Propagation(r0_I+r_error',v0_I+v_error',t0,tf);
    
    RI(:,i)=rf;
    h1 = plot3(rf(1),rf(2),rf(3),'*','Color','r','MarkerSize',4);
    hold on
    %grid on
    fprintf('%d ', i);
end
fprintf('\n');


%% C- Group
    
t0 = JD_C(2,1)*24*60*60;
cprintf('*black','Calculating C Error Distribution...\n');

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
    
    %Kepler Propagation
    [rf, vf, oef] = KeplerPropagation(r0_C'+r_error,v0_C'+v_error,t0,tf);
    
    %ode45 Propagation
    %[rf, vf, oef] = ode45_Propagation(r0_C+r_error',v0_C+v_error',t0,tf);
    
    RC(:,i)=rf;
    h2 = plot3(rf(1),rf(2),rf(3),'*','Color','b','MarkerSize',4);
    hold on
    grid on
    fprintf('%d ', i);
end
fprintf('\n');

legend([h1,h2],'I-Orbit','C-Orbit')
xlabel('I distance (km)')
ylabel('J distance (km)')
zlabel('K distance (km)')
title('3 Sigma Distribution of Two Overlapping Orbits')

hold off

%% 3-Sigma Histogram of Closest Approach 

AM= zeros(100,100);
for i=1:100
    for c=1:100
        AM(c,i)=norm(RI(:,i)-RC(:,c));
    end
end
minAM = min(AM,[],[1 2]);
maxAM = max(AM,[],[1 2]);
fprintf('The minimal found approach distance is: %g [km]',minAM)
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

figure(3)
bar(interval,app_num);
xlabel('Approach distance (km)')
ylabel('Number of approaches')
title('Distribution of closest approach distance for 100 samples per orbit each')
fprintf('\n');
