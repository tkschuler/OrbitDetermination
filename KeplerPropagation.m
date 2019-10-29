%--------------------------------------------------------------------------
%
%  OrbitElements.m
%
%  this function uses the kepeler time of flight information to determine
%  the future true anomaly and the new position and velocity vectors
% 
%  inputs:
%    r2        - 3x1 inital position vector of satellite                    km
%    v2        - 3x1 inital position vector of satellite                    km/s
%    t0        - inital time                                                s
%    tf        - final time                                                 s
%     
%  outputs:
%    rf       - 3x1 vector of final postion of satellite in ECI             km
%    vf       - 3x1 vector of final velocity of satellite in ECI            km/s
%    oef      - 6x1 vector of final classical orbital elements
%
% Last modified:   10/28/2019   T. Schuler
% 
% -------------------------------------------------------------------------
function [rf, vf, oef] = KeplerPropagation(r2,v2,t0,tf)

mu = 398600.4354;                         % [km^3/s^2]

[a,e,i,OMEGA,omega,f] = OrbitalElements(r2,v2);

f = wrapTo360(f);
p = a*(1-e^2);
n = sqrt(mu/a^3); 
T = (tf-t0);

E1 = atan(sqrt((1-e)/(1+e))*tand(f/2))*2; %radians  Convert True Anomaly to Eccentric Anomaly  Eq. 4.33

M1 = E1-e*sin(E1); %Calculate Mean Anomaly from Eccentric Anomaly
M2 = M1+n*T;

%Newton Raphson Method
E2=M2+e*sin(M2)+e^2/2*sin(2*M2); %initial guess
i = 0;
F = E2-e*sin(E2)-M2;

while (abs(F) > 1*10^-9)
    F = E2-e*sin(E2)-M2; %Eq. 4.51 to test if E2 is a solution (or root)
    G = 1-e*cos(E2); %Derrivative of F
    E2 = E2 -(F/G); %Newton's Method
    %fprintf('No. of Iterations : %d\n',i);
    i= i+1;
end

%New True Anomaly
f2 = atan(sqrt((1+e)/(1-e))*tan(E2/2))*2; %convert Eccentric Anomaly back to True Anomaly
f2 = rad2deg(f2); %convert to degrees

r = a*(1-e^2)/(1+e*cosd(f2)); %??? which one to use
v = mu*(2/r-1/a);


%% Propagation using new True Anomaly

%r2 = [-15635 4689 7407]
%v2 = [-4.6954 -2.3777 0.6497]

r0 = norm(r2);
v0 = norm(v2);


%I think we may need to be in the range of 0-360 to get the proper angle
%between? 
f2 = wrapTo360(f2);


r = p/(1+e*cosd(f2)); %propagated true anomaly

lf = 1-(r/p)*(1-cosd(f2-f));
lg = r*norm(r2)*sind(f2-f)/sqrt(p*mu);
lf_dot = sqrt(mu/p)*((1-cosd(f2-f))/p-1/r-1/norm(r2))*tand((f2-f)/2);
lg_dot = 1-norm(r2)/p*(1-cosd(f2-f));

rf = lf*r2 + lg*v2;
vf = lf_dot*r2 + lg_dot*v2;

[a,e,i,Omega,omega,f] = OrbitalElements(rf,vf);
oef = [a; e; i; Omega; omega; f];
end

