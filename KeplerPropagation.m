% Alternative method to propagate the Orbit
%

%{
function [rf_k, vf_k, oef_k] = KeplerPropagation(r0,v0,t0,tf)

mu = 3.986e5;                           % [km^3/s^2]
[a,e,i,OMEGA,omega,f] = OrbitalElements(r0,v0)
p = a*(1-e^2);
theta1 = f*(pi/180);                    % [rad]
ft = tf-t0;                             % [s]
n = sqrt(mu/a^3);               
E1 = atan(sqrt((1-e)/(1+e))*tan(theta1/2))*2; 
M1 = E1 - e*sin(E1);                    % Mean anomaly 1          
M2 = M1+n*ft;                           % Mean anomaly 2

%Newton Raphson for E2
i = 1;
E2(i) = M2;
f_E2(i) = E2(i) - e*sin(E2(i)) - M2;

while abs(f_E2(i)) > 1e-12
    E2(i+1)= E2(i) - f_E2(i)/(1-e*cos(E2(i)));
    i = i + 1;
    f_E2(i) = E2(i) - e*sin(E2(i)) - M2;    
end

E2 = E2(i);                             % [rad]
theta2 = acos((cos(E2)-e)/(1-e*cos(E2)));
f2 = theta2*(180/pi)                   % [�]
%gamma2 = atan((e*sin(theta2))/(1+e*cos(theta2)))*(180/pi);
r2 = (a*(1-e^2))/(1+e*cos(theta2));     % [km] NOT THE VECTOR
v2 = mu*(2/r2 - 1/a);                   % [km/s] NOT THE VECTOR
% Lagrangian coefficients, see Kluever p.133

lf = 1-(r2/p)*(1-cosd(f2-f));
lg = r2*norm(r0)*sind(f2-f)/sqrt(p*mu);
lf_dot = sqrt(mu/p)*((1-cosd(f2-f))/p-1/r2-1/norm(r0))*tand((f2-f)/2);
lg_dot = 1-norm(r0)/p*(1-cosd(f2-f));

rf_k = lf*r0 + lg*v0;
vf_k = lf_dot*r0 + lg_dot*v0;

[a,e,i,Omega,omega,f] = OrbitalElements(rf_k,vf_k);
oef_k = [a; e; i; Omega; omega; f];
end

%}



% Alternative method to propagate the Orbit
%
function [rf_k, vf_k, oef_k] = KeplerPropagation(r2,v2,t0,tf)

mu = 398600.4354;                         % [km^3/s^2]

[a,e,i,OMEGA,omega,f] = OrbitalElements(r2,v2);

f = wrapTo360(f);
p = a*(1-e^2);
n = sqrt(mu/a^3); 
T = (tf-t0)

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

rf_k = lf*r2 + lg*v2;
vf_k = lf_dot*r2 + lg_dot*v2;

[a,e,i,Omega,omega,f] = OrbitalElements(rf_k,vf_k);
oef_k = [a; e; i; Omega; omega; f];
end

