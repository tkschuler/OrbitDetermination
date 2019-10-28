% Alternative method to propagate the Orbit
%
function [rf_k, vf_k, oef_k] = KeplerPropagation(r0,v0,t0,tf,oe0)

mu = 3.986e5;                           % [km^3/s^2]
f = oe0(6,1);
a = oe0(1,1);
e = oe0(2,1);
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
f2 = theta2*(180/pi);                   % [°]
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