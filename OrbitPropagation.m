%--------------------------------------------------------------------------
%
%  OrbitPropagation.m
%
%  this function propagates an orbit forward in time using ode45 given an intial
%  position and velocity.  This function also plots the orbit
% 
%  inputs:
%    r0       - 3x1 vector of intial postion of satellite in ECI            km
%    v0       - 3x1 vector of intial velocity of satellite in ECI           km/s
%    t0       - Initial time                                                JD
%    t0       - Final   time                                                JD
%
%  outputs:
%    rf       - 3x1 vector of final postion of satellite in ECI             km
%    vf       - 3x1 vector of final velocity of satellite in ECI            km/s
%    oef      - 6x1 vector of final classical orbital elements
%
% Last modified:   10/21/2019   T. Schuler
% 
% -------------------------------------------------------------------------


function [rf, vf, oef] = OrbitPropagation(r0,v0,t0,tf)

%FOR NOW....
%t0 = 0;
%tf = 1;
dt = .01; %intervals
time_span = [t0:dt:tf];

%r0= [-6796; 4025; 3490]; 
%v0 = [-3.7817; -6.0146; 1.1418];

mu= 3.986004254*10^5;    %           Earth's Gravitational Constant
RE = 6378.137;           % km         Earth Radius

ICs = horzcat(r0',v0')';

tolerance = 1e-12;

options = odeset('RelTol',tolerance,'AbsTol',tolerance);

[T,Xt] = ode45(@(t,y)  dynEq(t,y,mu),time_span,ICs,options);

x= Xt(:,1);
y= Xt(:,2);
z= Xt(:,3);
vx= Xt(:,4);
vy= Xt(:,5);
vz= Xt(:,6);

%Extract Last item from iteration, which is the final time stamp
rf = [x(size(x,1),:,1); y(size(y,1),:,1); z(size(z,1),:,1)];
vf = [vx(size(vx,1),:,1); vy(size(vy,1),:,1); vz(size(vz,1),:,1)];

%Determine new orbital elements.  Only True anomaly should change for an
%orbit propagation 

[a,e,i,Omega,omega,f] = OrbitalElements(rf,vf);
oef = [a; e; i; Omega; omega; f];

% %Plot Orbits for time durration
%{
figure(1)
plot3(x,y,z,'Color', rand(1,3),'Linewidth',2);
hold on;
[ex, ey, ez] = sphere(25)
surf(ex*RE,ey*RE,ez*RE);
axis equal
%}

%% Helper function: 2 Body Kepelerian Dynamic Equation

function [dx] = dynEq(T, X , mu)

    x= X(1);    vx=X(4);
    y= X(2);    vy=X(5);
    z= X(3);    vz=X(6);

    r = sqrt(x^2 + y^2 + z^2);

    dx = [vx; vy; vz; -((mu/r^3)*x); -((mu/r^3)*y); -((mu/r^3)*z)];
end


end
