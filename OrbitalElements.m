%--------------------------------------------------------------------------
%
%  OrbitElements.m
%
%  this function calculates the 6 classical keplerian orbital elements from
%  and initial position and velocity in the ijk frame
% 
%  inputs:
%    r        - 3x1 position vector of satellite            km
%    v        - 3x1 position vector of satellite            km/s
%
%  outputs:
%    a        - semi-major axis                             km
%    e        - eccentricity           
%    i        - inclination                                 deg
%    OMEGA    - Longitude of the Ascending Node             deg
%    omega    - Argument of Periapsis                       deg
%    f        - True Anomaly                                deg
%
% Last modified:   10/28/2019   T. Schuler
% 
% -------------------------------------------------------------------------
function [a,e,i,OMEGA,omega,f] = OrbitalElements(r,v)

    mu = 3.986*10^5; %Gravitational Constant
    E = norm(v)^2/2 - mu/norm(r);
    a = -mu/(2*E);
    e_vec= 1/mu*((norm(v)^2-mu/norm(r))*r -dot(r,v)*v);
    e = norm(e_vec);

    %unit vectors
    K= [0;0;1];
    I = [1;0;0];
    J= [0;1;0];

    h= cross(r,v);
    i = acosd(dot(K,h)/norm(h));
    n = cross(K,h);

    cos_omega = dot(I,n)/norm(n);
    sin_omega = dot(J,n)/norm(n); %Double check Quadrant
    OMEGA = atan2d(sin_omega,cos_omega);
    omega = acosd(dot(n,e_vec)/(norm(n)*e));

    if e_vec(3) < 0
        omega = 360-omega;
    end

    if i == 0
        omega = omega+OMEGA; %Longitude of Periapsis
    end

    f = acosd(dot(e_vec,r)/(e*norm(r)));
    if dot(r,v) < 0
        f = 360-f;
    end 
end