%% Function to Determine Orbital Parameters
function [a,e,i,omega,w,theta] = OrbitalElements(r,v)

mu = 3.986*10^5; %Gravitational Constant

E = norm(v)^2/2 - mu/norm(r)

a = -mu/(2*E)

e_vec= 1/mu*((norm(v)^2-mu/norm(r))*r -dot(r,v)*v)
e = norm(e_vec)

%unit vectors
K= [0;0;1];
I = [1;0;0];
J= [0;1;0]

h= cross(r,v)

i = acosd(dot(K,h)/norm(h))

n = cross(K,h)

cos_omega = dot(I,n)/norm(n)
sin_omega = dot(J,n)/norm(n) %Double check Quadrant
omega = atan2d(sin_omega,cos_omega)

w = acosd(dot(n,e_vec)/(norm(n)*e))

if e_vec(3) < 0
    w = 360-w
end

if i == 0
    w = w+omega %Longitude of Periapsis
end

theta = acosd(dot(e_vec,r)/(e*norm(r))) 

if dot(r,v) < 0
    theta = 360-theta
end
    
end