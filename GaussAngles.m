%--------------------------------------------------------------------------
%
%  GaussAngles.m
%
%  this function uses the Gauss angles only method to find the kepelerian
%  elements of the orbit
% 
%  inputs:
%    lat      - latitude of site observation                                deg
%    lst      - 3x1 Local Siderial Time for Each Julian Date Observation    deg
%    alt      - altitude of site observation                                km
%    rho      - 3x1 site vector in ECEF                                     km
%    ra       - 3x1 vector of right ascension observations                  deg
%    dec      - 3x1 vector of declanation observations                      deg
%    JD       - 3x1 vector of Julian Date observations  
%    JD_prop  - Future Julian Date to propagate orbit to    
%
%  outputs:
%    r0       - 3x1 vector of intial postion of satellite in ECI            km
%    v0       - 3x1 vector of intial velocity of satellite in ECI           km/s
%    oe0      - 6x1 vector of intial classical orbital elements
%    rf       - 3x1 vector of final postion of satellite in ECI             km
%    vf       - 3x1 vector of final velocity of satellite in ECI            km/s
%    oef      - 6x1 vector of final classical orbital elements
%
% Last modified:   10/23/2019   T. Schuler
% 
% -------------------------------------------------------------------------

function [r0, v0, oe0, rf, vf, oef] = GaussAngles(lat,lst,alt,rho,ra,dec,JD,JD_Prop)

mu=  398600.4354;            % km^3/s^2  Earth's Gravitational Constant
RE = 6378.1366;              % km        Earth Radius
omega_E = 7.2921159e-5;      % rad/s     Earth's intertial Rotation Rate

T1 = (JD(1,1)-JD(2,1))*24*60*60
T3 = (JD(3,1)-JD(2,1))*24*60*60

%Line of site Unit Vectors
for i = 1:3
    L(i,1) = cosd(dec(i,1))*cosd(ra(i,1));
    L(i,2) = cosd(dec(i,1))*sind(ra(i,1));
    L(i,3) = sind(dec(i,1));
end

for i = 1:3
    rsite_eci(i,1) = cosd(lat)*cosd(lst(i,1));
    rsite_eci(i,2) = cosd(lat)*sind(lst(i,1));
    rsite_eci(i,3) = sind(lat);
end
rsite_eci = rsite_eci'*(RE+alt)
L
L = L'

% Parameters for finding middle range magnitude
a1 = T3/(T3-T1)
a3 = -T1/(T3-T1)

a1u = (T3*((T3-T1)^2-T3^2))/(6*(T3-T1))
a3u = (-T1*((T3-T1)^2-T1^2))/(6*(T3-T1))

% Determine Parameters for use in eigth-degree equation
M = inv(L)*rsite_eci;

d1 = M(2,1)*a1-M(2,2)+M(2,3)*a3     %km
d2 = M(2,1)*a1u+M(2,3)*a3u          %km s^2
L2 = L(:,2)
r_site_2 = rsite_eci(:,2)
C = L2'*r_site_2

%Solve for roots of 8th degree polynomial
eq_8 = [1 0 -(d1^2+2*C*d1+r_site_2'*r_site_2) 0 0 -2*mu*(C*d2+d1*d2) 0 0 -mu^2*d2^2]
r2 = roots(eq_8);
r2 = r2(real(r2)>0&imag(r2)==0)     %hopefully there's only one real root.
u= mu/r2^3                          % 1/s^2

%Find ci coefficients
ci_coefficients = [-a1-a1u*u; 1; -a3-a3u*u]*-1

%Initial guess of slant ranges. Solve for rho:
rho = M*(ci_coefficients*-1);
rho(1,1) = rho(1,1)/ci_coefficients(1,1);
rho(2,1) = rho(2,1)/ci_coefficients(2,1);
rho(3,1) = rho(3,1)/ci_coefficients(3,1);
rho

%Initial estimate of ranges to find position vectors
r1 = rho(1,1)*L(:,1)+rsite_eci(:,1)
r2 = rho(2,1)*L(:,2)+rsite_eci(:,2)
r3 = rho(3,1)*L(:,3)+rsite_eci(:,3)

r = [r1 r2 r3]'

%% Use Gibbs method to solve for velocity
%v2 = GIBBS(r1,r2,r3)
v2 = HGIBBS(r1,r2,r3,JD);
r0 = r2;
v0 = v2;

%% Determine classical orbital elements of this orbit
[a,e,i,Omega,omega,f] = OrbitalElements(r2,v2);
oe0 = [a; e; i; Omega; omega; f];

%% Iterate to get better solution

% Test Variables from Vallado
% RE = 6378.1363
% TU = sqrt(RE^3/mu) %time unit page 95 of Vallado
% r1 = [.117100; .970064; 1.324052] * RE;
% r2 = [-.078197; 0.370474; 1.574681] * RE;
% r3 = [-1.431677; -0.347732; 1.321131]* RE;
% v2 = [-.583968 -.510942 -.022494] *RE/TU

%From LEcture 10 slides
%r1 = [8004.7213; 2812.9960; 5408.8835]
%r2 = [6313.3958; 5247.5237; 6467.7250]
%r3 = [5272.0417; 6321.1254; 6810.4754]

%r1 = [8004.7213; 6313.3958; 5272.0417];
%r2 = [2812.9960; 5247.5237; 6321.1254];
%r3 = [5408.8835; 6467.7250; 6810.4754];


%% ITERATION CURRENTLY ISN'T WORKING
h = norm(cross(r2,v2))
p = h^2/mu

df12 = acosd(dot(r1,r2)/(norm(r1)*norm(r2)))
df32 = acosd(dot(r2,r3)/(norm(r2)*norm(r3)))

f1 = 1-(norm(r1)/p)*(1-cosd(df12))
f3 = 1-(norm(r3)/p)*(1-cosd(df32))

g1= (norm(r1)*norm(r2)*sind(df12))/(sqrt(mu*p))
g3= (norm(r3)*norm(r2)*sind(df32))/(sqrt(mu*p))

c1 = g3/(f1*g3-f3*g1)
c3 = -g1/(f1*g3-f3*g1)

%% Orbit Determination
t0 = JD(2,1)
tf = JD_Prop

[rf, vf, oef] = OrbitPropagation(r0, v0, t0, tf)

%% Visualization 
%{
figure
dt = datetime(JD(1,1),'convertfrom','juliandate')
[Y,M,D] = ymd(dt);
[h,m,s] = hms(dt);
utc1 = [Y M D h m s]
o1 = eci2lla(rsite_eci(:,1)',utc1)


dt2 = datetime(JD(2,1),'convertfrom','juliandate')
[Y,M,D] = ymd(dt2);
[h,m,s] = hms(dt2);
utc2 = [Y M D h m s]
o2 = eci2lla(rsite_eci(:,2)',utc2)

dt3 = datetime(JD(3,1),'convertfrom','juliandate')
[Y,M,D] = ymd(dt3);
[h,m,s] = hms(dt3);
utc3 = [Y M D h m s]
o3 = eci2lla(rsite_eci(:,3)',utc3)

%Site Location
geoplot([lat],[-110],'r*')
hold on
geoplot([o1(1) o2(1) o3(1)],[o1(2) o2(2) o3(2)],'g*')
hold on
geoplot([lat lat lat],[lst(1)-360 lst(2)-360 lst(3)-360],'b-*')
hold off
geolimits([45 62],[-180 180])

%}
end
