clear;
clc;

mu= 3.986004254*10^5 %Earth's Gravitational Constant
ER = 6378.137; %Earth Radius

%Tau, time diferences between observations
T1 = -8
T3 = 4

%Time differences
a1 = T3/(T3-T1)
a3 = -T1/(T3-T1)

%convert to seconds?
T1 = T1*60;
T3= T3 *60;

a1u = (T3*((T3-T1)^2-T3^2))/(6*(T3-T1))
a3u = (-T1*((T3-T1)^2-T1^2))/(6*(T3-T1))


%Observations
%First column is right ascention, second column is declanation
Obs =  [.939913 18.667717;
        45.025748 35.664741;
        67.886655 36.996583]

    %From Vallado Example on p. 401
%Obs =  [118.6780872 27.5777952;
%        162.5580465 30.0605002;
%        187.7907123 17.081080]

%Line of site Unit Vectors
for i = 1:size(Obs,1)
    L(i,1) = cosd(Obs(i,2))*cosd(Obs(i,1));
    L(i,2) = cosd(Obs(i,2))*sind(Obs(i,1));
    L(i,3) = sind(Obs(i,2));
end

L = L'

%Convert Latitude, Longitude, Altitude to ECEF vector
lla = [40 -110 2000];
rho = lla2ecef(lla); %ECEF Vector  May not be able to use this
rho = rho'/1000 %km


%% Find ECI

%clc;

rho_SEZ = [.073515; -.035223; .048402]

%t1 = datetime('1992-08-20 12:14:02');
t1 = datetime('1995-05-20 03:17:02');
format longG
jd = juliandate(t1)

gd = 40;

[GST, lst] = siderial_time(jd,-110)


SEZ2IJK = [sind(gd)*cosd(lst) -sind(lst) cosd(gd)*cosd(lst);
       sind(gd)*sind(lst)  cosd(lst) cosd(gd)*sind(lst);
       -cosd(gd)          0        sind(gd)]
   

%hard coded for now-------------------
r_site_ECI = [4054.881 3956.224 3905.073;
          2748.195 2888.232 2956.935;
          4074.237 4074.364 4074.430]
      
o1 = SEZ2IJK*rho_SEZ
      
L_inv = inv(L)

M = L_inv*r_site_ECI

d1 = M(2,1)*a1-M(2,2)+M(2,3)*a3 %km
d2 = M(2,1)*a1u-M(2,2)+M(2,3)*a3u

L2 = L(:,2)
r_site_2 = r_site_ECI(:,2)

C = L2'*r_site_2


eq_8 = [1 0 -(d1^2+2*C*d1+r_site_2'*r_site_2) 0 0 -2*mu*(C*d2+d1*d2) 0 0 -mu^2*d2^2] %Solve for roots of 8th degree polynomial
r2 = roots(eq_8);

r2 = r2(real(r2)>0&imag(r2)==0) %hopefully there's only one real root.

u= mu/r2^3

%Find ci coefficients
ci_coefficients = [-a1-a1u*u; 1; -a3-a3u*u]*-1

%Initial guess of slant ranges. Solve for rho:
rho = M*(ci_coefficients*-1);
rho(1,1) = rho(1,1)/ci_coefficients(1,1);
rho(2,1) = rho(2,1)/ci_coefficients(2,1);
rho(3,1) = rho(3,1)/ci_coefficients(3,1);
rho

%Initial estimate of ranges to find position vectors
r1 = rho(1,1)*L(:,1)+r_site_ECI(:,1)
r2 = rho(2,1)*L(:,2)+r_site_ECI(:,2)
r3 = rho(3,1)*L(:,3)+r_site_ECI(:,3)

r = [r1 r2 r3]'


%% GIBBS TEST

%r1 = [0 0 6378.137]/ER;
%r2 = [0 -4464.696 -5102.509]/ER;
%r3 = [0 5740.323 3189.068]/ER;

v2 = GIBBS(r1,r2,r3)

[a,e,i,omega,w,theta] = OrbitalElements(r2,v2)


