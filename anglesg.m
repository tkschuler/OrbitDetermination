%--------------------------------------------------------------------------
%
%  anglesg.m
%
%  this function solves the problem of orbit determination using three
%  optical sightings.
% 
%  inputs:
%    az1      - azimuth at t1               rad
%    az2      - azimuth at t2               rad
%    az3      - azimuth at t3               rad
%    el1      - elevation at t1             rad
%    el2      - elevation at t2             rad
%    el3      - elevation at t3             rad
%    Mjd1     - Modified julian date of t1
%    Mjd2     - Modified julian date of t2
%    Mjd3     - Modified julian date of t3
%    Rs1      - ijk site1 position vector   m
%    Rs2      - ijk site2 position vector   m
%    Rs3      - ijk site3 position vector   m
%
%  outputs:
%    r        - ijk position vector at t2   m
%    v        - ijk velocity vector at t2   m/s
%
% Last modified:   2015/08/12   M. Mahooti
% 
% -------------------------------------------------------------------------
function [r2, v2] = anglesg ( az1,az2,az3,el1,el2,el3,Mjd1,Mjd2,Mjd3,Rs1,Rs2,Rs3 )
global eopdata
SAT_Const
L1 = [cos(el1)*sin(az1); cos(el1)*cos(az1); sin(el1)]
L2 = [cos(el2)*sin(az2); cos(el2)*cos(az2); sin(el2)]
L3 = [cos(el3)*sin(az3); cos(el3)*cos(az3); sin(el3)]
[lon1, lat1, h1] = Geodetic(Rs1);
[lon2, lat2, h2] = Geodetic(Rs2);
[lon3, lat3, h3] = Geodetic(Rs3);
M1 = LTC(lon1, lat1);
M2 = LTC(lon2, lat2);
M3 = LTC(lon3, lat3);
% body-fixed system
Lb1 = M1'*L1;
Lb2 = M1'*L2;
Lb3 = M1'*L3;
% mean of date system (J2000)
Mjd_UTC = Mjd1;
[UT1_UTC, TAI_UTC, x_pole, y_pole] = IERS(eopdata, Mjd_UTC);
[UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);
Mjd_TT = Mjd_UTC + TT_UTC/86400;
Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
P = PrecMatrix(MJD_J2000,Mjd_TT);
N = NutMatrix(Mjd_TT);
T = N * P;
E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;
Lm1 = E'*Lb1;
Rs1 = E'*Rs1;
Mjd_UTC = Mjd2;
[UT1_UTC, TAI_UTC, x_pole, y_pole] = IERS(eopdata, Mjd_UTC);
[UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);
Mjd_TT = Mjd_UTC + TT_UTC/86400;
Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
P = PrecMatrix(MJD_J2000,Mjd_TT);
N = NutMatrix(Mjd_TT);
T = N * P;
E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;
Lm2 = E'*Lb2;
Rs2 = E'*Rs2;
Mjd_UTC = Mjd3;
[UT1_UTC, TAI_UTC, x_pole, y_pole] = IERS(eopdata, Mjd_UTC);
[UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);
Mjd_TT = Mjd_UTC + TT_UTC/86400;
Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
P = PrecMatrix(MJD_J2000,Mjd_TT);
N = NutMatrix(Mjd_TT);
T = N * P;
E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;
Lm3 = E'*Lb3;
Rs3 = E'*Rs3;
% geocentric inertial position
tau1 = (Mjd1-Mjd2)*86400;
tau3 = (Mjd3-Mjd2)*86400;
a1 = tau3/(tau3-tau1);
a3 =-tau1/(tau3-tau1);
b1 = tau3/(6*(tau3-tau1))*((tau3-tau1)^2-tau3^2);
b3 =-tau1/(6*(tau3-tau1))*((tau3-tau1)^2-tau1^2);
D = inv([Lm1,Lm2,Lm3])*[Rs1,Rs2,Rs2];
d1s = D(2,1)*a1-D(2,2)+D(2,3)*a3;
d2s = D(2,1)*b1+D(2,3)*b3;
Ccye = 2*dot(Lm2,Rs2);
poly(1)=  1.0;  % R2^8... polynomial
poly(2)=  0.0;
poly(3)=  -(d1s^2 + d1s*Ccye + (norm(Rs2))^2);
poly(4)=  0.0;
poly(5)=  0.0;
poly(6)=  -GM_Earth*(d2s*Ccye + 2*d1s*d2s);
poly(7)=  0.0;
poly(8)=  0.0;
poly(9)=  -GM_Earth^2*d2s^2;
rootarr = roots( poly );
bigr2= -99999990.0;
for j=1:8
    if ( rootarr(j) > bigr2 ) & ( isreal(rootarr(j)) )
        bigr2= rootarr(j);
    end  
end
u = GM_Earth/(bigr2^3);
C1 = a1+b1*u;
C2 = -1;
C3 = a3+b3*u;
temp = -D*[C1 C2 C3]';
rho1 = temp(1)/(a1+b1*u);
rho2 = -temp(2);
rho3 = temp(3)/(a3+b3*u);
rhoold1 = rho1;
rhoold2 = rho2;
rhoold3 = rho3;
rho2 = 99999999.9;
ll   = 0;
while ((abs(rhoold2-rho2) > 1e-12) && (ll <= 0 ))
    ll = ll + 1;
    rho2 = rhoold2;
    
    r1 = Rs1+rho1*Lm1;
    r2 = Rs2+rho2*Lm2;
    r3 = Rs3+rho3*Lm3;
    
    magr1 = norm(r1);
    magr2 = norm(r2);
    magr3 = norm(r3);
    
    [v2, theta,theta1,copa,error] = gibbs(r1,r2,r3);
    
    if ( (error ~= '          ok') & (copa < pi/180) )        
        [v2,theta,theta1,copa,error] = hgibbs(r1,r2,r3,Mjd1,Mjd2,Mjd3);
    end
    
    [p, a, e, i, Omega, omega, M] = elements ([r2,v2]);
    
    if ( ll <= 8 )
        u = GM_Earth/magr2^3;
        rdot= dot(r2,v2)/magr2;
        udot= (-3*GM_Earth*rdot)/(magr2^4);
        
        tausqr= tau1*tau1;
        f1=  1 - 0.5*u*tausqr -(1/6)*udot*tausqr*tau1 ...
            - (1/24) * u*u*tausqr*tausqr ...
            - (1/30)*u*udot*tausqr*tausqr*tau1;
        g1= tau1 - (1/6)*u*tau1*tausqr - (1/12) * udot*tausqr*tausqr ...
            - (1/120)*u*u*tausqr*tausqr*tau1 ...
            - (1/120)*u*udot*tausqr*tausqr*tausqr;
        tausqr= tau3*tau3;
        f3=  1 - 0.5*u*tausqr -(1/6)*udot*tausqr*tau3 ...
            - (1/24) * u*u*tausqr*tausqr ...
            - (1/30)*u*udot*tausqr*tausqr*tau3;
        g3= tau3 - (1/6)*u*tau3*tausqr - (1/12) * udot*tausqr*tausqr ...
            - (1/120)*u*u*tausqr*tausqr*tau3 ...
            - (1/120)*u*udot*tausqr*tausqr*tausqr;
    else
        
        theta  = angl( r1,r2 );
        theta1 = angl( r2,r3 );
        
        f1= 1 - ( (magr1*(1 - cos(theta)) / p ) );
        g1= ( magr1*magr2*sin(-theta) ) / sqrt( p );
        f3= 1 - ( (magr3*(1 - cos(theta1)) / p ) );
        g3= ( magr3*magr2*sin(theta1) ) / sqrt( p );
    end
    
    C1 = g3/(f1*g3-f3*g1);
    C2 = -1;
    C3 =-g1/(f1*g3-f3*g1);
    
    H1 = GM_Earth*tau3/12;
    H3 =-GM_Earth*tau1/12;
    H2 = H1-H3;
    
    G1 = -tau3/(tau1*(tau3-tau1));
    G3 = -tau1/(tau3*(tau3-tau1));
    G2 = G1-G3;
    
    D1 = G1+H1/magr1^3;
    D2 = G2+H2/magr2^3;
    D3 = G3+H3/magr3^3;
    
    temp = -[D1 D2 D3]*[C1 C2 C3]';
    rhoold1 = temp/(a1+b1*u);
    rhoold2 = -temp;
    rhoold3 = temp/(a3+b3*u);
    
    r1 = Rs1+rhoold1*Lm1;
    r2 = Rs2+rhoold2*Lm2;
    r3 = Rs3+rhoold3*Lm3;
    
end
r1 = Rs1+rho1*Lm1;
r2 = Rs2+rho2*Lm2;
r3 = Rs3+rho3*Lm3;