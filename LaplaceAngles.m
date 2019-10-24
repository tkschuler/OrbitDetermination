function [r0, v0, oe0, rf, vf, oef] = LaplaceAngles(lat,lst,alt,rho,ra,dec,JD,JD_Prop)

%https://oaktrust.library.tamu.edu/bitstream/handle/1969.1/ETD-TAMU-2011-12-10242/SCHAEPERKOETTER-THESIS.pdf?sequence=2&isAllowed=y

mu=  398600.4354             % km^3/s^2  Earth's Gravitational Constant
RE = 6378.1366;              % km        Earth Radius
omega_E = 7.2921159e-5       % rad/s     Earth's intertial Rotation Rate


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


t= JD(1,1) + 15;
t1 = JD(1,1);
t2 = JD(2,1);
t3 = JD(3,1);

L1 = L(:,1)
L2 = L(:,2)
L3 = L(:,3)

L_t = ((t-t2)*(t-t3)/(t1-t2)*(t1-t3))*L1 + ((t-t1)*(t-t3)/(t2-t1)*(t2-t3))*L2 + ((t-t1)*(t-t2)/(t3-t1)*(t3-t2))*L3
L_dot_t = (2*t-t2-t3)/((t1-t2)*(t1-t3))*L1 + (2*t-t1-t3)/((t2-t1)*(t2-t3))*L2 + (2*t-t1-t2)/((t3-t1)*(t3-t2))*L3
L_dotdot_t = 2/((t1-t2)*(t1-t3))*L1 + 2/((t2-t1)*(t2-t3))*L2 + 2/((t3-t1)*(t3-t2))*L3

r_site1 = rsite_eci(:,1)
r_site2 = rsite_eci(:,2)
r_site3 = rsite_eci(:,3)
%how to calculate angular rate?
%Idk if this is right...

T1 = (JD(1,1)-JD(2,1))*24*60*60
T3 = (JD(3,1)-JD(2,1))*24*60*60

r_site2_dot = -T3/(T1*(T1-T3))*r_site1 - (T3+t1)/(T1*T3)*r_site2 -T1/(T3*(T3-T1))*r_site3
r_site2_dotdot = 2/(T1*(T1-T3))*r_site1 + 2/(T1*T3)*r_site2 + 2/(T3*(T3-T1))*r_site3

%Assume we know r?
D = 2*det([L_t 2*L_dot_t L_dotdot_t])