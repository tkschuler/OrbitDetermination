%--------------------------------------------------------------------------
%
%  HGIBBS.m
%
%  this function uses Herricks-Gibbs method to determine velocity of satellite
%  given 3 observations that are "close"
% 
%  inputs:
%    r1       - 3x1 vector of ijk position of first observation             km
%    r2       - 3x1 vector of ijk position of second observation            km
%    r3       - 3x1 vector of ijk position of third observation             km
%    JD       - 3x1 vector of JDs of observations
%
%  outputs:
%    v2       - 3x1 vector of ijk velocity of second observation            km/s
%
%
% Last modified:   10/24/19 J.Hammes
% 
% -------------------------------------------------------------------------


function v2 = HGIBBS(r1,r2,r3,JD)

mu= 3.986004254*10^5; % Earth's Gravitational Constant

t31 = (JD(3,1)-JD(1,1))*24*60*60;
t32 = (JD(3,1)-JD(2,1))*24*60*60;
t21 = (JD(2,1)-JD(1,1))*24*60*60;
 
v2 = -t32*(1/(t21*t31)+mu/(12*norm(r1)^3))*r1 + (t32-t21)*(1/(t21*t32)+mu/(12*norm(r2)^3))*r2 + t21*(1/(t32*t31)+mu/(12*norm(r3)^3))*r3;
end