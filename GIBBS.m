%--------------------------------------------------------------------------
%
%  GIBBS.m
%
%  this function uses Gibbs method to determine velocity of satelite
%  given 3 observations
% 
%  inputs:
%    r1       - 3x1 vector of ijk position of first observation             km
%    r2       - 3x1 vector of ijk position of second observation            km
%    r3       - 3x1 vector of ijk position of third observation             km
%
%  outputs:
%    v2       - 3x1 vector of ijk velocity of second observation            km/s
%
%
% Last modified:   10/22/2019   T. Schuler
% 
% -------------------------------------------------------------------------


function v2 = GIBBS(r1,r2,r3)

mu= 3.986004254*10^5; % Earth's Gravitational Constant
ER = 6378.137; % Earth Radius

z12 = cross(r1,r2);
z23 = cross(r2,r3);
z31 = cross(r3,r1);

N = norm(r1)*z23 + norm(r2)*z31 + norm(r3)*z12;
D = z12 + z23 + z31;
S = (norm(r2)-norm(r3))*r1 + (norm(r3)-norm(r1))*r2 + (norm(r1)-norm(r2))*r3;
B = cross(D,r2);

L_g = sqrt(mu/(norm(N)*norm(D))); 

v2 = (L_g/norm(r2))*B+L_g*S;

end