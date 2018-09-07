% This function generates the Spherical Harmonics basis functions of degree
% L and order M, given a specified asteroid SH model and output gravity potential.
%
% EROS SH model can be found at:
% https://sbn.psi.edu/pds/resource/nearbrowse.html
%
% SYNTAX: [V Gravity]=spharmAst(r,phi,lam,ref,GM);
%
% INPUTS:
%
% r   - distance from COM
% phi   - latitude
% lam   - longitude
% ref   - reference radius
% GM   - GM value of asteroid
%
% OUTPUTS:
% 
% Potential   - Gravity Potential
%
% NOTE: It is very important to keep the various usages of THETA and PHI
% straight.  For this function THETA is the Azimuthal/Longitude/Circumferential 
% coordinate and is defined on the interval [0,2*pi], whereas PHI is the 
% Altitude/Latitude/Elevation and is defined on the interval [-pi/2,pi/2].
%
% SUGIMOTO 2013/08/05

function Potential=spharmAst(r,phi,lam,ref,GM)

% if nargin==0 % Automatically take Eros case
%     r=16; % distance from the COM (km)
%     phi=0; % latitude (rad)
%     lam=0; % longitude (rad)
%     ref=16; % reference radius (km)
%     GM = 4.46275472004E-04; % GM value of Eros (km3/s2)
% end

EROS = csvread('SH_EROS_NORMAL.csv');
% EROS = csvread('SH_EROS.csv');
l = EROS(:,1);
m = EROS(:,2);
Clm = EROS(:,3);
Slm = EROS(:,4);

rho = r/ref;
Plm=[];

for i=1:15
    Plm = [Plm;legendre(i,sin(phi),'norm')];
end

V=GM/r; % C00-term's potential

for i = 1:numel(EROS(:,1))
    V=V+(GM/ref)/rho^(l(i)+1)*Plm(i)*(Clm(i)*cos(m(i)*lam)+Slm(i)*sin(m(i)*lam));
end
Potential=V;

return