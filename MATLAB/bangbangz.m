%EignAxis Slew - Large Angle Maneuver
clear all; close all; clc;
load inertia.mat
D = diag([20000 20000 100000]);
%% Plan optimal trajectory and control
qd = [0 0 0 1]'; %Desired attitude
% q0 = [1 2 3 4]';
r = [0 0 1]'; r = r/norm(r);
theta_error = pi*.99;
phi0 = theta_error*r;
q0 = phi2q(phi0);
q0 = q0/norm(q0);
phi0 = q2phi(qmult(q0)*qd);
theta0 = norm(phi0);
tf = 24;
dt = .2;
t = 0:dt:tf;

% Now I have an optimal control plan/trajectory

%% Simulation
[a, e, inc, RAAN, w, M0, E, anom, revs_per_day] =...
    TLE_Reader('CRS-14_TLE.txt');

load inertia.mat
load drag_props.mat
load covariance.mat
J11 = D(1,1);
J22 = D(2,2);
J33 = D(3,3);
rN1 = [1 0 0]';
rN2 = [0 1 0]';
rN3 = [0 0 1]';
rN = [rN1 rN2 rN3];

mu = 398600;
J2 = 1.081874*10^-3;
Re = 6378137*10^-3;
tol = 1e-6;

q0 = qmult(q0)*phi2q([0 0 .01]');
% q0 = q0;
qd = [0 0 0 1]';
% q0 = qconj(qd);
% w0 = om0';
w0 = [0 0 0]';

[r_eci, v_eci] = OE2ECI(a, e, inc, RAAN, w, anom, mu);

initCond = [r_eci; v_eci ; w0 ; q0];

uMax = 1500;
uMin = -1500;
u = zeros(3,length(t));
for ii = 1:length(t)/2
    u(:,ii) = [0 0 uMin]';
end
for ii = floor(length(t)/2):length(t)
    u(:,ii) = [0 0 uMax]';
end


rho = 0;
h = dt;
y = zeros(length(initCond),length(t));
phi = zeros(3,length(t));
y(:,1) = initCond;
phi(:,1) = q2phi(qmult(qconj(qd))*q0);
theta = zeros(1,length(t));
theta(1) = norm(phi(:,1));
wn = zeros(3,length(t));
wn(:,1) = w0;
for ii = 2:length(t)    
    phi(:,ii) = q2phi(qmult(qconj(qd))*y(10:13,ii));
    theta(ii) = norm(phi(:,ii));
    k1 = h*rk4(t(ii-1),y(:,ii-1),D,c,n,A,u(:,ii));
    k2 = h*rk4(t(ii-1)+h/2, y(:,ii-1)+k1/2,D,c,n,A,u(:,ii));
    k3 = h*rk4(t(ii-1)+h/2, y(:,ii-1)+k2/2,D,c,n,A,u(:,ii));
    k4 = h*rk4(t(ii-1)+h, y(:,ii-1)+k3,D,c,n,A,u(:,ii));
    y(:,ii) = y(:,ii-1) + (k1+2*k2+2*k3+k4)/6;
end

q = y(10:13,:);
w = y(7:9,:);

vecplot(t,u);
vecplot(t,q);

save('BangBang.mat','q','u','t','w','D','dt','q0','qd');
function drdt = rk4(t,r,J,c,n,A,u)
mu = 398600;
rECI = r(1:3);
vECI = r(4:6);
omB = r(7:9);
qB = r(10:13);

tauGG = cross(3*mu/(dot(rECI,rECI))^(5/2)*rECI, q2Q(qB)*J/(1000^2)*rECI);
tauGG = tauGG*1000;
tauGGb = q2Q(qB)'*tauGG;
[D, tauD] = drag(c,n,A,r(4:6),rECI,qB);
Dn = q2Q(qB)*D/1000;
M = 6000;

drdt(1:6,1) = [vECI; -mu*rECI/norm(rECI)^3 - Dn/M];
drdt(7:9,1) = J\(tauD + tauGGb + u - cross(omB, J*omB));
drdt(10:13,1) = 1/2*qmult(qB)*[omB;0];
end

function [a, e, inc, RAAN, w, M, E, anom, revs_per_day] =...
    TLE_Reader(filename)
% OE2ECI Converts orbital elements to r, v in inertial frame
%
%   Notes: 
%       In cases of equatorial and/or circular orbits, it is assumed
%       that valid orbital elements are provided as inputs (ie. there is 
%       no back-end validation)
%
% Inputs:
%      filename - Two Line Element (TLE) .txt file
%
% Outputs:
%      a - semi-major axis of orbit [km]
%      e - eccentricity of orbit
%      inc - inclination of orbit [deg]
%      RAAN - right ascension of the ascending node [deg]
%      w - argument of periapsis [deg]
%      M - mean anomaly [deg]
%      E - eccentric anomaly [deg]
%      anom - true anomaly [deg]
%      revs_per_day - mean motion 

fid = fopen(filename, 'rb');
L1c = fscanf(fid,'%21c%',1);
L2c = fscanf(fid,'%71c%',1);
L3c = fscanf(fid,'%71c%',1);
fprintf(L1c);
fprintf(L2c);
fprintf([L3c,'\n']);
fclose(fid);

fid = fopen(filename, 'rb');
L1 = fscanf(fid,'%24c%*s',1);
L2 = fscanf(fid,'%d %6d %*c%5d%*3c%f%f%f%5d%*c%*d%5d%*c%*d%d%5d',[1,9]);
L3 = fscanf(fid,'%d%6d%f%f%f%f%f%f%f',[1,8]);
fclose(fid);

date  = L2(1,4);              % Epoch Date and Julian Date Fraction
Db    = L2(1,5);             % Ballistic Coefficient
inc   = L3(1,3);             % Inclination [deg]
RAAN  = L3(1,4);             % Right Ascension of the Ascending Node [deg]
e     = L3(1,5)/1e7;         % Eccentricity 
w     = L3(1,6);             % Argument of periapsis [deg]
M     = L3(1,7);             % Mean anomaly [deg]
n     = L3(1,8);             % Mean motion [Revs per day]

% Orbital elements
mu = 398600; %  Standard gravitational parameter for the earth
revs_per_day = n;
n = revs_per_day*2*pi/(24*3600);
a = (mu/n^2)^(1/3);     
E = M2E(M*pi/180,e,10^-6);  %[rad]
anom = E2anom(E, e);        %[rad]
E = E*180/pi;               %[deg]
anom = anom*180/pi;         %[deg]

% Six orbital elements 
OE = [a e inc RAAN w M];

%Date Conversion
if (date<57000)
    epoch = datestr(datenum(round(date/1000)+2000,1,0)+mod(date,1000));
else
    epoch = datestr(datenum(round(date/1000)+1900,1,0)+mod(date,1000));
end

fprintf('\n a[km]     e        inc[deg]   RAAN[deg]   w[deg]     M[deg]')
fprintf('\n %4.2f   %4.4f   %4.4f    %4.4f    %4.4f   %4.4f', OE);
fprintf('\n\nEpoch: %s UTC', epoch)
fprintf('\n\n---------- End of TLE Import ----------\n')

end

function [r_eci, v_eci] = OE2ECI(a, e, i, Om, w, anom, mu)
% OE2ECI Converts orbital elements to r, v in inertial frame
%
%   Notes: 
%       In cases of equatorial and/or circular orbits, it is assumed
%       that valid orbital elements are provided as inputs (ie. there is 
%       no back-end validation)
%
% Inputs:
%      a - semi-major axis of orbit [km]
%      e - eccentricity of orbit
%      i - inclination of orbit [deg]
%      Om - right ascension of the ascending node [deg]
%      w - argument of periapsis [deg]
%      anom - true anomaly [deg]
%      mu - central body gravitational parameters [km^3/s^2]
%
% Outputs:
%   r_eci - 3x1 vector of radius in ECI frame [km]
%   v_eci - 3x1 vector of velocity in ECI frame [km/s]

n = sqrt(mu/a^3);      % rad/s

E = anom2E(deg2rad(anom), e);    % rad

% Compute radius and velocity of orbit in perifocal coordinates
rPeri = [        a*(cos(E) - e);
         a*sqrt(1 - e^2)*sin(E);
                              0];

vPeriComp = [             -sin(E);
             sqrt(1 - e^2)*cos(E);
                                0];
vPeri = (a*n)/(1 - e*cos(E))*vPeriComp;

% Develop rotation matrix depending on orbit shape/inclination
if i == 0 && e ~= 0         % Equatorial + elliptical
    rotPeri2ECI = rotz(w);
elseif e == 0 && i ~= 0     % Circular + inclined
    rotPeri2ECI = rotz(Om)*rotx(i);
elseif i == 0 && e == 0     % Equatorial + circular
    rotPeri2ECI = 1;
else                        % Elliptical + inclined
    rotPeri2ECI = rotz(Om)*rotx(i)*rotz(w);
end
    
% Rotate vectors into ECI frame
r_eci = rotPeri2ECI*rPeri;
v_eci = rotPeri2ECI*vPeri;
end










