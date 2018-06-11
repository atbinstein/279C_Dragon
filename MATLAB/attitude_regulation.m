clear;
clc;
close all;

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

% nQb0 = [0 1 0; -1 0 0; 0 0 1]; %rotation to align +y body axis with +x ECI axis
% q0 = Qtoq(nQb0);
qd = [0 0 0 1]'; %Desired attitude
r = rand(3,1); r = r/norm(r);
theta_error = .1;
phi0 = theta_error*r;
q0 = phi2q(phi0); 
% q0 = qconj(qd);
% w0 = om0';
w0 = [.1 .1 .1]';

Kp = 10000*ones(3,1);
Kd = 100000*ones(3,1);
    

[r_eci, v_eci] = OE2ECI(a, e, inc, RAAN, w, anom, mu);

initCond = [r_eci; v_eci ; w0 ; q0];

orbitCount = 20;
stopTime = orbitCount/revs_per_day*24*60*60;
stepSize = 3;
tspan = [0:stepSize:stopTime];
% opts  = odeset('reltol', tol, 'abstol', tol);
rho = 0;

h = .1;
t_final = 500;
t = 0:h:t_final;
y = zeros(length(initCond),length(t));
phi = zeros(3,length(t));
y(:,1) = initCond;
phi(:,1) = q2phi(qmult(qconj(qd))*q0);
theta = zeros(1,length(t));
theta(1) = norm(phi(:,1))*180/pi;
u = zeros(3,length(t));
yn = zeros(size(y));
yn(:,1) = y(:,1);
xf = zeros(7,length(t));
xf(1:4,1) = y(10:13,1);
beta = zeros(3,length(t));
P = zeros(6,6,length(t));
P(:,:,1) = [(.05*pi/180)^2*eye(3) zeros(3);
       zeros(3) (.05*pi/180)^2*eye(3)];
wn = zeros(3,length(t));
wn(:,1) = w0;
for ii = 2:length(t)
    [yn(1:4,ii), ~] = addNoise(y(7:13,ii-1),Est,'StarTracker',rN,beta(:,ii-1));
    [yn(5:7,ii), ~] = addNoise(y(7:13,ii-1),Em,'Magno',rN,beta(:,ii-1));
    [yn(8:13,ii), ~] = addNoise(y(7:13,ii-1),Egps,'GPS',rN,beta(:,ii-1));
    [wn(:,ii), beta(:,ii)] = addNoise(y(7:13,ii-1),W,'Gyro',rN,beta(:,ii-1));
    
    [xf(:,ii), P(:,:,ii)] = mekfPID(xf(:,ii-1), P(:,:,ii-1), W, Vv, rN, wn(:,ii-1), yn(:,ii), h);
    
    
    u(:,ii) = Kp.*phi(:,ii-1) + Kd.*wn(:,ii-1);
    k1 = h*rk4(t(ii-1),y(:,ii-1),D,c,n,A,u(:,ii));
    k2 = h*rk4(t(ii-1)+h/2, y(:,ii-1)+k1/2,D,c,n,A,u(:,ii));
    k3 = h*rk4(t(ii-1)+h/2, y(:,ii-1)+k2/2,D,c,n,A,u(:,ii));
    k4 = h*rk4(t(ii-1)+h, y(:,ii-1)+k3,D,c,n,A,u(:,ii));
    y(:,ii) = y(:,ii-1) + (k1+2*k2+2*k3+k4)/6;
    phi(:,ii) = q2phi(qmult(qconj(qd))*xf(1:4,ii));
    theta(ii) = norm(phi(:,ii))*180/pi;
end
figure;
plot(t,theta,'LineWidth',2)
title('Error')
RMS = sum(theta)/length(theta)

y = y';
% figure
% hold on
% plot3(y(:,1),y(:,2),y(:,3),'c','LineWidth',1)
% earthPlot(1)
% axis equal
% hold off
% xlabel('x [km]')
% ylabel('y [km]')
% zlabel('z [km]')

figure
subplot(3,1,1)
plot(t,y(:,7),'b','LineWidth',1)
title('Angular Velocity Components')
subplot(3,1,2)
plot(t,y(:,8),'b','LineWidth',1)
subplot(3,1,3)
plot(t,y(:,9),'b','LineWidth',1)

for ii = 1:size(y,1)
    h(ii) = norm(I*y(ii,7:9)');
    ii;
end
figure;
plot(t,h,'b','LineWidth',2)
title('Angular Momentum')

% figure
% [X,Y,Z] = sphere(50);
% C(:,:,1) = .85.*ones(size(X,1));
% C(:,:,2) = .85.*ones(size(X,1));
% C(:,:,3) = .85.*ones(size(X,1));
% surf(X*Hr*.99, Y*Hr*.99, Z*Hr*.99, C, 'FaceAlpha',.5)
% lighting phong
% shading interp
% hold on;

% for ii = 1:size(y,1)
%     hr(ii,:) = (I*y(ii,10:12)')';
% end
% plot3(hr(:,1),hr(:,2),hr(:,3),'.','LineWidth',2)
% axis equal

figure
subplot(4,1,1);
plot(t,y(:,10),'b','LineWidth',1)
xlabel('t')
ylabel('q(1)')
title('q(1)')

subplot(4,1,2);
plot(t,y(:,11),'b','LineWidth',1)
xlabel('t')
ylabel('q(2)')
title('q(2)')

subplot(4,1,3);
plot(t,y(:,12),'b','LineWidth',1)
xlabel('t')
ylabel('q(3)')
title('q(3)')

subplot(4,1,4);
plot(t,y(:,13),'b','LineWidth',1)
xlabel('t')
ylabel('q(4)')
title('q(4)')

figure
subplot(3,1,1)
plot(t,phi(1,:))
title('Axis Angle')
subplot(3,1,2)
plot(t,phi(2,:))
subplot(3,1,3)
plot(t,phi(3,:))


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
drdt(7:9,1) = -J\(tauD + tauGGb + u - cross(omB, J*omB));
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