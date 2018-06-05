%% Arthur Binstein
% AA279C PSet 3
clear all;close all;clc;

%% Sensor Specs
%Star Tracker specs
sxa = 18; %arcsec
sya = 18; %arcsec
sza = 122; %arcsec

sx = sxa/2*(4.85*10^-5);
sy = sya/2*(4.85*10^-5);
sz = sza/2*(4.85*10^-5);

Est = zeros(3);
Est(1,1) = sx^2;
Est(2,2) = sy^2;
Est(3,3) = sz^2;
Est;
mu = [0 0 0]';

%Magno Specs
gx = .01*pi/180; %rad
gy = gx;
gz = gx;

Em = (gx^2)*eye(3);

rtrue0 = [1 1 0]';
rtrue0 = rtrue0/norm(rtrue0);
thtrue0 = 70*pi/180;
q0true = [rtrue0.*sin(thtrue0/2);cos(thtrue0/2)];
w0 = pi/180*[1 10 1]';

%Gyro Specs
Gyro_noise = .001*pi/180; %rad/sec/sqrt(Hz)
% Gyro_noise = .1*pi/180;
Gyro_WW = .001*pi/180;  %rad/sqrt(hr)
% Gyro_WW = .1*pi/180;
Egu = Gyro_noise^2*eye(3);
Egb = Gyro_WW^2*eye(3);

%GPS Specs
ggps = .01*pi/180;
Egps = ggps^2*eye(3);

Vv = [Est zeros(3) zeros(3) zeros(3);
     zeros(3) Em zeros(3) zeros(3);
     zeros(3) zeros(3) Egps zeros(3);
     zeros(3) zeros(3) zeros(3) Egps];
 
 W = [Egu zeros(3);
      zeros(3) Egb];
  
 save('covariance.mat','Est','Em','Egu','Egb','Egps','Vv','W')

%% Quaternion Propogation
x0 = [w0' q0true'];
tspan = [0:.1:100];
load inertia.mat
[t,y] = ode45(@(t,r) orb_prop(t,r,I, Egu, Egb), tspan, x0);
q_true = y(:,4:7)';
whist = y(:,1:3)';

for ii = 1:size(q_true,2)
     q_true(:,ii) = q_true(:,ii)/norm(q_true(:,ii));
end

for ii = 1:size(y,1)
%     y(ii,4:7);
    noisy_phi(:,ii) = chol(Est)*randn(3,1);
    noisy_phi(:,ii);
    qq = phi2q(noisy_phi(:,ii));
    qq = qq/norm(qq);
    noisy_quat(:,ii) = qmult(q_true(:,ii))*qq;
    norm(noisy_quat(:,ii));
end

% n_vecs = 3;
% rN = 2*rand(3,n_vecs)-1;
% for ii = 1:n_vecs
%    rN(:,ii) = rN(:,ii)/norm(rN(:,ii));
% end
% 
% rN1 = rN(:,1);
% rN2 = rN(:,2);
% rN3 = rN(:,3);
rN1 = [1 0 0]';
rN2 = [0 1 0]';
rN3 = [0 0 1]';

rB10 = q2Q(q0true)'*rN1;
rB20 = q2Q(q0true)'*rN2;
rB30 = q2Q(q0true)'*rN3;

rB1hist = rB10;
rB2hist = rB20;
rB3hist = rB30;
rB1 = rB10;
rB2 = rB20;
rB3 = rB30;
for ii = 2:size(y,1)
    rB1 = q2Q(q_true(:,ii))'*rN1;
    rB2 = q2Q(q_true(:,ii))'*rN2;
    rB3 = q2Q(q_true(:,ii))'*rN3;
    rB1hist = [rB1hist rB1];
    rB2hist = [rB2hist rB2];
    rB3hist = [rB3hist rB3];
end

rB1noisy =[];
rB2noisy =[];
rB3noisy =[];
for ii = 1:size(y,1)
    rB1noisy = [rB1noisy rB1hist(:,ii) + chol(Em)*randn(3,1)];
    rB1noisy(:,ii) = rB1noisy(:,ii)/norm(rB1noisy(:,ii));
    rB2noisy = [rB2noisy rB2hist(:,ii) + chol(Egps)*randn(3,1)];
    rB2noisy(:,ii) = rB2noisy(:,ii)/norm(rB2noisy(:,ii));
    rB3noisy = [rB3noisy rB3hist(:,ii) + chol(Egps)*randn(3,1)];
    rB3noisy(:,ii) = rB3noisy(:,ii)/norm(rB3noisy(:,ii));
end

btrue = zeros(size(whist));
for ii = 1:size(y,1)-1
    btrue(:,ii+1) = btrue(:,ii) + chol(Egb)*randn(3,1);
end
wnoisy = zeros(size(whist));
for ii = 1:size(y,1)
    wnoisy(:,ii) = whist(:,ii) + btrue(:,ii) + chol(Egu)*randn(3,1);
end




%% Whaba's Problem
% SVD
B = zeros(3,3);
error = [];
qsvdhist = [];
trQ = [];
for ii = 1:size(y,1)
    B = zeros(3,3);
    B = 1/gx*rB1noisy(:,ii)*rN1' + 1/ggps*rB2noisy(:,ii)*rN2' ...
        + 1/ggps*rB3noisy(:,ii)*rN3';

    [U E V] = svd(B');
    E = eye(3);
    det(U);
    det(V);
    E(3,3) = det(U)*det(V);
    Qest = (U*E*V');
    qest = Qtoq(Qest);
    for jj = 1:4
        qest(jj) = sign(q_true(jj,ii))*abs(qest(jj));
    end
    qsvdhist  = [qsvdhist qest];
    qe = qmult(qconj(q_true(:,ii)))*qest;
    trQ = [trQ trace(Qest)];
    e = q2phi(qe);
    error = [error 180/pi*norm(e)];
end

%Davenport q-method
error2 = [];
qdavest = [];
for ii = 1:size(y,1)
    z = zeros(3,1);
    B = zeros(3,3);
    z = 1/gx*cross(rB1noisy(:,ii),rN1) + 1/ggps*cross(rB2noisy(:,ii),rN2) ...
        + 1/ggps*cross(rB3noisy(:,ii),rN3);
    
    B = 1/gx*rB1noisy(:,ii)*rN1' + 1/ggps*rB2noisy(:,ii)*rN2' ...
        + 1/ggps*rB3noisy(:,ii)*rN3';
    K = [B + B' - trace(B)*eye(3) z ; z' trace(B)];

    [V,d] = eig(K);
    for jj = 1:4
        V(jj,end) = sign(q_true(jj,ii))*abs(V(jj,end));
    end
    
    q = V(:,end)/norm(V(:,end));
    qdavest = [qdavest q];
    qe = qmult(qconj(q))*q_true(:,ii);
    phie = q2phi(qe);
    error2 = [error2 180/pi*norm(phie)];
end

% figure;
% plot(t,qdavest(1,:),t,q_true(1,:))
% title('Davenport Quaternion vs True Quaternion');
% figure;
% plot(t,q_true(1,:),t,noisy_quat(1,:))
% title('true vs noisy')

% %% Plots
% figure;
% plot(t,y(:,4),t,noisy_quat(1,:));
% figure;
% plot(t,y(:,5),t,noisy_quat(2,:));
% figure;
% plot(t,y(:,6),t,noisy_quat(3,:));
% figure;
% plot(t,y(:,7),t,noisy_quat(4,:));
% 
figure;
plot(t,qsvdhist(1,:),t,q_true(1,:))
title('SVD Quaternion vs. True Quaternion')
xlabel('time (s)')
ylabel('error (degrees)')
legend('Estimated','True')

figure;
plot(t,qdavest(1,:),t,q_true(1,:))
title('Davenport Quaternion vs True Quaternion');
% 
figure;
plot(t,q_true(1,:),t,noisy_quat(1,:))
title('true vs noisy quaternion')
% 
% figure;
% plot(t, rB3hist(1,:),t, rB3noisy(1,:))

figure;
plot(t,error)
title('SVD Quaternion Error')
xlabel('time (s)')
ylabel('error (degrees)')

% figure;
% plot(t,error2)
% figure;
% plot(t,btrue(1,:))
% figure;
% plot(t,rB1hist(1,:)-rB1noisy(1,:))
% load mekf_truth
dt = .1;
save('mekf_inputs_ps3.mat','noisy_quat','rB1noisy','rB2noisy','rB3noisy', ...
    'rN1','rN2','rN3','dt','wnoisy','Vv','W','q0true')
save('mekf_truth_ps3.mat','q_true','btrue')

%% Function Calls
function dqdt = orb_prop(t,r,I, Egu, Egb)
r(4:7,1) = r(4:7,1)/norm(r(4:7,1));
% norm(r(4:7,1));
dqdt(1:3,1) = -I\cross(r(1:3), I*r(1:3));
dqdt(4:7,1) = 1/2*qmult(r(4:7))*[r(1:3);0];
% dqdt(8:10,1) = sqrt(Egb)*randn(3,1)*sqrt(dt);
end

