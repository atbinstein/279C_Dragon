clear all; close all; clc;

load mekf_inputs_ps3
dt = .1;
whist = wnoisy;
rN = [rN1, rN2, rN3];
yhist = [noisy_quat; rB1noisy; rB2noisy; rB3noisy];
% yhist = [q0true; rB1hist; rB2hist; rB3hist];
q0 = [0 0 0 1]';
V = Vv;

x0 = [q0; 0; 0; 0]; %Initialize with no bias
P0 = [(.05*pi/180)^2*eye(3) zeros(3);
       zeros(3) (.05*pi/180)^2*eye(3)]; %10 deg. and 10 deg/sec 1-sigma uncertainty

[xhist,Phist] = mekf(x0,P0,W,V,rN,whist,yhist,dt);

load mekf_truth_ps3

%Calculate error quaternions
e = zeros(3,size(q_true,2));
for k = 1:size(q_true,2)
    q_true(:,k);
    xhist(1:4,k);
    e(:,k) = q2phi(qmult(qconj(q_true(:,k)))*xhist(1:4,k));
end

% e = zeros(size(qtrue));
% for k = 1:size(qtrue,2)
%     e(:,k) = qmult(qconj(qtrue(:,k)), xhist(1:4,k));
% end
%     Q = q2Q(q);
%     E = q2Q(q_true(:,ii))'*Q;
%     e = unhat(logm(E));
%     error2 = [error2 180/pi*norm(e)];


%----- Plots -----%
figure(2);
subplot(4,1,1);
plot(q_true(1,:));
hold on;
plot(xhist(1,:));
title('Attitude');
legend('True', 'Estimated');
subplot(4,1,2);
plot(q_true(2,:));
hold on;
plot(xhist(2,:));
subplot(4,1,3);
plot(q_true(3,:));
hold on;
plot(xhist(3,:));
subplot(4,1,4);
plot(q_true(4,:));
hold on;
plot(xhist(4,:));

figure(3);
subplot(3,1,1);
plot((360/pi)*e(1,:));
hold on
plot((360/pi)*sqrt(squeeze(Phist(1,1,:))),'r');
plot(-(360/pi)*sqrt(squeeze(Phist(1,1,:))),'r');
title('Attitude Error');
subplot(3,1,2);
plot((360/pi)*e(2,:));
hold on
plot((360/pi)*sqrt(squeeze(Phist(2,2,:))),'r');
plot(-(360/pi)*sqrt(squeeze(Phist(2,2,:))),'r');
ylabel('degrees');
subplot(3,1,3);
plot((360/pi)*e(3,:));
hold on
plot((360/pi)*sqrt(squeeze(Phist(3,3,:))),'r');
plot(-(360/pi)*sqrt(squeeze(Phist(3,3,:))),'r');
% 
figure(4);
subplot(3,1,1);
plot(xhist(5,:)-btrue(1,:));
hold on
plot(2*sqrt(squeeze(Phist(4,4,:))),'r');
plot(-2*sqrt(squeeze(Phist(4,4,:))),'r');
title('Bias Error');
subplot(3,1,2);
plot(xhist(6,:)-btrue(2,:));
hold on
plot(2*sqrt(squeeze(Phist(5,5,:))),'r');
plot(-2*sqrt(squeeze(Phist(5,5,:))),'r');
subplot(3,1,3);
plot(xhist(7,:)-btrue(3,:));
hold on
plot(2*sqrt(squeeze(Phist(6,6,:))),'r');
plot(-2*sqrt(squeeze(Phist(6,6,:))),'r');
figure(5);
hold on
plot(btrue(1,:))
plot(xhist(5,:))

