%% Arthur Binstein
% Eng 205
% Homework 7 Problem 3

clear all; clc; close all;
J = [20000 0 0; 0 40000 0; 0 0 60000];

s = tf('s')

Y = 1;
Kp = 1;
t = 0:.01:5;
u = [zeros(1,100) 100*ones(1,50) zeros(1,351)];

K = Kp*(1 + Y*s)

H1 = 1/(J(1,1)*s^2);

y1 = lsim(H1*K,u,t)
y1 = 180/pi*y1;
% plot(t,u);
hold on;
plot(t,y1)


% figure
% rlocus(H1*K)
% figure
% T = feedback(H1*K,1);
% step(T)