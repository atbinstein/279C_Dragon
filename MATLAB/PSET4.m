%% Arthur Binstein
% AA279C SADC
% Problem Set 4
clear all; clc; close all;

%% Thruster Jacobian
u = 400; %Newtons
x = [1 0 0]';
y = [0 1 0]';
z = [0 0 1]';

rCOM = [-.11 2.298 0]';
r1 = [-1.224 3.25 0]' - rCOM;
r2 = r1; r3 = r1; r4 = r1;          %-x thruster pod
r5 = [0 3.25 1.224]' - rCOM; 
r6 = r5; r7 = r5; r8 = r5;          %+z thruster pod
r9 = [1.224 3.25 0]' - rCOM;
r10 = r9; r11 = r9; r12 = r9;       %+x thruster pod
r13 = [0 3.25 -1.224]' - rCOM; 
r14 = r13; r15 = r13; r16 = r13;    %-z thruster pod
r = [r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 r11 r12 r13 r15 r16];
%-x pod normal vectors
a1 = y;
a2 = y;
a3 = z;
a4 = -y;
%+z pod normal vectors
a5 = y;
a6 = y;
a7 = -x;
a8 = -y;
%+x pod normal vectors
a9 = y;
a10 = y;
a11 = -z;
a12 = -y;
%-z pod normal vectors
a13 = y;
a14 = y;
a15 = x;
a16 = -y;
a = [a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16];

Bt = [];
for ii = 1:size(r,2)
    Bt = [Bt cross(r(:,ii),a(:,ii))];
end

%% Gyrostat Jacobian
Bg = [x y z x+y+z];

%% LEO Drag

% Centroids of panels
% Top
c1 = [0 4.7 0]' - rCOM;
% Pressure Chamber
c2 = [1.46 2.94 0]' - rCOM; %+x
c3 = [0 2.94 -1.46]' - rCOM; %-z
c4 = [-1.46 2.94 0]' - rCOM; %-x
c5 = [0 2.94 1.46]' - rCOM; %+z
%Cargo Trunk
c6 = [1.65 1.5 0]' - rCOM; %+x
c7 = [0 1.5, -1.65]' - rCOM; %-z
c8 = [-1.65 1.5 0]' - rCOM; %-x
c9 = [0 1.5 1.65]' - rCOM; %+z
%Solar Panels
c10 = [5.85 1.5 0]';
c11 = c10;
c12 = [-5.85 1.5 0]';
c13 = c12;
c = [c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13];

% Normal vectors
% Top
n1 = y;
% Pressure Chamber
n2 = [.97 .26 0]'; %+x
n3 = [0 .26 -.97]'; %-z
n4 = [-.97 .26 0]'; %-x
n5 = [0 .26 .97]'; %+z
% Cargo Trunk
n6 = x;
n7 = -z;
n8 = -x;
n9 = z;
%Solar panels
n10 = z;
n11 = -z;
n12 = z;
n13 = -z;
n = [n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 n11 n12 n13];

% Areas
% Top
A1 = 3.94;
% Pressure Chamber
A2 = 9.42;
A3 = A2;
A4 = A2;
A5 = A2;
% Cargo Trunk
A6 = 11.1;
A7 = A6;
A8 = A6;
A9 = A6;
% Solar Panel
A10 = 24;
A11 = 24;
A12 = A10;
A13 = A10;
A = [A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13];

save('drag_props.mat','c','n','A')



