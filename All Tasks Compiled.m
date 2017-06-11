% Daniela Attalla, Alexandra Tang, Pallav Rawat
close all;clear all; clc;

%% Task 1
L = 110*10^-6;                  % H
R = 1.45;                       % Ohm
Km = 43.5*10^-3;                % Nm/A
Kemf = 4.555*10^(-3)*60/(2*pi);   
J = 60*10^-7;                   %kg*m^2
Cv = 1.3*10^-6*60/(2*pi);       % Damping coefficient
n0 = 5500*2*pi/60;              % No-load speed

% State-space model
A = [-R/L -Kemf/L; Km/J 0];
B = [1/L; 0];
C = [0 1];
D = [0];

sys = ss(A,B,C,D);

% Transfer function, plant
[n,d] = ss2tf(A,B,C,D);
tf(n, d);

% Step response, plant
% figure()
% step(sys)

%% Task 3
Kp = 2;
Kd = Cv*n0/5;             % dynamic friction/5, we decided on 5 bc we wanted
                          % a faster response

%% Task 4
% figure()
% rlocus(sys)
% figure()
% pzmap(sys), sgrid

% TF, simplified
num = [Km/(R*J)];
den = [1 (Km*Kemf)/(R*J)];
sys2 = tf(num, den);

% figure()
% pzmap(sys2), sgrid

%% Task 5
Mh = 718*10^-3;         % Stall torque [Nm]
t2 = 50;                % Torque nutrunner (15-50 Nm)
n = sqrt(t2/Mh);        % Gear ratio

%% Task 6
% Gear inertia constants
rg = 0.0185;
rgr = 0.0163;
b = 0.006;
bc1 = 0.013;
bc2 = 0.0049;
Cgr = rg/rgr;

% System constants
dt = 1.5;
dj = dt;
kt = 739;

% Intertia
rog = 7850; %kg/m3

Jg1 = rog*pi*rgr^4*(9*b*n^2 + bc1*n^2 - 36*b*n + 52*b)/(32*(n-1)^4); 
Jg2 = rog*pi*rgr^4*(9*b*n^2 + bc2*n^2 - 36*b*n + 52*b)/(32*(n-1)^4); 
Jmg = J + Jg1 + Jg2/n;

l = (58+13.1+7+9.55)*10^-3;
dia = (10+11+14+19.6+27.7+9.55)/6*10^-3;
Jout = pi*rog*l*dia^4/32;

A= [0 0 1 0;...
    0 0 0 1; ...
    -kt/((n^4)*Jmg) kt/((n^2)*Jmg) 1/Jmg*(-Km*Kemf/R-dt/(n^4)) dt/(n^2*Jmg);...
    kt/((n^2)*Jout) -kt/Jout dt/((n^2)*Jout) -(dt+dj)/Jout];
B=[0 0 Km/(R*Jmg) 0]';
C=eye(4);
D=zeros(length(C),length(B(1,:)));

sys3 = ss(A,B,C,D);
%figure()
%step(sys3)

%% Task 7
n2 = n^2;

%% Task 8
Esteel = 200*10^9;
d = 8/1000;
p = 1.25/1000;
ds = d-(13*sqrt(3)/24)*p;
Acs = (pi/4)*(ds)^2;
%Acs = d^2*pi/4;
L = 0.1;
kj = Acs*Esteel/L; % translational
kr = kj*(p/(2*pi))^2;

A= [0 0 1 0;...
    0 0 0 1; ...
    -kt/((n^4)*Jmg) kt/((n^2)*Jmg) 1/Jmg*(-Km*Kemf/R-dt/(n^4)) dt/(n^2*Jmg);...
    kt/((n^2)*Jout) (-kt-kr)/Jout dt/((n^2)*Jout) -(dt+dj)/Jout];

%% Task 9
T = 25;     %target tightening torque [Nm]

%% Task 10 
mu = 0.2;   % friction coefficient

%% Task 11
ko = 2400;
do = 38;

mo = 3.9;   % mass of operator
mb =  0.85; % mass of battery
mt = 1.7;   % mass of tool

ro = 0.3;
rb = 0.498;
rt = 0.230;

J1 = mo*ro^2 + mb*rb^2 + mt*rt^2;

A = [0 0 0 1 0 0;...
    0 0 0 0 1 0; ...
    0 0 0 0 0 1; ...
   (kt-ko)/J1 kt/(J1*n^2) -kt/J1 (dt-do)/J1 dt/(J1*n^2) -dt/J1;...
   -kt/(Jmg*n^2) -kt/((n^4)*Jmg) kt/((n^2)*Jmg) -dt/(Jmg*n^2) 1/Jmg*(-Km*Kemf/R-dt/(n^4)) dt/(n^2*Jmg);...
   kt/Jout kt/((n^2)*Jout) (-kt-kr)/Jout dt/Jout dt/((n^2)*Jout) -(dt+dj)/Jout];
     
B = [0 0 0 0 Km/(R*Jmg) 0]';

C = eye(6);

D = [0 0 0 0 0 0]';
