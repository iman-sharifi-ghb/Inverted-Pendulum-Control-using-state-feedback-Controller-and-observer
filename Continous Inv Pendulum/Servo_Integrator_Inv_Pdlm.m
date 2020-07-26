clc;
clear
close all
%% System Equations
global M m l g
M = 5;
m = 1;
l = 0.5;
g = 9.81;
[A,B,C,D]=State_Space();
%%
dt = 0.01;
T = 50;
Abar = [A zeros(4,1);-C 0];
Bbar = [B;0];
Bb= [zeros(4,1);1];
CC = [C 0];
% DESIGN STATE FEEDBACK CONTROLLER
desired_poles = [-2+1j -2-1j -5 -5 -8];
K = acker(Abar,Bbar,desired_poles);
Acl = Abar-Bbar*K;
%% Plots
in = input("Inter 0 or 1:");
init = [0 0 5*3.14/180 -2*3.14/180 0];
options = odeset('RelTol',1e-2,'AbsTol',1e-4);
tspan = 0:dt:T;
if in==0
    [t,X] = ode45(@(t,x) linear_ode(t,x,Abar,Bbar,Bb,K),tspan,init,options);
else 
    [t,X] = ode45(@(t,x) nonlinear_ode(t,x,K),tspan,init,options);  
end
Yr = 2*sign(sin(0.5*t));
plot(t,Yr,t,X(:,1))
title('X')
xlabel('Time');ylabel('X')

function dX = linear_ode(t,X,A,B,BB,K)
    Yr = 2*sign(sin(0.5*t));
    u = -K*X;
    dX = A*X + B*u + BB*Yr;
end
function dX = nonlinear_ode(t,X,K)
    global M m l g
    Yr = 2*sign(sin(0.5*t));
    u = -K*X;
    dX = [X(2);
        1/(m+M-m*cos(X(3)))*(u-m*l*(X(4))^2*sin(X(3))-m*g*sin(X(3)));
        X(4);
        (g*sin(X(3))-(1/(m+M-m*cos(X(3)))*(u-m*l*(X(4))^2*sin(X(3))-m*g*sin(X(3))))*cos(X(3)))/l;
        Yr - X(1)];
end
