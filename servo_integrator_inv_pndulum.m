clc;
clear
close all
%% System Equations
global M m l g
M = 5;
m = 1;
l = 0.5;
g = 9.81;
[A,B,C,D]=state_space();
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
function [A,B,C,D]=state_space()
    syms x1 x2 x3 x4 u
    global m M l g
    dx1 = x2;
    dx2 = 1/(m+M-m*cos(x3))*(u-m*l*(x4)^2*sin(x3)-m*g*sin(x3));
    dx3 = x4;
    dx4 = (g*sin(x3)-(1/(m+M-m*cos(x3))*(u-m*l*(x4)^2*sin(x3)-m*g*sin(x3)))*cos(x3))/l;

    x = [x1;x2;x3;x4];
    dx = [dx1;dx2;dx3;dx4];

    A = jacobian(dx,x);
    A = simplify(A);
    B = jacobian(dx,u);
    B = simplify(B);

    A = subs(A,[x1,x2,x3,x4,u],[0,0,0,0,0]);
    B = subs(B,[x1,x2,x3,x4,u],[0,0,0,0,0]);

    A = vpa(A,6);
    B = vpa(B,6);
    A = double(A);
    B = double(B);
    C = [1 0 0 0];
    D = 0;
end


