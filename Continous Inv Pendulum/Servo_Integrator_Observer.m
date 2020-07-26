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
%% DESIGN STATE_FEEDBACK CONTROLLER
dt = 0.01;
T = 30;
Abar = [A zeros(4,1);-C 0];
Bbar = [B;0];
BB= [zeros(4,1);1];
CC = [C 0];
%% DESIGN STATE FEEDBACK CONTROLLER
desired_poles = [-2+1j -2-1j -2 -2 -4];
K = acker(Abar,Bbar,desired_poles);
Acl = Abar-Bbar*K;
%% DESIGN STATE_OBSERVER
rank(obsv(A,C))
Desired_poles = [-5 -5 -5 -5];
L = acker(A',C',Desired_poles)';

%% SOLVE EQUATIONS
in = input("Inter 0 or 1:");
init = [0 0 5*3.14/180 -2*3.14/180 0 0 0 0 0];
options = odeset('RelTol',1e-2,'AbsTol',1e-4);
tspan = 0:dt:T;
if in==0
    [t,X] = ode45(@(t,x) linear_ode(t,x,A,B,Abar,Bbar,BB,K,L),tspan,init,options);
else 
    [t,X] = ode45(@(t,x) nonlinear_ode(t,x,K,L),tspan,init,options);  
end
%% PLOTs
plot(t,X(:,1),'b',t,X(:,5)/2.03,'r',t,0.3*sign(sin(0.5*t)),'g')
title('X')
xlabel('Time');ylabel('X')
%% FUNCTIONS
function dx = nonlinear_ode(t,XX,K,L)
    Yr = 0.3*sign(sin(0.5*t));
    X = XX(1:5);
    Xhat=XX(6:9);
    u = -K*X;
    
    global M m l g
    dX = [  X(2);
            1/(m+M-m*cos(X(3)))*(u-m*l*(X(4))^2*sin(X(3))-m*g*sin(X(3)));
            X(4);
            (g*sin(X(3))-(1/(m+M-m*cos(X(3)))*(u-m*l*(X(4))^2*sin(X(3))-m*g*sin(X(3))))*cos(X(3)))/l;
            Yr-X(1)];

    dXhat = [Xhat(2);
            1/(m+M-m*cos(Xhat(3)))*(u-m*l*(Xhat(4))^2*sin(Xhat(3))-m*g*sin(Xhat(3)));
            Xhat(4);
            (g*sin(Xhat(3))-(1/(m+M-m*cos(Xhat(3)))*(u-m*l*(Xhat(4))^2*sin(Xhat(3))-m*g*sin(Xhat(3))))*cos(Xhat(3)))/l];
    dXhat = dXhat + L*(Yr-Xhat(1));
    dx = [dX;dXhat];
end
function dx = linear_ode(t,XX,A,B,Abar,Bbar,BB,K,L)
    Yr = 0.3*sign(sin(0.5*t));
    
    X = XX(1:5);
    Xhat=XX(6:9);
    u = -K*X;
    dX = Abar*X + Bbar*u + [zeros(4,1);Yr];
    Yhat = Xhat(1);
    dXhat = A*Xhat + B*u + L*(Yr-Yhat);
    dx=[dX;dXhat];
end
