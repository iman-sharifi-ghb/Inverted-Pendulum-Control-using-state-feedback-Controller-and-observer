clc;
clear
close all
%% System Equations
global M m l g
M = 5;m = 1;
l = 0.5;g = 9.81;
[A,B,C,D]=State_Space();
%DESIGN STATE_FEEDBACK CONTROLLER
rank(ctrb(A,B))
desired_poles = [-2+1j -2-1j -5 -5];
K = acker(A,B,desired_poles);
%% LINEAR ODE45
% init = [0.5 1 -10*3.14/180 -5*3.14/180];
init = [1 -0.5 20*3.14/180 -10*3.14/180];
options = odeset('RelTol',1e-2,'AbsTol',1e-4);
tspan = 0:0.01:10;
[t,X] = ode45(@(t,x) linear_ode(t,x,A,B,K),tspan,init,options);

[t,XX] = ode45(@(t,x) nonlinear_ode(t,x,K),tspan,init,options);

%% ANIMATION
t1 = 0:0.01:5;
Animate_Inv_Pdlm(X,t1)
%% plots
% subplot(2,2,1);plot(t,X(:,1));title('X');
% subplot(2,2,2);plot(t,X(:,2));title('X-dot');
% subplot(2,2,3);plot(t,X(:,3)/3.14*180);title('Teta');
% subplot(2,2,4);plot(t,X(:,4)/3.14*180);title('Teta-dot');
% figure;
% subplot(2,2,1);plot(t,XX(:,1));title('X');
% subplot(2,2,2);plot(t,XX(:,2));title('X-dot');
% subplot(2,2,3);plot(t,XX(:,3)/3.14*180);title('Teta');
% subplot(2,2,4);plot(t,XX(:,4));title('Teta-dot');
% figure;
% plot(t,X(:,1),t,XX(:,1),'r')
% legend('Linear Equation','NonLinear Equation')
% title('Compare X')
% xlabel('Time');ylabel('X')
% figure;
% plot(t,X(:,3),t,XX(:,3),'r')
% legend('Linear Equation','NonLinear Equation')
% title('Compare Teta')
% xlabel('Time');ylabel('Teta')
%%
function dx = nonlinear_ode(t,x,K)
    global M m l g
    u = -K*x;
    dx = [x(2);
        1/(m+M-m*cos(x(3)))*(u-m*l*(x(4))^2*sin(x(3))-m*g*sin(x(3)));
        x(4);
        (g*sin(x(3))-(1/(m+M-m*cos(x(3)))*(u-m*l*(x(4))^2*sin(x(3))-m*g*sin(x(3))))*cos(x(3)))/l];
end
function dx = linear_ode(t,x,A,B,K)
    u = -K*x;
    dx = A*x + B*u;
end
