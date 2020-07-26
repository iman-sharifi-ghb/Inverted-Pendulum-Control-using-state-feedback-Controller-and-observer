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
T = 10;
t = 0:dt:T;
u = sign(t);
rank(ctrb(A,B))
desired_poles = [-2+1j -2-1j -5 -5];
K = acker(A,B,desired_poles);
disp(K)
Acl = A-B*K;
% DESIGN STATE_OBSERVER
rank(obsv(A,C))
Desired_poles = [-10 -10 -15 -15];
L = acker(A',C',Desired_poles)';
disp(L)
%%
% init=[0 0 0 0];
% Ahat = A - L*C;
% Bhat = [B L];
% Chat = C;
% Dhat = D;
% sys2 = ss(Ahat,Bhat,Chat,Dhat);
% u=0;
% for i=1:length(t)
%     [Yhat(i,1), ~, Xhat(i,1)]=lsim(sys2, [u(i,1)' Y(i,1)], t(i) ,init);
%     u(i,1)=-K*Xhat(i,:);
% end
%% Plots
in = input("Inter 0 or 1:");
init = [0 0 -5*3.14/180 -2*3.14/180 0 0 0 0];
% init = [0 0 20*3.14/180 -10*3.14/180 0 0 0 0];
options = odeset('RelTol',1e-2,'AbsTol',1e-4);
tspan = 0:0.01:10;
if in==0
    [t,XX] = ode45(@(t,x) linear_ode(t,x,A,B,C,K,L),tspan,init,options);
else 
    [t,XX] = ode45(@(t,x) nonlinear_ode(t,x,C,K,L),tspan,init,options);  
end
X = XX(:,1:4);
% subplot(2,2,1);plot(t,X(:,1));title('X');
% subplot(2,2,2);plot(t,X(:,2));title('X-dot');
% subplot(2,2,3);plot(t,X(:,3)/3.14*180);title('Teta');
% subplot(2,2,4);plot(t,X(:,4)/3.14*180);title('Teta-dot');
% figure
Xhat = XX(:,5:8);
% subplot(2,2,1);plot(t,Xhat(:,1));title('X');
% subplot(2,2,2);plot(t,Xhat(:,2));title('X-dot');
% subplot(2,2,3);plot(t,Xhat(:,3)/3.14*180);title('Teta');
% subplot(2,2,4);plot(t,Xhat(:,4)/3.14*180);title('Teta-dot');
% figure
plot(t,Xhat(:,1));figure
plot(t,X(:,1),t,Xhat(:,1))
legend('state','observer')
title('X')
xlabel('Time');ylabel('X')
figure
plot(t,X(:,2),t,Xhat(:,2))
legend('state','observer')
title('X-dot')
xlabel('Time');ylabel('X-dot')
figure
plot(t,X(:,3)*180/3.14,t,Xhat(:,3)*180/3.14)
legend('state','observer')
title('Teta')
xlabel('Time');ylabel('Teta')
figure
plot(t,X(:,4)*180/3.14,t,Xhat(:,4)*180/3.14)
legend('state','observer')
title('Teta-dot')
xlabel('Time');ylabel('Teta-dot')

function dx = linear_ode(t,XX,A,B,C,K,L)
    X = XX(1:4);
    Xhat=XX(5:8);
    u = -K*Xhat;
    dX = A*X + B*u;
    Y = C*X;
    Yhat = C*Xhat;
    dXhat = A*Xhat + B*u + L*(Y-Yhat);
    dx=[dX;dXhat];
end
function dx = nonlinear_ode(t,XX,C,K,L)
    global M m l g
    X = XX(1:4);
    Xhat=XX(5:8);
    u = -K*Xhat;
    Y = C*X;
    Yhat = C*Xhat;
    dX = [X(2);
        1/(m+M-m*cos(X(3)))*(u-m*l*(X(4))^2*sin(X(3))-m*g*sin(X(3)));
        X(4);
        (g*sin(X(3))-(1/(m+M-m*cos(X(3)))*(u-m*l*(X(4))^2*sin(X(3))-m*g*sin(X(3))))*cos(X(3)))/l];
%     dXhat = A*Xhat + B*u + L*C*(X-Xhat);
    dXhat = [Xhat(2);
        1/(m+M-m*cos(Xhat(3)))*(u-m*l*(Xhat(4))^2*sin(Xhat(3))-m*g*sin(Xhat(3)));
        Xhat(4);
        (g*sin(Xhat(3))-(1/(m+M-m*cos(Xhat(3)))*(u-m*l*(Xhat(4))^2*sin(Xhat(3))-m*g*sin(Xhat(3))))*cos(Xhat(3)))/l];
    dXhat = dXhat + L*(Y-Yhat);
    dx = [dX;dXhat];
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
 


