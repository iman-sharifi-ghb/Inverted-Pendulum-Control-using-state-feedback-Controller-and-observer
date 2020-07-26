clc;
clear
close all
%% System Equations
global M m l g
M = 5;m = 1;
l = 0.5;g = 9.81;
[A,B,C,D]=state_space();
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
animate_inv_pendulum(X,t1)
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
function animate_inv_pendulum(X,t)
    a = 2;
    b = 1;
    r = a/8;
    L = 3;
    R = 2*r;
    for i=1:length(t)
        xg = X(i,1);
        yg = 0;
        teta = X(i,3);
        %% cart
        x1 = xg-a/2:0.01:xg+a/2;
        y1 = (yg+b/2)*ones(1,length(x1));
        y2 = yg-b/2:0.01:yg+b/2;
        x2 = (xg+a/2)*ones(1,length(y2));
        x3 = (xg-a/2)*ones(1,length(y2));
        y3 = y2;
        x4 = x1;
        y4 =(yg-b/2)*ones(1,length(x1));
        %% wheels
        x5 = xg-a/4; y5=yg-b/2-r;
        x6 = xg+a/4; y6=yg-b/2-r;
        t1 = 0:0.01:2*pi;
        x55 = x5+r*cos(t1);y55 = y5+r*sin(t1);
        x66 = x6+r*cos(t1);y66 = y6+r*sin(t1);
        %% pendulum
        x7 = xg;y7=yg+b/2;
        x8 = x7+L*sin(teta);y8=y7+L*cos(teta);
        x9 = x8+R*cos(t1);y9 = y8+R*sin(t1);
        %% plot
        plot(x1,y1,'b-',x2,y2,'b-',x3,y3,'b-',x4,y4,'b-',x5,y5,'go',x6,y6,'go',...
            x55,y55,'g-',x66,y66,'g-',x7,y7,'bo',x8,y8,'ro',[x7 x8],[y7 y8],'g-',x9,y9,'r-','linewidth',3)
        axis([-2*a 2*a -1*a 2*a])
        grid on
        axis equal
%         text(a, 2*a,'Dr MEHRI KABIR')
        text(-2*a+0.1, 2*a,['Time= ' num2str(t(i)) ' (s)'])
        text(-2*a+0.1, 2*a-0.5,['X= ' num2str(X(i,1)) ' (m)'])
        text(-2*a+0.1, 2*a-1,['Theta= ' num2str(X(i,3)*180/3.14) ' (deg)'])
        pause (0.01)
%         hold off
    end
end
