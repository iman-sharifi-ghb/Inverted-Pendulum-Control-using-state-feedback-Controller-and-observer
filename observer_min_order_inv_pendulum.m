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
%% CONTROLLERR
dt = 0.01;
T = 10;
t = 0:dt:T;
u = sign(t);
rank(ctrb(A,B))
desired_poles = [-2+1j -2-1j -5 -5];
K = acker(A,B,desired_poles);
disp(K)
Acl = A-B*K;
%% DESIGN STATE_OBSERVER
rank(obsv(A,C))
Desired_poles = [-10 -10 -15 -15];
L = acker(A',C',Desired_poles)';
disp(L)
%% DESIGN REDUCE-ORDER OBSERVER
Aaa = A(1,1);Aab = A(1,2:4);
Aba = A(2:4,1);Abb = A(2:4,2:4);
Ba = B(1);Bb = B(2:4);
rank(obsv(Abb,Aab))
Des_poles = [-5 -5 -5];
Lr = acker(Abb',Aab',Des_poles)';
%% SOLVE EQUATIONS
in = input("Inter 0 or 1:");
init = [0 0 -5*3.14/180 -2*3.14/180 0 0 0 0];
options = odeset('RelTol',1e-2,'AbsTol',1e-4);
tspan = 0:0.01:10;
if in==0
    [t,XX] = ode45(@(t,x) linear_ode_order_reduction(t,x,A,B,K,Lr),tspan,init,options);
else 
    [t,XX] = ode45(@(t,x) nonlinear_ode_reduce_order(t,x,K,Lr),tspan,init,options);  
end

%% PLOTS
X = XX(:,1:4);
subplot(2,2,1);plot(t,X(:,1));title('X');
subplot(2,2,2);plot(t,X(:,2));title('X-dot');
subplot(2,2,3);plot(t,X(:,3)/3.14*180);title('Teta');
subplot(2,2,4);plot(t,X(:,4)/3.14*180);title('Teta-dot');
figure
Xhat = XX(:,5:8);
subplot(2,2,1);plot(t,Xhat(:,1));title('X');
subplot(2,2,2);plot(t,Xhat(:,2));title('X-dot');
subplot(2,2,3);plot(t,Xhat(:,3)/3.14*180);title('Teta');
subplot(2,2,4);plot(t,Xhat(:,4)/3.14*180);title('Teta-dot');
figure
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
%% FUNCTIONS
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
function dx = linear_ode_order_reduction(t,XX,A,B,K,Lr)
    Aaa = A(1,1);Aab = A(1,2:4);
    Aba = A(2:4,1);Abb = A(2:4,2:4);
    Ba = B(1);Bb = B(2:4);
    
    Xa = XX(1);
    Xb = XX(2:4);
    Xhat_b = XX(6:8);
    u = -K*[Xa;Xhat_b];
    
    dXa = Aaa*Xa + Aab*Xb + Ba*u;
    dXb = Aba*Xa + Abb*Xb + Bb*u;
    dXhat_a = dXa;
    dXhat_b = Aba*Xa + Abb*Xhat_b + Bb*u + Lr*(Aab*Xb - Aab*Xhat_b);

    dx=[dXa;dXb;dXhat_a;dXhat_b];
end
function dx = nonlinear_ode_reduce_order(t,XX,K,Lr)
    global M m l g
    Xa = XX(1);
    Xb = XX(2:4);
    Xhat_b=XX(6:8);
    u = -K*[Xa;Xhat_b];
    
    dXa = Xb(1);
    dXb =[1/(m+M-m*cos(Xb(2)))*(u-m*l*(Xb(3))^2*sin(Xb(2))-m*g*sin(Xb(2)));
         Xb(3);
        (g*sin(Xb(2))-(1/(m+M-m*cos(Xb(2)))*(u-m*l*(Xb(3))^2*sin(Xb(2))-m*g*sin(Xb(2))))*cos(Xb(2)))/l];
%     dXhat = A*Xhat + B*u + L*C*(X-Xhat);
    dXhat_a = dXa;
    dXhat_b=[
        1/(m+M-m*cos(Xhat_b(2)))*(u-m*l*(Xhat_b(3))^2*sin(Xhat_b(2))-m*g*sin(Xhat_b(2)));
        Xhat_b(3);
        (g*sin(Xhat_b(2))-(1/(m+M-m*cos(Xhat_b(2)))*(u-m*l*(Xhat_b(3))^2*sin(Xhat_b(2))-m*g*sin(Xhat_b(2))))*cos(Xhat_b(2)))/l];
    dXhat_b = dXhat_b + Lr*(Xb(1)-Xhat_b(1));
    dx = [dXa;dXb;dXhat_a;dXhat_b];
    
end
 


