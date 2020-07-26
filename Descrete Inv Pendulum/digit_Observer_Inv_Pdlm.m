clc;
clear
close all
%% System Equations
global M m l g
M = 3;m = 1;
l = 0.5;g = 9.81;
[A,B,C,D]=inv_pendulum_ss();
rank(ctrb(A,B))
h = 0.2;
% while h<1
%     G = expm(A*h);
%     la=eig(G);
%     if la(1)<1 && la(2)<1 && la(3)<1 && la(4)<1
%         break
%     end
%     h=h+0.01;
% end
G = expm(A*h);
%syms tav
% H = int(G,tav,[0 h])*B; H = vpa(H,4);
H=(h*eye(4)+1/2*A*h^2+1/6*A^2*h^3+1/24*A^3*h^4+1/120*A^4*h^5)*B;
rank(ctrb(G,H))

des_poles = [0 0 0 0];
K=acker(G,H,des_poles);

rank(obsv(G,C))
des_poles =[0.1 0.1 -0.1 -0.1];
L=acker(G',C',des_poles)';

T = 10;
dt = 0.01;
t=0;
k = 1;k1= 0;
N = floor(h/dt);
X(:,k) = 0.1*[0;0;2*3.14/180;-1*3.14/180];
Xh(:,k) = [0;0;0;0];
% %% LINEAR SYSTEM
% while t<T
%     if mod(k,N)==1
%        k1 = k1+1;
%        u(k)=-K*Xh(:,k1);
%        Y=C*X(:,k);
%        Yh=C*Xh(:,k1);
%        Xh(:,k1+1)=G*Xh(:,k1)+H*u(k)+L*(Y-Yh);  
%     else
%        u(k)=u(k-1);
%     end
%     X(:,k+1)=X(:,k)+(A*X(:,k)+B*u(k))*dt;
%     k = k+1;
%     t=t+dt;    
% end
%% NONLINEAR SYSTEM
while t<T
    if mod(k,N)==1
       k1 = k1+1;
       u(k)=-K*Xh(:,k1);
       Y=C*X(:,k);
       Yh=C*Xh(:,k1);
       Xh(:,k1+1)=G*Xh(:,k1)+H*u(k)+L*(Y-Yh);  
    else
       u(k)=u(k-1);
    end
    SYS =[X(2,k);
          1/(m+M-m*cos(X(3,k)))*(u(k)-m*l*X(4,k)^2*sin(X(3,k))-m*g*sin(X(3,k)));
          X(4,k);
         (g*sin(X(3,k))-(1/(m+M-m*cos(X(3,k)))*(u(k)-m*l*(X(4,k))^2*sin(X(3,k))-m*g*sin(X(3,k))))*cos(X(3,k)))/l];

    X(:,k+1)=X(:,k)+(SYS)*dt;
    k = k+1;
    t=t+dt;    
end

%% plots
Time=0:dt:T+dt;
Time1=0:h:T+h;
plot(Time,X(1,:),Time1,Xh(3,:));title('X');
xlabel('Time');ylabel('X(t)');
legend('WITHOUT observer','WITH observer');
figure;plot(Time,X(3,:),Time1,Xh(3,:));title('Teta');
xlabel('Time');ylabel('Teta(t)');
legend('WITHOUT observer','WITH observer');
figure;plot(Time(1:end-1),u);title('Control Effort');
xlabel('Time');ylabel('u(t)');
%% FUNCTIONS
function [A,B,C,D]=inv_pendulum_ss()
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

    A = vpa(A,6);B = vpa(B,6);
    A = double(A);B = double(B);
    C = [1 0 0 0];D = 0;
end