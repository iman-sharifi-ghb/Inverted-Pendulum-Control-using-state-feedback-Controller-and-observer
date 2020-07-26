clc;
clear;close all
%   X(k+1)=GX(k)+Hu(k)
%    y(k)=CX(k)+Du(k)
%%  Servo Type 1 + Observer
global M m l g
M = 3;m = 1;
l = 0.5;g = 9.81;

[A,B,C,D]=State_Space();
h = 0.2;
G = expm(A*h);
%syms tav
% H = int(G,tav,[0 h])*B; H = vpa(H,4);
H=(h*eye(4)+1/2*A*h^2+1/6*A^2*h^3+1/24*A^3*h^4+1/120*A^4*h^5)*B;
rank(ctrb(G,H))
rank(obsv(G,C))

rank(obsv(G,C))
des_poles = [0.09 -0.11 0.11 -0.09];
L=acker(G',C',des_poles)';

mu=[0 0 0 0];
K=acker(G,H,mu);
eig(G-H*K)
Gcl_0 = -C/(G-H*K)*H; 
X(:,1)=0.1*[0.5;0;5*pi/180;-1*pi/180];
N=100;
dt = 0.0001;
%% DIGITAL SYSTEM
% for k=1:N-1
%     yd(k)=2*sign(sin(0.2*k));
%     u(k)=-K*X(:,k)+yd(k)/Gcl_0;
%     X(:,k+1)=G*X(:,k)+H*u(k);
% end
%% LINEAR SYSTEM
% for k=1:N-1
%     yd(k)=2*sign(sin(0.2*k));
%     u(k)=-K*X(:,k)+yd(k)/Gcl_0;
%     X(:,k+1)=X(:,k)+(A*X(:,k)+B*u(k))*dt;
% end
%% NONLINEAR SYSTEM
for k=1:N-1
    yd(k)=2*sign(sin(0.2*k));
    u(k)=-K*X(:,k)+yd(k)/Gcl_0;
    SYS =[ X(2,k);
          1/(m+M-m*cos(X(3,k)))*(u(k)-m*l*X(4,k)^2*sin(X(3,k))-m*g*sin(X(3,k)));
          X(4,k);
         (g*sin(X(3,k))-(1/(m+M-m*cos(X(3,k)))*(u(k)-m*l*(X(4,k))^2*sin(X(3,k))-m*g*sin(X(3,k))))*cos(X(3,k)))/l];
    X(:,k+1)=X(:,k)+(SYS)*dt;
end
%% PLOTS
Time=1:N;
% X(1,:)=X(1,:)+1.93*sign(sin(0.202*Time));
plot(Time,X(1,:),'b',Time(1:N-1),yd,'r');title('X');
xlabel('Time');ylabel('X(t)');
figure;
% u=u+12*sign(sin(0.3.*Time(1:length(u))));
plot(u);title('Control Effort');
xlabel('Time');ylabel('u(t)');
