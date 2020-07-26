clc;
clear
close all
%% System Equations
global M m l g
M = 5;m = 1;
l = 0.5;g = 9.81;
[A,B,C,D]=State_Space();
h = 0.2;
G = expm(A*h);
%syms tav
%H = int(G,tav,[0 h])*B; H = vpa(H,4);
H=(h*eye(4)+1/2*A*h^2+1/6*A^2*h^3+1/24*A^3*h^4+1/120*A^4*h^5)*B;
rank(ctrb(G,H))

des_poles = [0.5 0.5 -0.5 -0.5];
K=acker(G,H,des_poles);

T = 10;
dt = 0.0001;t=0;
k = 1;
N = floor(h/dt);
X(:,k) = [0;0;2*3.14/180;-1*3.14/180];
%% LINIEAR SYSTEM
% while t<T
%     if mod(k,N)==1
%        u(k)=-K*X(:,k);
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
       u(k)=-K*X(:,k);
    else
       u(k)=u(k-1);
    end
    SYS =[ X(2,k);
          1/(m+M-m*cos(X(3,k)))*(u(k)-m*l*X(4,k)^2*sin(X(3,k))-m*g*sin(X(3,k)));
          X(4,k);
         (g*sin(X(3,k))-(1/(m+M-m*cos(X(3,k)))*(u(k)-m*l*(X(4,k))^2*sin(X(3,k))-m*g*sin(X(3,k))))*cos(X(3,k)))/l];
    
    X(:,k+1)=X(:,k)+ SYS*dt;
    k = k+1;
    t=t+dt;    
end

%% PLOTS
Time=0:dt:T+dt;
subplot(2,2,1);plot(Time,X(1,:));title('X');
xlabel('Time');ylabel('X(t)');
subplot(2,2,2);plot(Time,X(2,:));title('X-dot');
xlabel('Time');ylabel('X-dot(t)');
subplot(2,2,3);plot(Time,X(3,:));title('Teta');
xlabel('Time');ylabel('Teta(t)');
subplot(2,2,4);plot(Time,X(4,:));title('Teta-dot');
xlabel('Time');ylabel('Teta-dot(t))');
figure;plot(Time(1:end-1),u);title('Control Effort');
xlabel('Time');ylabel('u(t)');

figure;plot(Time,X(1,:));title('X');
xlabel('Time');ylabel('X(t)');
figure;plot(Time,X(2,:));title('X-dot');
xlabel('Time');ylabel('X-dot(t)');
figure;plot(Time,X(3,:));title('Teta');
xlabel('Time');ylabel('Teta(t)');
figure;plot(Time,X(4,:));title('Teta-dot');
xlabel('Time');ylabel('Teta-dot(t)');
