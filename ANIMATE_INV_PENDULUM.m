clc;
clear;close all

a = 2;
b = 1;
xg = 0;
yg = 0;
r = a/8;
L = 3;
t = pi/6;
R = 2*r;
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
teta = 0:0.01:2*pi;
x5 = x5+r*cos(teta);y5 = y5+r*sin(teta);
x6 = x6+r*cos(teta);y6 = y6+r*sin(teta);
%% pendulum
x7 = xg;y7=yg+b/2;
x8 = x7+L*sin(t);y8=y7+L*cos(t);
teta = 0:0.01:2*pi;
x9 = x8+R*cos(teta);y9 = y8+R*sin(teta);

plot(x1,y1,'b-',x2,y2,'b-',x3,y3,'b-',x4,y4,'b-',x5,y5,'g-',x6,y6,'g-',...
    x7,y7,'bo',[x7 x8],[y7 y8],'g-',x9,y9,'r-','linewidth',3)
axis([-2*a 2*a -1*a 2*a])
axis equal
% text(-2*a, 2*a,['time= ' num2str(t(i)) ' (s)'])
% pause (.1)
% hold off

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
        x5 = x5+r*cos(t1);y5 = y5+r*sin(t1);
        x6 = x6+r*cos(t1);y6 = y6+r*sin(t1);
        %% pendulum
        x7 = xg;y7=yg+b/2;
        x8 = x7+L*sin(teta);y8=y7+L*cos(teta);
        x9 = x8+R*cos(t1);y9 = y8+R*sin(t1);
        %% plot
        plot(x1,y1,'b-',x2,y2,'b-',x3,y3,'b-',x4,y4,'b-',x5,y5,'g-',x6,y6,'g-',...
            x7,y7,'bo',[x7 x8],[y7 y8],'g-',x9,y9,'r-','linewidth',3)
        axis([-2*a 2*a -1*a 2*a])
        axis equal
        text(-2*a, 2*a,['time= ' num2str(t(i)) ' (s)'])
        pause (.1)
%         hold off
    end
end