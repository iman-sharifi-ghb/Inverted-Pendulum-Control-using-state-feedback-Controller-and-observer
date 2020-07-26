function Animate_Inv_Pdlm(X,t)
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
