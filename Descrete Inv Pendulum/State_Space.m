function [A,B,C,D]=State_Space()
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

