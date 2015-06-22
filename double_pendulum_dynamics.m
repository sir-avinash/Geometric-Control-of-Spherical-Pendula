function dx = double_pendulum_dynamics(t,x)

d.g = -9.8;
d.e3 = [0 0 1]';
d.m = 1;
d.l = 0.2;

m1 = d.m;
m2 = d.m;
l1 = d.l;
l2 = d.l;
g=d.g;
e3=d.e3;

q1 = x(1:3);
q2 = x(4:6);
w1 = x(7:9);
w2 = x(10:12);

dq1 = cross2(w1,q1);
dq2 = cross2(w2,q2);

J = [(m1+m2)*l1^2*eye(3) -m2*l1*l2*hat(q1)*hat(q2);-m2*l1*l2*hat(q2)*hat(q1) m2*l2^2*eye(3)];
C = [-m2*l1*l2*norm(w2)^2*hat(q1)*q2;-m2*l1*l2*norm(w1)^2*hat(q2)*q1];
G = [(m1+m2)*g*l1*hat(q1)*e3;m2*g*l1*hat(q2)*e3];

[u1,u2] = double_pendulum_controller(t,q1,q2,dq1,dq2,w1,w2,J,G,C);
U =  [u1;u2];% [hat(q1)*u1;hat(q2)*u2];

q2'*u2

dw = J\(U-G-C);

dx =[dq1;dq2;dw];
% dx = dx';
end