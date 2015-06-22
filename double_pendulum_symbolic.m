dim = [3,1];
e3=[0 0 1]';
q1 = sym('q1',dim);
q2 = sym('q2',dim);
w1 = sym('w1',dim);
w2 = sym('w2',dim);
u1 = sym('u1',dim);
u2 = sym('u2',dim);

syms m1 m2 l1 l2 g

q1d = cross(w1,q1);
q2d = cross(w2,q2);

J = [(m1+m2)*l1^2*eye(3) -m2*l1*l2*hat(q1)*hat(q2);-m2*l1*l2*hat(q2)*hat(q1) m2*l2^2*eye(3)];
C = [-m2*l1*l2*norm(w2)^2*hat(q1)*q2;-m2*l1*l2*norm(w1)^2*hat(q2)*q1];
G = [(m1+m2)*g*l1*hat(q1)*e3;m2*g*l1*hat(q2)*e3];
U = [0*hat(q1)*u1;hat(q2)*u2];

%  to check rank of J
% Jbar = double(subs(J,{q1(1),q1(2),q1(3),q2(1),q2(2),q2(3),m1,m2,l1,l2},{0,0,1,0,0,1,1,1,1,1}))

wd = J\(U-G-C);
% wd = simplify(wd);
w1d = wd(1:3,:)
w2d = wd(4:6,:)
