function u = spherical_pendulum_controller(t,q,dq,w,d)
kq = 10;
kw = 10.1;

alpha = 0.5;
beta = 2;

syms ta
qd = [cos(alpha*ta)*cos(beta*ta); cos(alpha*ta)*sin(beta*ta); sin(alpha*ta)];
dqd = diff(qd,ta);
wd = cross(qd,dqd);
ddqd = diff(dqd,ta);
dwd = cross(qd,ddqd);

qd = double(subs(qd,ta,t));
wd = double(subs(wd,ta,t));
% dqd = double(subs(dqd,ta,t));
dwd = double(subs(dwd,ta,t));


eq = cross(qd,q);
ew = w + (hat(q)^2)*wd;
% u = (d.m*d.l^2)*(-kw*ew -kq*eq - (hat(q)^2)*dwd -dq*dot(q,wd)-(d.g/d.l)*cross(q,d.e3));
u = (d.m*d.l^2)*(1*(-kw*ew -kq*eq + dwd -dq*dot(q,wd))-(d.g/d.l)*cross(q,d.e3))

[dot(ew,q) dot(w,q)]

end
