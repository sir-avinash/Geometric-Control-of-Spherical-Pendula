function dx = spherical_pendulum_dynamics(t,x)

d.g = 9.8;
d.e3 = [0 0 1]';
d.m = 1;
d.l = 0.2;

q = x(1:3);
w = x(4:6);

dq = cross(w,q);
u = spherical_pendulum_controller(t,q,dq,w,d);
dw = (d.g/d.l)*cross(q,d.e3) + (1/(d.m*d.l^2))*u;

dx =[dq;dw];
% dx = dx';
end