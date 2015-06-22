function [u1,u2] = double_pendulum_controller(t,q1,q2,dq1,dq2,w1,w2,J,G,C)


% kq2 = 1050;
% kw2 = 800;

% Double Input Case
kq2 = 100;
kw2 = 100.1;
% 
% 
% alpha = 0.5;
% beta = 2;
% 
% % syms ta
% % qd1 = [cos(alpha*ta)*cos(beta*ta); cos(alpha*ta)*sin(beta*ta); sin(alpha*ta)];
% % qd2 = [cos(alpha*ta)*cos(beta*ta); cos(alpha*ta)*sin(beta*ta); sin(alpha*ta)];
% 
% % dqd1 = diff(qd1,ta);
% % wd1 = simplify(cross2(qd1,dqd1));
% % ddqd1 = diff(dqd1,ta);
% % dwd1 = simplify(cross2(qd1,ddqd1));
% 
% % qd = double(subs(qd,ta,t));
% % wd = double(subs(wd,ta,t));
% % dwd = double(subs(dwd,ta,t));
% % 
%  qd1 = [cos(alpha*t)*cos(beta*t); cos(alpha*t)*sin(beta*t); sin(alpha*t)];
%  qd2 = [cos(alpha*t)*cos(beta*t); cos(alpha*t)*sin(beta*t); sin(alpha*t)];
%  wd1 = [sin(2*t)/2 - sin(3*t)/2 + sin(t)/2;20*sin(t/2)^4 - 4*sin(t/2)^2 - 16*sin(t/2)^6 - 1/2;cos(t) + 1];
%  wd2 = [sin(2*t)/2 - sin(3*t)/2 + sin(t)/2;20*sin(t/2)^4 - 4*sin(t/2)^2 - 16*sin(t/2)^6 - 1/2;cos(t) + 1];
%  dwd1 = [2*sin(t/2)^2*(24*sin(t/2)^4 - 32*sin(t/2)^2 + 9);
%              sin(2*t) - (3*sin(3*t))/2 + sin(t)/2;
%                                             -sin(t)];
% dwd2 = [2*sin(t/2)^2*(24*sin(t/2)^4 - 32*sin(t/2)^2 + 9);
%              sin(2*t) - (3*sin(3*t))/2 + sin(t)/2;
%                                             -sin(t)];
% 

%% Single Input Case

R = [1 0 0;0 1 0;0 0 -1];
qd2 = q1; %R*q1;%Ry(-pi/3)*Rx(-pi/4)*q1; %Ry(-pi/3)*Rx(-pi/8)*[0 0 1]'; %q1;
% wd2 = w1;

% qd2 = R2*q1; %[0;0;1];%
dq1 = cross2(w1,q1);
wd2 = cross2(R*q1,R*dq1);

dw = J\(-G-C);
dw1 = dw(4:6);
q1dd = cross2(dw1,q1)+cross2(w1,cross2(w1,q1));
dwd2 = cross2(R*q1,R*q1dd); 
% dwd2 = dw1;

eq2 = cross2(qd2,q2);
% eq2 = cross2(qd2,q2);
ew2 = w2 + (hat(q2)^2)*wd2;
% ew2 = w2 + (hat(q2)^2)*wd2;
%v = -[kw1*eye(3) zeros(3);zeros(3) kw2*eye(3)]*[ew1;ew2] -[kq1*eye(3) zeros(3);zeros(3) kq2*eye(3)]*[eq1;eq2] - [(hat(q1)^2)*dwd1;(hat(q2)^2)*dwd2] -[dq1*dot(q1,wd1);dq2*dot(q2,wd2)];

%% single input case
v2 = -kw2*ew2 - kq2*eq2 - (hat(q2)^2)*dwd2 - dq2*dot(q2,wd2);
v1 = zeros(3,1);
v =[v1;v2];

u = J*(1*v + (J\G) + (J\C));
u1 = u(1:3);
u2 = u(4:6);

end

