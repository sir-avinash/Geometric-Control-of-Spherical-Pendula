function double_pendulum_main
% R = [1 0 0;0 1 0;0 0 -1];
% q10 = [0 0 -1]';
% q20 = R*q10;%[0 0 1]';
% 
% w10 = [0 -2 0]';
% % w20 = [0  0]';
% dq10 = cross2(w10,q10);
% w20 = cross2(R*q10,R*dq10);

q10 = Rx(-1*pi/12)*Ry(-1*pi/8)*[0 0 1]';
w10 = [0 1 0]';
pert = Rx(0*pi/12)*Ry(1*pi/8);%inv([1 0 0;0 1 0;0 0 -1]);
[q20,w20]=double_pendulum_perturbation(q10,w10,pert);
x0 = [q10;q20;w10;w20];

[T,X]=ode45(@double_pendulum_dynamics,[0 1],x0);


 
% alpha = 0.5;
% beta = 2;
figure
for i = 1:size(X,1)
    q1 = X(i,1:3)';
    q2 = X(i,4:6)';
    w1 = X(i,7:9)';
    w2 = X(i,10:12)';
%     ta = T(i);
%     q1d = [cos(alpha*ta)*cos(beta*ta) cos(alpha*ta)*sin(beta*ta) sin(alpha*ta)]';
%     q2d = [cos(alpha*ta)*cos(beta*ta) cos(alpha*ta)*sin(beta*ta) sin(alpha*ta)]';

% [xs,ys,zs]=sphere;
%  alpha = 0.01;
%  surf(2*xs,2*ys,2*zs,'FaceColor',[0 0 0],'FaceAlpha',alpha);
%  colormap([0 0 0])
%  hold on
hplot = plot3(0,0,0,'r-o');
set(hplot,'XDATA',[0,q1(1),NaN,q1(1),q1(1)+q2(1),NaN],'YDATA',[0,q1(2),NaN,q1(2),q1(2)+q2(2),NaN],'ZDATA',[0,q1(3),NaN,q1(3),q1(3)+q2(3),NaN]); 
hold on
[qd2,wd2]=double_pendulum_constraints(q1,w1,zeros(3,1));
hplot2 = plot3(0,0,0,'b-x');
set(hplot2,'XDATA',[0,q1(1),NaN,q1(1),q1(1)+qd2(1),NaN],'YDATA',[0,q1(2),NaN,q1(2),q1(2)+qd2(2),NaN],'ZDATA',[0,q1(3),NaN,q1(3),q1(3)+qd2(3),NaN]); 
hold off 
axis([-2 2 -2 2 -2 2])
drawnow

%% to plot error dynamics %%%

% qd2 = R*q1;
% dq1 = cross2(w1,q1);
% wd2 = cross2(R*q1,R*dq1);
eq2(:,i) = cross2(qd2,q2);
ew2(:,i) = w2 + (hat(q2)^2)*wd2;


end

%% Error Dynamics (eq) and (ew)
%%%% Q Error plots for swing leg (eq2) and trunk (eq3)
figure(3)
plot(T',eq2')
axis([0 1 -1 1])
% axis('tight');
xlabel('time');
ylabel('eq2');
legend('X','Y','Z');

%%%%% Omega error plots for swing leg (ew2) and trunk (ew3)
figure(5)
plot(T',ew2')
axis([0 1 -1 1])
% axis('tight');
xlabel('time');
ylabel('ew2');
legend('X','Y','Z');

end




% hold off

% figure(2)
% plot(T,X(:,4),T,X(:,5),T,X(:,6));


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

% q2'*u2

dw = J\(U-G-C);
dx =[dq1;dq2;dw];
% dx = dx';
end

function [u1,u2] = double_pendulum_controller(t,q1,q2,dq1,dq2,w1,w2,J,G,C)

% % nargin
% % kq2 = 1050;
% % kw2 = 800;

% Double Input Case
eps = 1/sqrt(25);
kq2 = 10/eps^2;
kw2 = 10.1/eps^2;
% % 
% % 
% % alpha = 0.5;
% % beta = 2;
% % 
% % % syms ta
% % % qd1 = [cos(alpha*ta)*cos(beta*ta); cos(alpha*ta)*sin(beta*ta); sin(alpha*ta)];
% % % qd2 = [cos(alpha*ta)*cos(beta*ta); cos(alpha*ta)*sin(beta*ta); sin(alpha*ta)];
% % 
% % % dqd1 = diff(qd1,ta);
% % % wd1 = simplify(cross2(qd1,dqd1));
% % % ddqd1 = diff(dqd1,ta);
% % % dwd1 = simplify(cross2(qd1,ddqd1));
% % 
% % % qd = double(subs(qd,ta,t));
% % % wd = double(subs(wd,ta,t));
% % % dwd = double(subs(dwd,ta,t));
% % % 
% %  qd1 = [cos(alpha*t)*cos(beta*t); cos(alpha*t)*sin(beta*t); sin(alpha*t)];
% %  qd2 = [cos(alpha*t)*cos(beta*t); cos(alpha*t)*sin(beta*t); sin(alpha*t)];
% %  wd1 = [sin(2*t)/2 - sin(3*t)/2 + sin(t)/2;20*sin(t/2)^4 - 4*sin(t/2)^2 - 16*sin(t/2)^6 - 1/2;cos(t) + 1];
% %  wd2 = [sin(2*t)/2 - sin(3*t)/2 + sin(t)/2;20*sin(t/2)^4 - 4*sin(t/2)^2 - 16*sin(t/2)^6 - 1/2;cos(t) + 1];
% %  dwd1 = [2*sin(t/2)^2*(24*sin(t/2)^4 - 32*sin(t/2)^2 + 9);
% %              sin(2*t) - (3*sin(3*t))/2 + sin(t)/2;
% %                                             -sin(t)];
% % dwd2 = [2*sin(t/2)^2*(24*sin(t/2)^4 - 32*sin(t/2)^2 + 9);
% %              sin(2*t) - (3*sin(3*t))/2 + sin(t)/2;
% %                                             -sin(t)];
% % 

%% Single Input Case

% R = [1 0 0;0 1 0;0 0 -1];
% qd2 = q1; %R*q1;%Ry(-pi/3)*Rx(-pi/4)*q1; %Ry(-pi/3)*Rx(-pi/8)*[0 0 1]'; %q1;
% % wd2 = w1;
% 
% % qd2 = R2*q1; %[0;0;1];%
% dq1 = cross2(w1,q1);
% wd2 = cross2(R*q1,R*dq1);

dw = J\(-G-C);
dw1 = dw(4:6);  
% ddq1 = cross2(dw1,q1)+cross2(w1,cross2(w1,q1));
% dwd2 = cross2(R*q1,R*ddq1); 
% dwd2 = dw1;
[qd2,wd2,dwd2]=double_pendulum_constraints(q1,w1,dw1);
eq2 = cross2(qd2,q2);
% eq2 = cross2(qd2,q2);
ew2 = w2 + (hat(q2)^2)*wd2;
% ew2 = w2 + (hat(q2)^2)*wd2;
%v = -[kw1*eye(3) zeros(3);zeros(3) kw2*eye(3)]*[ew1;ew2] -[kq1*eye(3) zeros(3);zeros(3) kq2*eye(3)]*[eq1;eq2] - [(hat(q1)^2)*dwd1;(hat(q2)^2)*dwd2] -[dq1*dot(q1,wd1);dq2*dot(q2,wd2)];

%% single input case
v2 = -kw2*ew2 - kq2*eq2 ;
v1 = zeros(3,1);
v =[v1;v2];
extra = dwd2 - dq2*dot(q2,wd2)+ 0*q2*(dot(q2,dwd2)+dot(dq2,wd2));

u = J*(1*v + 1*[zeros(3,1);extra] +(J\G) + (J\C));
u1 = u(1:3);
u2 = u(4:6);

end


function [q20,w20]=double_pendulum_perturbation(q10,w10,pertq,pertw)
if nargin<4
    pertw = eye(3);
end
dw1=zeros(3,1);
[q20,w20]=double_pendulum_constraints(q10,w10,dw1);
q20=pertq*q20;
w20=pertw*w20;
end

function [qd2,wd2,dwd2]=double_pendulum_constraints(q1,w1,dw1)
R = [1 0 0;0 1 0;0 0 -1];  %% virtual constraints
dq1 = cross2(w1,q1);
qd2 = (R*q1);
dqd2 = R*dq1; 
wd2 = cross2(qd2,dqd2);
ddq1 = cross2(dw1,q1)+cross2(w1,cross2(w1,q1));
dwd2 = cross2(qd2,R*ddq1); 
end
