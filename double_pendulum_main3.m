function double_pendulum_main3

%%%%% Defining output y = q-H(q) => q2=R*q1   %%%%%%%

addpath('C:\Users\Avinash\Dropbox\CMU_lap\Misc\My Toolbox');

%%%%%% Declaring Constants %%%%%%%%%%
d.g = 9.8;
d.e3 = [0 0 1]';
d.m = 1;
d.l = 0.2;
d.R = [1 0 0;0 1 0;0 0 -1];

%%%%%% Initial Conditions %%%%%%%%%%%

q10 = Rx(-0*pi/180)*Ry(-20*pi/180)*[0 0 1]';
w10 = [0 -2 0]';
dq10 = hat(w10)*q10;

pertq = Rx(0*pi/180)*Ry(0*pi/180);
[q20,dq20]=double_pendulum_initial(q10,w10,d,pertq);

x0 = [q10;q20;dq10;dq20];
options = odeset('RelTol',1e-9,'AbsTol',1e-12); 
[T,X]=ode45(@double_pendulum_dynamics,[0 1],x0,options,d);


figure
for i = 1:size(X,1)
    q1 = X(i,1:3)';
    q2 = X(i,4:6)';
    w1 = X(i,7:9)';
    w2 = X(i,10:12)';
    dq1 = hat(w1)*q1;
    dq2 = hat(w2)*q2;

hplot = plot3(0,0,0,'r-o');
set(hplot,'XDATA',[0,q1(1),NaN,q1(1),q1(1)+q2(1),NaN],'YDATA',[0,q1(2),NaN,q1(2),q1(2)+q2(2),NaN],'ZDATA',[0,q1(3),NaN,q1(3),q1(3)+q2(3),NaN]); 
hold on

%%%%%  the desired configurations
[q2d,dq2d]=double_pendulum_constraints(q1,w1,d);

hplot2 = plot3(0,0,0,'b-x');
set(hplot2,'XDATA',[0,q1(1),NaN,q1(1),q1(1)+q2d(1),NaN],'YDATA',[0,q1(2),NaN,q1(2),q1(2)+q2d(2),NaN],'ZDATA',[0,q1(3),NaN,q1(3),q1(3)+q2d(3),NaN]); 
hold off 
axis([-2 2 -2 2 -2 2])
drawnow

% %% to plot error dynamics %%%
y(:,i) = q2-q2d;
dy(:,i) = dq2-dq2d;

[norm(y(:,i)), norm(dy(:,i))]
% eq2(:,i) = cross2(q2d,q2);
% ew2(:,i) = w2 + (hat(q2)^2)*w2d;
% % [norm(eq2(:,i)),norm(ew2(:,i))] 
% [norm(w1),norm(w2)]
end

% Output Dynamics (y) and (dy)
%%% y %%%%%%%%
figure(3)
plot(T',y')
% axis([0 1 -1 1])
% axis('tight');
xlabel('time');
ylabel('y');
legend('X','Y','Z');

%%%% dy %%%%%%
figure(5)
plot(T',dy')
% axis([0 1 -1 1])
% axis('tight');
xlabel('time');
ylabel('dy');
legend('X','Y','Z');

% % Error Dynamics (eq) and (ew)
% %%% Q Error plots for swing leg (eq2) and trunk (eq3)
% figure(3)
% plot(T',eq2')
% axis([0 1 -1 1])
%  axis('tight');
% xlabel('time');
% ylabel('eq2');
% legend('X','Y','Z');
% 
% %%%% Omega error plots for swing leg (ew2) and trunk (ew3)
% figure(5)
% plot(T',ew2')
% axis([0 1 -1 1])
% axis('tight');
% xlabel('time');
% ylabel('ew2');
% legend('X','Y','Z');

end

function dx = double_pendulum_dynamics(t,x,d)



m1 = d.m;
m2 = d.m;
l1 = d.l;
l2 = d.l;
g=d.g;
e3=d.e3;
R=d.R;

q1 = x(1:3);
q2 = x(4:6);
w1 = x(7:9);
w2 = x(10:12);

dq1 = cross2(w1,q1);
dq2 = cross2(w2,q2);

J = [(m1+m2)*l1^2*eye(3) -m2*l1*l2*hat(q1)*hat(q2);-m2*l1*l2*hat(q2)*hat(q1) m2*l2^2*eye(3)];
C = [-m2*l1*l2*norm(w2)^2*hat(q1)*q2;-m2*l1*l2*norm(w1)^2*hat(q2)*q1];
G = [(m1+m2)*g*l1*hat(q1)*e3;m2*g*l1*hat(q2)*e3];

B = [zeros(3);eye(3)];

u = double_pendulum_controller(q1,q2,dq1,dq2,J,G,C,B,d);


% U =  [zeros(3,1);u]; % [hat(q1)*u1;hat(q2)*u2];
dw = J\(B*u-G-C);
% dw1 = J\(-G-C);
% % norm(dw-dw1)
dx =[dq1;dq2;dw];


%%%%%%%   Testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rank(inv(J))
% q2'*u2
%test = [rank(inv(J(1:3,1:3))),rank(inv(J(1:3,4:6))),rank(inv(J(4:6,1:3))),rank(inv(J(4:6,4:6)))]
% rank(J(1:3,4:6))
% dw1=dw(1:3,1);
% norm(dw-[dw1;dw2d])
% dw2 = J\(-G-C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

function u = double_pendulum_controller(q1,q2,dq1,dq2,J,G,C,B,d)

R=d.R;
w1=hat(q1)*dq1;
w2=hat(q2)*dq2;

%%%% Partial Feedback Linearation  %%%%%%
Fq = J\(-C-G);
Gq = J\B;

Fq1 = Fq(1:3); Fq2 = Fq(4:6); 
Gq1 = Gq(1:3,:); Gq2 = Gq(4:6,:);

[q2d,dq2d]=double_pendulum_constraints(q1,w1,d);

y = q2-q2d;
tm = eye(3);%(-hat(q2)^2);
dy = dq2-tm*dq2d; %Lfh
[y,dy,dq2+(hat(q2)^2)*dq2d] 
Lf2h = -hat(q2)*Fq2 + tm*R*hat(q1)*Fq1 + hat(w2)^2*q2 + tm*R*hat(dq1)*w1 -0*(hat(dq2)*hat(q2)*R*hat(q1)*w1 + hat(q2)*hat(dq2)*R*hat(q1)*w1);
LgLfh = -hat(q2)*Gq2 + tm*R*hat(q1)*Gq1;
eps=1/sqrt(1);
Ky = 10/eps^2;
Kdy = 10.1/eps;
v = -Ky*y -Kdy*dy;
u = pinv(LgLfh)*(0*v-Lf2h);

end

function [q20,dq20,w20]=double_pendulum_initial(q10,w10,d,pertq,pertw)
if nargin<5
    pertw = eye(3);
end
[q20,dq20,w20]=double_pendulum_constraints(q10,w10,d);
q20=pertq*q20;
w20=pertw*w20;
end

function [q2d,dq2d,w2d]=double_pendulum_constraints(q1,w1,d)
R = d.R;  %% virtual constraints
dq1 = cross2(w1,q1);
q2d = (R*q1);
dq2d = R*dq1; 
w2d = cross2(q2d,dq2d);
end

