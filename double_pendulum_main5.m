function double_pendulum_main5

addpath('C:\Users\Avinash\Dropbox\CMU_lap\Misc\My Toolbox');
d.g = 9.8;
d.e3 = [0 0 1]';
d.m = 1;
d.l = 0.2;

q10 = Rx(-50*pi/180)*Ry(-70*pi/180)*d.e3;
w10 = [0 -2.5 0]';
w10 = trans_map(q10,w10);
pert = Rx(0*pi/180)*Ry(0*pi/180);

[q20,w20]=double_pendulum_perturbation(q10,w10,pert);

x0 = [q10;q20;w10;w20];
options = odeset('RelTol',1e-9,'AbsTol',1e-12); 

[T,X]=ode45(@double_pendulum_dynamics,[0 1],x0,options,d);
 Fs=120;
 [T,X] = even_sample(T,X,Fs);

figure
for i = 1:size(X,1)
    q1 = X(i,1:3)';
    q2 = X(i,4:6)';
    w1 = X(i,7:9)';
    w2 = X(i,10:12)';

hplot = plot3(0,0,0,'r-o');
set(hplot,'XDATA',[0,q1(1),NaN,q1(1),q1(1)+q2(1),NaN],'YDATA',[0,q1(2),NaN,q1(2),q1(2)+q2(2),NaN],'ZDATA',[0,q1(3),NaN,q1(3),q1(3)+q2(3),NaN]); 
hold on
q2d=double_pendulum_constraints(q1,w1);
hplot2 = plot3(0,0,0,'b-x');
set(hplot2,'XDATA',[0,q1(1),NaN,q1(1),q1(1)+q2d(1),NaN],'YDATA',[0,q1(2),NaN,q1(2),q1(2)+q2d(2),NaN],'ZDATA',[0,q1(3),NaN,q1(3),q1(3)+q2d(3),NaN]); 
hold off 
view(-45,75)
axis([-2 2 -2 2 -2 2])

drawnow

%% to plot error dynamics %%%
[J,C,G,B]=double_pendulum_model(q1,q2,w1,w2,d);
dq2=hat(w2)*q2;
[u,eq2c,ew2c] = double_pendulum_controller(q1,q2,dq2,w1,w2,J,G,C,B);
eq2(:,i) = eq2c;
ew2(:,i) = ew2c;

normeq(i) = norm(eq2(:,i));
normew(i) = norm(ew2(:,i));
doteq(i) = dot(eq2(:,i),q2);
dotew(i) = dot(ew2(:,i),q2);
end

% Error Dynamics (eq) and (ew)
%%% Q Error plots for swing leg (eq2) and trunk (eq3)
figure(3)
% plot(T',eq2')
plot(T,normeq,T,doteq);
axis([0 1 -1 1])
 axis('tight');
xlabel('time');
ylabel('eq2');
% legend('X','Y','Z');
legend('norm(eq_2)','dot(q_2,eq_2)');
%%%% Omega error plots for swing leg (ew2) and trunk (ew3)
figure(5)
% plot(T',ew2')
plot(T,normew,T,dotew);
axis([0 1 -1 1])
axis('tight');
xlabel('time');
ylabel('ew2');
% legend('X','Y','Z');
legend('norm(ew_2)','dot(q_2,ew_2)');

end

function dx = double_pendulum_dynamics(t,x,d)

q1 = x(1:3);
q2 = x(4:6);
w1 = x(7:9);
w2 = x(10:12);

dq1 = cross2(w1,q1);
dq2 = cross2(w2,q2);

[J,C,G,B]=double_pendulum_model(q1,q2,w1,w2,d);
u = double_pendulum_controller(q1,q2,dq2,w1,w2,J,G,C,B);
dw = J\(B*u-G-C);

dx =[dq1;dq2;dw];
end

function [u,eq2,ew2] = double_pendulum_controller(q1,q2,dq2,w1,w2,J,G,C,B)
dq1=hat(w1)*q1;

Fq = J\(-C-G);
Gq = J\B;

Fq1 = Fq(1:3); Fq2 = Fq(4:6); 
Gq1 = Gq(1:3,:); Gq2 = Gq(4:6,:);

%%%%%% Defining error dynamics
tmp=-hat(q2)^2; 
eps = 0.5;
kq2 = 10/eps^2;
kw2 = 10.1/eps;
[q2d,w2d,R]=double_pendulum_constraints(q1,w1);
eq2 = cross2(q2d,q2);
ew2 = w2 - tmp*w2d; 

v=-kq2*eq2-kw2*ew2;
Lf2h = tmp*Fq2 - tmp*(hat(q2d)^2*Fq1) - tmp*hat(w2d)*w1 + hat(q2)*hat(dq2)*w2d + 0*hat(dq2)*hat(q2)*w2d;
LgLfh = tmp*Gq2 - tmp*(hat(q2d)^2*Gq1);
% Lf2h = Fq2 + tmp*R*Fq1 + dot(q2,w2d)*dq2 + dot(dq2,w2d)*q2;
% LgLfh = Gq2 + tmp*R*Gq1;
u = pinv(LgLfh)*(0*v-Lf2h);

% dot(q2,Fq2)
% % -norm(hat(q2d)^2*Gq1)-norm(R*Gq1)
end

function [q20,w20]=double_pendulum_perturbation(q10,w10,pertq,pertw)
if nargin<4
    pertw = eye(3);
end
[q20,w20,R]=double_pendulum_constraints(q10,w10);
q20=pertq*q20;
w20=pertw*w20;
w20 = trans_map(q20,w20);
end

function [q2d,w2d,R]=double_pendulum_constraints(q1,w1)
R = [1 0 0;0 1 0;0 0 -1];  %% virtual constraints
dq1 = cross2(w1,q1);
q2d = (R*q1);
dq2d = R*dq1; 
w2d = cross2(q2d,dq2d);
% w2d = R*w1;
end

function [J,C,G,B]=double_pendulum_model(q1,q2,w1,w2,d)
m1 = d.m;
m2 = d.m;
l1 = d.l;
l2 = d.l;
g=d.g;
e3=d.e3;
J = [(m1+m2)*l1^2*eye(3) -m2*l1*l2*hat(q1)*hat(q2);-m2*l1*l2*hat(q2)*hat(q1) m2*l2^2*eye(3)];
C = [-m2*l1*l2*norm(w2)^2*hat(q1)*q2;-m2*l1*l2*norm(w1)^2*hat(q2)*q1];
G = [(m1+m2)*g*l1*hat(q1)*e3;m2*g*l1*hat(q2)*e3];
B= [zeros(3);eye(3)];
end

function [Et, Ex] = even_sample(t, x, Fs)

% Obtain the process related parameters
N = size(x, 2);    % number of signals to be interpolated
M = size(t, 1);    % Number of samples provided
t0 = t(1,1);       % Initial time
tf = t(M,1);       % Final time
EM = (tf-t0)*Fs;   % Number of samples in the evenly sampled case with
% the specified sampling frequency
Et = linspace(t0, tf, EM)';

% Using linear interpolation (used to be cubic spline interpolation)
% and re-sample each signal to obtain the evenly sampled forms
for s = 1:N,
	Ex(:,s) = interp1(t(:,1), x(:,s), Et(:,1));
end
end
