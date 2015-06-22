function spherical_pendulum_main()

q0 = [0 0 1]';
w0 = [0 0 0]';
x0 = [q0;w0];

[T,X]=ode45(@spherical_pendulum_dynamics,[0 20],x0);

% [xs,ys,zs]=sphere;
figure
% meshc(xs,ys,zs)
% colormap([0 0 0])
% hold on
alpha = 0.5;
beta = 2;
for i = 1:size(X,1)
    q = X(i,1:3)';
    w = X(i,4:6)';
    ta = T(i);
    qd = [cos(alpha*ta)*cos(beta*ta) cos(alpha*ta)*sin(beta*ta) sin(alpha*ta)]';
    dqd= [-(cos(2*ta)*sin(ta/2))/2 - 2*cos(ta/2)*sin(2*ta), 2*cos(2*ta)*cos(ta/2) - (sin(2*ta)*sin(ta/2))/2, cos(ta/2)/2];
%     dqd = diff(qd,ta);
    wd = cross(qd,dqd);
    eq = cross(qd,q);
    ew = w + (hat(q)^2)*wd';
hplot = plot3(0,0,0,'r-o');
set(hplot,'XDATA',[0,q(1),NaN],'YDATA',[0,q(2),NaN],'ZDATA',[0,q(3),NaN]); 
hold on
hplot2 = plot3(0,0,0,'b-x');
set(hplot2,'XDATA',[0,qd(1),NaN],'YDATA',[0,qd(2),NaN],'ZDATA',[0,qd(3),NaN]);
hold off
axis([-1 1 -1 1 -1 1])
drawnow
cross2(q,qd)

dot1(i) = dot(w,q);
dot2(i) = dot(ew,q);
end
hold off
% hold off
figure(2)
% plot(T,X(:,4),T,X(:,5),T,X(:,6));
plot(T,dot1,T,dot2)