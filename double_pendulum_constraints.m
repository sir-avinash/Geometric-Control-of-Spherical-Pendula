function [qd2,wd2,dwd2]=double_pendulum_constraints(q1,w1,dw1,pert)

R = [1 0 0;0 1 0;0 0 -1];  %% virtual constraints

if nargin<4
    pert = eye(3);
end

dq1 = cross2(w1,q1);
qd2 = pert*(R*q1);
dqd2 = pert*R*dq1; 
wd2 = cross2(qd2,dqd2);

q1dd = cross2(dw1,q1)+cross2(w1,cross2(w1,q1));
dwd2 = cross2(qd2,pert*R*q1dd); 

end