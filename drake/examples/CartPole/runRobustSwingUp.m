close all

p = CartPolePlant();
v = CartPoleVisualizer(p);

N = 81;
D = (1/2^2); %This corresponds to +/- 2 Newtons of disturbance force
Q = blkdiag((1/.01^2)*eye(2), (1/.1^2)*eye(2));
R = (1/10^2);
Qf = (1/.001^2)*eye(4);

%Standard DIRTRAN
[utraj1,xtraj1]=swingUpTrajectory(p);

%Robust DIRTRAN
[utraj2,xtraj2,z,prog] = p.robustSwingUpTrajectory(N,D,Q,R,Qf);

%Closed-loop simulations
p.friction_on = 0;
c1 = tvlqr(p,xtraj1,utraj1,Q,R,Qf);
p.friction_on = 1;
clsys = feedback(c1,p);
[~,xcl1] = clsys.simulate([utraj1.tspan(1) utraj1.tspan(2)], xtraj1.eval(0));
v.playback(xcl1);

p.friction_on = 0;
c2 = tvlqr(p,xtraj2,utraj2,Q,R,Qf);
p.friction_on = 1;
clsys = feedback(c2,p);
[~,xcl2] = clsys.simulate([utraj2.tspan(1) utraj2.tspan(2)+1], xtraj2.eval(0));
v.playback(xcl2);

%Plots
figure(2);
t = xtraj1.tspan(1):.05:xtraj1.tspan(2);
xhist = xtraj1.eval(t);
xhist_cl = xcl1.eval(t);
subplot(2,1,1);
plot(t,xhist(1,:));
hold on
plot(t,xhist_cl(1,:));
subplot(2,1,2);
plot(t,xhist(2,:));
hold on
plot(t,xhist_cl(2,:));

figure(3);
uhist = utraj1.eval(t);
plot(t,uhist);
