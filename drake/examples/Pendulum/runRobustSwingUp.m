close all

%b = 0.1; %damping constant;
%popts.torque_limit = 1;
%p = PendulumPlant(b,popts);
p = PendulumPlant();
v = PendulumVisualizer();

N = 21;
D = 50; %This corresponds to +/-.2 uncertainty in mass (20%)

options.integration_method = DirtranTrajectoryOptimization.FORWARD_EULER;
[utraj1,xtraj1,z1,prog1] = p.swingUpTrajectory(N,options);

%Plot nominal trajectory
h1 = z1(prog1.h_inds);
t1 = [0; cumsum(h1)];
x1 = z1(prog1.x_inds);
u1 = z1(prog1.u_inds);
figure(2);
subplot(2,1,1);
plot(t1,u1);
ylabel('u_{nominal}');

figure(3);
subplot(2,1,1);
plot(t1,x1(1,:));
hold on
plot(t1,x1(2,:));
ylabel('x_{nominal}');

figure(1);
[utraj2,xtraj2,z2,prog2] = p.robustSwingUpTrajectory(N,D);

%Plot robust trajectory
h2 = z2(prog2.h_inds);
t2 = [0; cumsum(h2)];
x2 = z2(prog2.x_inds);
dx = z2(prog2.dx_inds);
u2 = z2(prog2.u_inds);
du = z2(prog2.du_inds);
w = z2(prog2.w_inds);

figure(2);
subplot(2,1,2);
plot(t2,u2);
ylabel('u_{robust}');

figure(3);
subplot(2,1,2);
plot(t1,x2(1,:));
hold on
plot(t1,x2(2,:));
ylabel('x_{robust}');

figure(4);
subplot(3,1,1);
plot(t2,x2(1,:));
hold on;
plot(t2,x2(2,:));
ylabel('x');
subplot(3,1,2);
plot(t2,dx(1,:));
hold on;
plot(t2,dx(2,:));
ylabel('\delta x');
subplot(3,1,3);
plot(t2,x2(1,:)+dx(1,:));
hold on;
plot(t2,x2(2,:)+dx(2,:));
ylabel('x + \delta x');

figure(5);
subplot(3,1,1);
plot(t2,u2);
ylabel('u');
ylim([-3 3]);
subplot(3,1,2);
plot(t2,du);
ylim([-3 3]);
ylabel('\delta u');
subplot(3,1,3);
plot(t2,u2+du);
ylim([-3 3]);
ylabel('u + \delta u');

figure(6);
plot(t2,w);
ylabel('w')

% open-loop playback
% olsys = cascade(utraj2,p);
% [~,xol] = olsys.simulate([0 4], [0 0]');
% v.playback(xol);

% closed-loop
Q = [3 0; 0 1];
R = .1;
p = p.setMass(1);
c = tvlqr(p,xtraj1,utraj1,Q,R,Q);
p = p.setMass(1.2);
clsys = feedback(p,c);
[~,xcl] = clsys.simulate([0 4], [0 0]');
v.playback(xcl);

keyboard

p = p.setMass(1);
c = tvlqr(p,xtraj2,utraj2,Q,R,Q);
p = p.setMass(1.2);
clsys = feedback(p,c);
[~,xcl] = clsys.simulate([0 4], [0 0]');
v.playback(xcl);



