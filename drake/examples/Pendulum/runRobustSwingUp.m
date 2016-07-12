close all

%b = 0.1; %damping constant;
%popts.torque_limit = 1;
%p = PendulumPlant(b,popts);
p = PendulumPlant();
v = PendulumVisualizer();

N = 41;
Q = [10 0; 0 1];
R = .1;
Qf = 100*eye(2);

options.integration_method = DirtranTrajectoryOptimization.MIDPOINT;
[utraj1,xtraj1,z1,prog1] = p.swingUpTrajectory(N,options);

D = (2/.2^2); %This corresponds to +/-.2 uncertainty in mass (20%)
[utraj2,xtraj2,z2,prog2] = p.robustSwingUpTrajectory(N,D,Q,R,Qf);

%closed-loop
p = p.setMass(1);
c = tvlqr(p,xtraj1,utraj1,Q,R,Qf);
p = p.setMass(1.2);
p.limit_torque = 1;
clsys = feedback(p,c);
[~,xcl1] = clsys.simulate([0 4], [0 0]');
v.playback(xcl1);

p = p.setMass(1);
c = tvlqr(p,xtraj2,utraj2,Q,R,Qf);
p = p.setMass(1.2);
p.limit_torque = 1;
clsys = feedback(p,c);
[~,xcl2] = clsys.simulate([0 4], [0 0]');
v.playback(xcl2);
v.playbackSWF(xcl2, 'swing2.swf');

%Write movie files
setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
setenv('PATH', [getenv('PATH') ':/Library/TeX/texbin']);
v.playbackSWF(xcl1, 'swing1.swf');
v.playbackSWF(xcl2, 'swing2.swf');

%Plots
h1 = z1(prog1.h_inds);
t1 = [0; cumsum(h1)];
x1 = z1(prog1.x_inds);
u1 = z1(prog1.u_inds);

h2 = z2(prog2.h_inds);
t2 = [0; cumsum(h2)];
x2 = z2(prog2.x_inds);
dx = z2(prog2.dx_inds);
u2 = z2(prog2.u_inds);
du = z2(prog2.du_inds);
w = z2(prog2.w_inds);

figure(2);
subplot(2,1,1);
plot(t1(1:end-1),u1(1:end-1));
ylabel('u_{dirtran}');
xlim([0 t2(end)]);
ylim([-3.5 3.5]);
subplot(2,1,2);
plot(t2(1:end-1),u2(1:end-1));
xlim([0 t2(end)]);
ylim([-3.5 3.5]);
ylabel('u_{robust}');

figure(3);
subplot(2,1,1);
plot(t1,x1(1,:));
hold on
plot(t1,x1(2,:));
ylabel('x_{dirtran}');
xlim([0 t2(end)]);
l = legend('$\theta$', '$\dot{\theta}$');
set(l,'Interpreter','latex')
subplot(2,1,2);
plot(t2,x2(1,:));
hold on
plot(t2,x2(2,:));
ylabel('x_{robust}');
xlim([0 t2(end)]);

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
ylabel('w');
