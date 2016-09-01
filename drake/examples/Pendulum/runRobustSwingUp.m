close all

%b = 0.1; %damping constant;
%popts.torque_limit = 1;
%p = PendulumPlant(b,popts);
p = PendulumPlant();
v = PendulumVisualizer();

N = 31;

options.integration_method = DirtranTrajectoryOptimization.MIDPOINT;
[utraj1,xtraj1,z1,prog1] = p.swingUpTrajectory(N,options);

D = (1/.1^2); %This corresponds to +/-.1 uncertainty in mass (10%)
[utraj2,xtraj2,z2,prog2] = p.robustSwingUpTrajectory(N,D,options);

%closed-loop
Q = [100 0; 0 10];
R = 1;
Qf = 1000*eye(2);
      
%p = p.setMass(1);
c = tvlqr(p,xtraj1,utraj1,Q,R,Qf);
%p = p.setMass(1.2);
p.limit_torque = 1;
clsys = feedback(p,c);
[~,xcl1] = clsys.simulate([0 3], [0 0]');
v.playback(xcl1);

%p = p.setMass(1);
c = tvlqr(p,xtraj2,utraj2,Q,R,Qf);
%p = p.setMass(1.2);
p.limit_torque = 1;
clsys = feedback(p,c);
[~,xcl2] = clsys.simulate([0 4], [0 0]');
v.playback(xcl2);

%Write movie files
%v.playbackAVI(xcl1, 'swing1.avi');
%v.playbackAVI(xcl2, 'swing2.avi');
% setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
% setenv('PATH', [getenv('PATH') ':/Library/TeX/texbin']);
% v.playbackSWF(xcl1, 'swing1.swf');
% v.playbackSWF(xcl2, 'swing2.swf');

%Plots
h1 = z1(prog1.h_inds);
t1 = [0; cumsum(h1)];
x1 = z1(prog1.x_inds);
u1 = z1(prog1.u_inds);

h2 = z2(prog2.h_inds);
t2 = [0; cumsum(h2)];
x2 = z2(prog2.x_inds);
u2 = z2(prog2.u_inds);

figure(2);
subplot(2,1,1);
plot(t1(1:end-1),u1);
ylabel('u_{dirtran}');
xlim([0 t2(end)]);
ylim([-3.5 3.5]);
subplot(2,1,2);
plot(t2(1:end-1),u2);
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
