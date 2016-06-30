close all

%b = 0.1; %damping constant;
%popts.torque_limit = 1;
%p = PendulumPlant(b,popts);
p = PendulumPlant();
v = PendulumVisualizer();

N = 21;
D = 200; %This corresponds to +/-0.1 uncertainty in mass

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
% x1_e = zeros(size(x1));
% x1_e(:,1) = x1(:,1);
% for k = 1:(length(t1)-1)
%     x1_e(:,k+1) = x1_e(:,k) + h1(k)*p.dynamics(0,x1_e(:,k),u1(k));
% end
% figure();
% plot(t1,x1_e(1,:));
% hold on;
% plot(t1,x1_e(2,:));
% 
% x2_e = zeros(size(x2));
% x2_e(:,1) = x2(:,1);
% for k = 1:(length(t2)-1)
%     x2_e(:,k+1) = x2_e(:,k) + h2(k)*p.dynamics(0,x2_e(:,k),u2(k));
% end
% figure();
% plot(t2,x2_e(1,:));
% hold on;
% plot(t2,x2_e(2,:));

% xtraj1_ol = PPTrajectory(foh(t1,x1_e));
% xtraj1_ol = xtraj1_ol.setOutputFrame(p.getStateFrame);
v.playback(xtraj1);

% xtraj2_ol = PPTrajectory(foh(t2,x2_e));
% xtraj2_ol = xtraj2_ol.setOutputFrame(p.getStateFrame);
v.playback(xtraj2);




