close all

%b = 0.1; %damping constant;
%popts.torque_limit = 1;
%p = PendulumPlant(b,popts);
p = PendulumPlant();
v = PendulumVisualizer();

N = 51;
Q = [10 0; 0 1];
R = .1;
Qf = 100*eye(2);

options.integration_method = DirtranTrajectoryOptimization.MIDPOINT;
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

D = (2/.2^2); %This corresponds to +/-.2 uncertainty in mass (20%)
[utraj2,xtraj2,z2,prog2] = p.robustSwingUpTrajectory(N,D,Q,R,Qf);

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

%Get linearized dynamics along trajectory
% Ak = zeros(2,2,N-1);
% Bk = zeros(2,1,N-1);
% for k = 1:(N-1)
%     [~,dx1] = prog2.forward_robust_dynamics_fun(h2(k),x2(:,k),u2(k),0);
%     Ak(:,:,k) = dx1(:,2:3);
%     Bk(:,:,k) = dx1(:,4);
% end

%Calculate LQR gains
% S = Qf;
% K = zeros(1,2,N);
% for k = (N-1):-1:1
%     K(:,:,k) = (Bk(:,:,k)'*S*Bk(:,:,k)+R)\(Bk(:,:,k)'*S*Ak(:,:,k));
%     S = Q + K(:,:,k)'*R*K(:,:,k) + (Ak(:,:,k) - Bk(:,:,k)*K(:,:,k))'*S*(Ak(:,:,k) - Bk(:,:,k)*K(:,:,k));
% end

% %open-loop playback
% xsim = zeros(2,N);
% xsim(:,1) = x2(:,1);
% for k = 1:(N-1)
%     xsim(:,k+1) = xsim(:,k) + h2(k)*p.dynamics_w(0,xsim(:,k),u2(k)-K(:,:,k)*(xsim(:,k)-x2(:,k)),.4);
% end
% xol = PPTrajectory(foh(t2,xsim));
% xol = xol.setOutputFrame(p.getStateFrame);
% v.playback(xol);

%closed-loop
p = p.setMass(1);
c = tvlqr(p,xtraj1,utraj1,Q,R,Qf);
p = p.setMass(1.2);
p.limit_torque = 1;
clsys = feedback(p,c);
[~,xcl] = clsys.simulate(utraj1.tspan, [0 0]');
v.playback(xcl);

p = p.setMass(1);
c = tvlqr(p,xtraj2,utraj2,Q,R,Qf);
p = p.setMass(1.2);
p.limit_torque = 1;
clsys = feedback(p,c);
[~,xcl] = clsys.simulate(utraj2.tspan, [0 0]');
v.playback(xcl);



