p = PendulumPlant();
v = PendulumVisualizer();

N = 21;
D = 200; %This corresponds to +/-0.1 uncertainty in mass

%Turn on if statement on line 67 of DirectTrajectoryOptimization.m
%keyboard 

% [utraj1,xtraj1] = p.swingUpTrajectory(N);
% 
% %Plot nominal trajectory
% t1 = utraj1.tspan(1):.1:utraj1.tspan(2);
% uhist = utraj1.eval(t1);
% xhist = xtraj1.eval(t1);
% figure(2);
% plot(t1,uhist);
% figure(3);
% plot(t1,xhist(1,:));
% hold on
% plot(t1,xhist(2,:));

%Turn off if statement on line 67 of DirectTrajectoryOptimization.m
keyboard 

[~,~,z,prog] = p.robustSwingUpTrajectory(N,D);

%Plot robust trajectory
h = z(prog.h_inds);
t2 = [0; cumsum(h)];
x = z(prog.x_inds);
dx = z(prog.dx_inds);
u = z(prog.u_inds);
du = z(prog.du_inds);
w = z(prog.w_inds);
figure(2);
hold on
plot(t2,u);
figure(3);
hold on
plot(t2,x(1,:));
plot(t2,x(2,:));
figure(4);
plot(t2,dx(1,:));
hold on;
plot(t2,dx(2,:));
figure(5);
plot(t2,du);
figure(6);
plot(t2,w);
