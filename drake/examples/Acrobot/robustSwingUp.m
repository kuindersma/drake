function robustSwingUp(N,M,lb,ub)
%% runs trajectory optimization and animates open-loop playback

p = AcrobotPlant;

K = 0.2*[-100,-100,-30,-30]; % PD gains

if 1
 [utraj,xtraj,z,prog] = robustSwingUpTrajectory(p,N,M,lb,ub,K);
%     xs = z(prog.x_inds);
%     dxs = z(prog.dx_inds);
%     us = z(prog.u_inds);
%     dus = z(prog.du_inds);
%     hs = z(prog.h_inds);
% 
% %     xs_ = xs;
% %     xs__ = xs;
% %     for i=1:N-1
% %       xs_(:,i+1) = xs_(:,i) + hs(i) * p.dynamics_w(0,xs_(:,i), us(i) + dus(i), lb);
% %       xs__(:,i+1) = xs__(:,i) + hs(i) * p.dynamics_w(0,xs__(:,i), us(i) + K*(xs__(:,i)-xs(:,i)), lb);
% %     end
% 
%     figure(12);
%     plot(xs(1,:),xs(2,:),'b.-','MarkerSize',10);
%     hold on;
%     plot(xs(1,:)+dxs(1,:),xs(2,:)+dxs(2,:),'g.-','MarkerSize',10);
% %     plot(xs_(1,:),xs_(2,:),'r.-','MarkerSize',10);
% %       plot(xs__(1,:),xs__(2,:),'m.-','MarkerSize',10);
%     hold off;
% 
%     keyboard
elseif 1
 
  [utraj,xtraj] = swingUpTrajectory(p,N);


else
  load nominal_20.mat

  xtraj = xtraj.setOutputFrame(p.getStateFrame);
  utraj = utraj.setOutputFrame(p.getInputFrame);
end



% [utraj,xtraj] = robustSwingUpTrajectory(p,N,M,dist);
% [utraj,xtraj] = robustSwingUpTrajectory(p,N,1,0);

% Q = diag([100 100 10 10]);
% R = 0.01;
% Qf = 2*Q;

% show plan
v = AcrobotVisualizer(p);
% v.playback(xtraj,struct('slider',true));
% keyboard;

% % run LQR on nominal plan
% c = tvlqr(p,xtraj,utraj,Q,R,Qf);
% sys = feedback(p,c);
% traj=simulate(sys,utraj.tspan,xtraj.eval(xtraj.tspan(1)));
% v.playback(traj,struct('slider',true));
% 
% keyboard;
 
% c = tvlqr(p,xtraj,utraj,Q,R,Qf);

% use fixed linear PD
nX = p.getNumStates;
nU = p.getNumInputs;
% create linear feedback controller
iframe = CoordinateFrame([p.getStateFrame.name,' - x0(t)'],nX,p.getStateFrame.prefix);
p.getStateFrame.addTransform(AffineTransform(p.getStateFrame,iframe,eye(nX),-xtraj));
iframe.addTransform(AffineTransform(iframe,p.getStateFrame,eye(nX),xtraj));
oframe = CoordinateFrame([p.getInputFrame.name,' + u0(t)'],nU,p.getInputFrame.prefix);
oframe.addTransform(AffineTransform(oframe,p.getInputFrame,eye(nU),utraj));
p.getInputFrame.addTransform(AffineTransform(p.getInputFrame,oframe,eye(nU),-utraj));
D = ConstantTrajectory(K);
c = AffineSystem([],[],[],[],[],[],[],D,[]);
c = setInputFrame(c,iframe);
c = setOutputFrame(c,oframe);

sys = feedback(p,c);
for i=1:1
  traj=simulate(sys,xtraj.tspan,xtraj.eval(xtraj.tspan(1)));
  v.playback(traj,struct('slider',true));
end

keyboard
end
