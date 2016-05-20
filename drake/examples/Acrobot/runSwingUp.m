function runSwingUp(N,M,dist)
%% runs trajectory optimization and animates open-loop playback

p = AcrobotPlant;

% [utraj,xtraj] = swingUpTrajectory(p,N);

[utraj,xtraj] = robustSwingUpTrajectory(p,N,M,dist);
% [utraj,xtraj] = robustSwingUpTrajectory(p,N,1,0);

Q = diag([100 100 10 10]);
R = 0.01;
Qf = 2*Q;

% show plan
v = AcrobotVisualizer(p);
v.playback(xtraj,struct('slider',true));
keyboard;

% % run LQR on nominal plan
% c = tvlqr(p,xtraj,utraj,Q,R,Qf);
% sys = feedback(p,c);
% traj=simulate(sys,utraj.tspan,xtraj.eval(xtraj.tspan(1)));
% v.playback(traj,struct('slider',true));
% 
% keyboard;
 
p2 = AcrobotPlant(); % noisy inputs
v2 = AcrobotVisualizer(p2);
xtraj = xtraj.setOutputFrame(p2.getStateFrame);
utraj = utraj.setOutputFrame(p2.getInputFrame);
c = tvlqr(p2,xtraj,utraj,Q,R,Qf);
p2 = p2.addDisturbanceBound(dist);
sys = feedback(p2,c);
for i=1:10
  traj=simulate(sys,utraj.tspan,xtraj.eval(xtraj.tspan(1)));
  v2.playback(traj,struct('slider',true));

  keyboard
end

keyboard
end
