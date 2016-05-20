function runSwingUp(dist,M)
%% runs trajectory optimization and animates open-loop playback

p = AcrobotPlant;


% [utraj,xtraj] = swingUpTrajectory(p);

[utraj,xtraj] = robustSwingUpTrajectory(p,dist,M);

% Q = diag([100 100 10 10]);
% R = 0.01;
% Qf = 2*Q;
% 
v = AcrobotVisualizer(p);
% c = tvlqr(p,xtraj,utraj,Q,R,Qf);
% sys = feedback(p,c);
% traj=simulate(sys,utraj.tspan,xtraj.eval(xtraj.tspan(1)));
% v.playback(traj,struct('slider',true));
% 
% keyboard
% 
% p2 = AcrobotPlant(disturbances); % noisy inputs
% v2 = AcrobotVisualizer(p2);
% xtraj = xtraj.setOutputFrame(p2.getStateFrame);
% utraj = utraj.setOutputFrame(p2.getInputFrame);
% c = tvlqr(p2,xtraj,utraj,Q,R,Qf);
% sys = feedback(p2,c);
% %      sys = cascade(utraj,p);
v.playback(xtraj,struct('slider',true));
% keyboard
% for i=1:10
%   traj=simulate(sys,utraj.tspan,xtraj.eval(xtraj.tspan(1)));
%   v2.playback(traj,struct('slider',true));
% 
%   keyboard
% end
% keyboard
end
