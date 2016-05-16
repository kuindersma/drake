function runSwingUp()
%% runs trajectory optimization and animates open-loop playback

p = AcrobotPlant;
v = AcrobotVisualizer(p);


% [utraj,xtraj] = swingUpTrajectory(p);

[utraj,xtraj] = robustSwingUpTrajectory(p);

Q = diag([100 100 10 10]);
R = 0.01;
Qf = 2*Q;

c = tvlqr(p,xtraj,utraj,Q,R,Qf);

sys = feedback(p,c);
%      sys = cascade(utraj,p);
% v.playback(xtraj,struct('slider',true));
% keyboard
xtraj=simulate(sys,utraj.tspan,xtraj.eval(xtraj.tspan(1)));
v.playback(xtraj,struct('slider',true));

end
