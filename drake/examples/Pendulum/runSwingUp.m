function runSwingUp()
% runs trajectory optimization and animates 

pd = PendulumPlant;
pv = PendulumVisualizer();
[utraj,xtraj] = pd.swingUpTrajectory();

if (0) % open-loop playback
  sys = cascade(utraj,pd);
  xtraj=simulate(sys,utraj.tspan,[0;0]);
end
pv.playback(xtraj);

end
