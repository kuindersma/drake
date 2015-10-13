function runDDP()

p = PendulumPlantDT;
pv = PendulumVisualizer();

T = 3.0;
N = T/p.dt; % num time steps;
xtraj = zeros(2,N);
utraj = zeros(1,N);

options.visualizer = pv;
options.enable_visualizer = true;

ddp = DifferentialDynamicProgramming(p,false,options);
ddp = ddp.addRunningCost(@p.cost);
ddp = ddp.addFinalCost(@p.finalCost);
[xtraj,utraj] = ddp.solveTraj(xtraj,utraj,N,p.dt,30);

xtraj = ddp.reconstructStateTrajectory(xtraj,0,p.dt);
utraj = ddp.reconstructInputTrajectory(utraj,0,p.dt);


if (0) % open-loop playback
  sys = cascade(utraj,p);
  xtraj = simulate(sys,[0,T],xtraj.eval(0));
end
pv.playback(xtraj,struct('slider',true));

