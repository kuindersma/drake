function runDDP()

p = PendulumPlantDT;
pv = PendulumVisualizer();

T=1.0;
N = T/p.dt; % num time steps;
xtraj = zeros(2,N);
utraj = zeros(1,N);

options.visualizer = pv;
options.enable_visualizer = true;

ddp = DifferentialDynamicProgramming(p,false,options);
ddp = ddp.addRunningCost(@p.cost);
ddp = ddp.addFinalCost(@p.finalCost);
[xtraj,utraj] = ddp.solveTraj(xtraj,utraj,N,p.dt,20);

xtraj = ddp.reconstructStateTrajectory(xtraj,0,p.dt);
utraj = ddp.reconstructInputTrajectory(utraj,0,p.dt);

sys = cascade(utraj,p);

traj = simulate(sys,[0,T],xtraj.eval(0));
pv.playback(traj,struct('slider',true));

