function runSwingUp

p = CartDoublePendulumPlant();

[utraj,xtraj]=swingUpTrajectory(p);

v = CartDoublePendulumVisualizer(p);
v.playback(xtraj);

