classdef DifferentialDynamicProgramming
  %DIFFERENTIALDYNAMICPROGRAMMING A class for performing DDP or iLQR
  %
  % This class assumes that there are a fixed number (N) time steps, and
  % that the trajectory is discreteized into timesteps h (N-1), state x
  % (N), and control input u (N)

  properties (SetAccess = protected)
    options % options structure
    plant   % the robot
    compute_second_derivatives % if true: DDP, false: iLQR
    final_cost_function
    running_cost_function
    visualizer
  end
  
  methods
    function obj = DifferentialDynamicProgramming(plant,compute_second_derivatives,options)
      % @param plant the robot
      % @param N the number of time samples
      % @param options 

      if nargin < 3
        options = struct();
      end

      if ~plant.isTI
        error('Drake:DifferentialDynamicProgramming:UnsupportedPlant','Only time-invariant plants are currently supported');
      end

      if ~isfield(options,'visualizer')
        options.visualizer=[];
      end
      if ~isfield(options,'enable_visualizer')
        options.enable_visualizer=false;
      end
      if ~isfield(options,'visualizer_period')
        options.visualizer_period=10;
      end
      
      obj.options = options;
      obj.plant = plant;
      obj.compute_second_derivatives = compute_second_derivatives;
    end
    
    function obj = addFinalCost(obj,final_cost_function)
      obj.final_cost_function = final_cost_function;
    end
    
    function obj = addRunningCost(obj,running_cost_function)
      obj.running_cost_function = running_cost_function;
    end

    function obj = addVisualizer(obj,visualizer)
      obj.options.visualizer = visualizer;
    end

    function obj = setVisualizerProperties(obj,enable_visualizer,visualizer_period)
      obj.options.enable_visualizer = enable_visualizer;
      obj.options.visualizer_period = visualizer_period;
    end
    
    function [xtraj,utraj] = solveTraj(obj,xtraj0,utraj0,N,h,iter)
      xtraj=xtraj0;
      utraj=utraj0;

      p=obj.plant;
      nx=size(xtraj,1);
      nu=size(utraj,1);

    	for j=1:iter
        [V,Vx,Vxx] = feval(obj.final_cost_function,xtraj(:,N));
        K = cell(1,N-1);
        k = cell(1,N-1);

        for i=N-1:-1:1
          x = xtraj(:,i); u = utraj(:,i);

          x_ind = 1:nx;
          u_ind = nx+(1:nu);
          
          [~,df,d2f] = p.update(i*p.dt,x,u);
          Fx = full(df(:,1+x_ind)); % +1 to skip time argument
          Fu = full(df(:,1+u_ind));

          [~,dg,d2g] = feval(obj.running_cost_function,x,u);
          
          gx = dg(x_ind);
          gu = dg(u_ind);
          Gxx = d2g(x_ind,x_ind);
          Gxu = d2g(x_ind,u_ind);
          Guu = d2g(u_ind,u_ind);
          Gux = d2g(u_ind,x_ind);
          
          Qx = gx + Vx*Fx;
          Qu = gu + Vx*Fu;
          Qxx = Gxx + Fx'*Vxx*Fx;
          Qux = Gux + Fu'*Vxx*Fx;
          Quu = Guu + Fu'*Vxx*Fu;

          K{i} = -Quu\Qux;
          k{i} = -Quu\Qu;

          V = V - 0.5*k{i}'*Quu*k{i};
          Vx = Qx - K{i}*Quu*k{i};
          Vxx = Qxx - K{i}'*Quu*K{i}; 
        end

        xtraj_new = xtraj;
        utraj_new = utraj;

        for i=1:N-1
          utraj_new(:,i) = utraj(:,i) + k{i} + K{i}*(xtraj_new(:,i)-xtraj(:,i));
          xtraj_new(:,i+1) = p.update(i*p.dt,xtraj_new(:,i),utraj_new(:,i));
        end

        norm(xtraj-xtraj_new)
        
        xtraj=xtraj_new;
        utraj=utraj_new;
        
        if obj.options.enable_visualizer && ~mod(j,obj.options.visualizer_period)
          xtraj_pp = obj.reconstructStateTrajectory(xtraj,0,h);
          obj.options.visualizer.playback(xtraj_pp);
        end
      end

    end
    
    function utraj = reconstructInputTrajectory(obj,utraj,t0,h)
      % default behavior is to use first order holds, but this can be
      % re-implemented by a subclass.
      N = size(utraj,2);
      ts = linspace(t0,N*h,N);

      utraj = PPTrajectory(foh(ts,utraj));
      utraj = utraj.setOutputFrame(obj.plant.getInputFrame);
    end
    
    function xtraj = reconstructStateTrajectory(obj,xtraj,t0,h)
      % default behavior is to use first order holds, but this can be
      % re-implemented by a subclass.
      N = size(xtraj,2);
      ts = linspace(t0,N*h,N);

      xtraj = PPTrajectory(foh(ts,xtraj));
      xtraj = xtraj.setOutputFrame(obj.plant.getStateFrame);
    end
  end
  
end
