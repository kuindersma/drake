classdef SimpleTIQPController < DrakeSystem
  methods
  function obj = SimpleTIQPController(r,controller_data,options)
    % @param r rigid body manipulator instance
    % @param controller_data 
    % @param options structure for specifying objective weights, slack
    % bounds, etc.
    
    if nargin>3
      typecheck(options,'struct');
    else
      options = struct();
    end
    
    input_frame = r.getStateFrame;
    output_frame = r.getInputFrame();

    obj = obj@DrakeSystem(0,0,input_frame.dim,output_frame.dim,true,true);
    obj = setInputFrame(obj,input_frame);
    obj = setOutputFrame(obj,output_frame);

    obj.robot = r;
    obj.numq = getNumPositions(r);
    obj.controller_data = controller_data;
    
    if isfield(options,'solver') 
      % 0: fastqp, fallback to gurobi barrier (default)
      % 1: gurobi primal simplex with active sets
      typecheck(options.solver,'double');
      sizecheck(options.solver,1);
      assert(options.solver==0 || options.solver==1);
    else
      options.solver = 0;
    end
    obj.solver = options.solver;
    
    obj.gurobi_options.outputflag = 0; % not verbose
    if options.solver==0
      obj.gurobi_options.method = 2; % -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier
    else
      obj.gurobi_options.method = 0; % -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier
    end
    obj.gurobi_options.presolve = 0;
    % obj.gurobi_options.prepasses = 1;

    if obj.gurobi_options.method == 2
      obj.gurobi_options.bariterlimit = 20; % iteration limit
      obj.gurobi_options.barhomogeneous = 0; % 0 off, 1 on
      obj.gurobi_options.barconvtol = 5e-4;
    end
  end
    
  function u=output(obj,t,~,x)
    ctrl_data = obj.controller_data.data;
      
    r = obj.robot;
            
    Q = ctrl_data.Q;
    R = ctrl_data.R;
    S = ctrl_data.S;
    x0 = ctrl_data.x0;
    u0 = ctrl_data.u0;
    B = ctrl_data.B;
    
    x_bar = x - x0;

    lb = r.umin; 
    ub = r.umax;
          
    %----------------------------------------------------------------------
    % QP cost function ----------------------------------------------------
    %
    % min_{u}: ubar*R*ubar + 2*xbar'S*B
    Hqp = R;
    fqp = x_bar'*S*B - u0'*R;

    model.Q = sparse(Hqp);
    model.A = sparse(zeros(1,r.getNumInputs));
    model.rhs = 0;
    model.sense = '=';
    model.lb = lb;
    model.ub = ub;
    model.obj = fqp;
    if obj.gurobi_options.method==2
      % see drake/algorithms/QuadraticProgram.m solveWGUROBI
      model.Q = .5*model.Q;
    end
    result = gurobi(model,obj.gurobi_options);
    u = result.x;
  end
  end

  properties (SetAccess=private)
    robot; % to be controlled
    numq;
    controller_data; % shared data handle that holds S, h, foot trajectories, etc.
    W_kdot; % angular momentum cost term weight matrix
    w_qdd; % qdd objective function weight vector
    w_grf; % scalar ground reaction force weight
    w_slack; % scalar slack var weight
    slack_limit; % maximum absolute magnitude of acceleration slack variable values
    Kp_ang; % proportunal gain for angular momentum feedback
    Kp_accel; % gain for support acceleration constraint: accel=-Kp_accel*vel
    rfoot_idx;
    lfoot_idx;
    gurobi_options = struct();
    solver=0;
    use_mex;
    mex_ptr;
    lc;
    eq_array = repmat('=',100,1); % so we can avoid using repmat in the loop
    ineq_array = repmat('<',100,1); % so we can avoid using repmat in the loop
    use_bullet;
    using_flat_terrain; % true if using DRCFlatTerrain
    jlmin;
    jlmax;
    output_qdd = false;
    body_accel_input_weights; % array of doubles, negative values signal constraints
    n_body_accel_inputs; % scalar
  end
end
