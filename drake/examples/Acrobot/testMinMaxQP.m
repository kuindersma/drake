
n=2;
d=2;

Q = eye(n);
R = eye(d);

q = 0*ones(n,1);
r = 0*ones(d,1);

% min_x max_y x'Qx + q'x + y'Ry + r'y
%
%  s.t.    x\in X,   y\in Y

% transform to

% min_x,l x'Qx + q'x + l
%
%  s.t.    x\in X
%        l >= y'Ry + r'y   for all y at the vertices of Y

ix = zeros(n,n+1); ix(:,1:n)=eye(n);
il = zeros(1,n+1); il(end)=1;

ymax = 5*ones(d,1);
ymin = 2*ones(d,1);

Ain = repmat(il,4,1);
bin = zeros(4,1);

bin(1) = (ymin'*R*ymin+r'*ymin);
bin(2) = (ymax'*R*ymax+r'*ymax);
y=[ymin(1);ymax(2)];
bin(3) = (y'*R*y+r'*y);
y=[ymax(1);ymin(2)];
bin(4) = (y'*R*y+r'*y);

xmax = ones(n,1);
xmin = -ones(n,1);

Ain_fqp = full([-Ain; -eye(n)*ix; eye(n)*ix]);
bin_fqp = [-bin; -xmin; xmax];
 
H = blkdiag(Q,0);
c = [q;1];

%[alpha,info_fqp] = fastQPmex(H,c,Ain_fqp,bin_fqp,[],[],[]);


gurobi_options.outputflag = 0; % verbose flag
gurobi_options.method = 1; % -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier

model.Q = sparse(H);
model.obj = c;
model.A = sparse(Ain);
model.rhs = bin;
model.sense = repmat('>',length(bin),1);
model.lb = [xmin;-inf];
model.ub = [xmax; inf];
result = gurobi(model,gurobi_options);
alpha = result.x;

active_set = find(abs(Ain_fqp*alpha - bin_fqp)<1e-6);




% 
% xx = linspace(-1,1,100);
% yy = linspace(-1,1,100);
% 
% figure(1);
% plot(xx, xx.'*Q*xx + q'.*xx,'b-');
