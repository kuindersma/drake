function [x, lambda] = qcqp(A,b,D,xguess)
%QCQP solves the following problem:
% min (1/2)*x'*A*x + b'*x
% s.t. (1/2)*x'*D*x <= 1
% where D > 0 defines an ellipsoidal region and A can be indefinite.
% The algorithm is a primal-dual interior point method and uses several
% heuristics to attempt to handle non-convex cases.

%Algorithm parameters
mu = 10;
alpha = 0.1;
beta = 0.6;
tol = 1e-12;

if nargin == 4 && ~isempty(xguess)
    x = xguess;
    lambda = -mean((A*x+b)./(D*x));
    
    %check residual
    if ~any((A+lambda*D)*x + b > tol)
        return
    end
else
    %Get feasible initial lambda using Gershgorin Circle Theorem
    [Amin, ind] = min(diag(A));
    lambda = max(1, 1.1*sum(abs(A(ind,:))));
    
    %Heuristic initial guess for x. If it looks like the minimum eigenvalue of
    %A is negative, start near the boundary. Otherwise, start at the origin.
    x = zeros(size(b));
    if Amin < 0
        if b(ind) ~= 0
            x(ind) = -sign(b(ind))/sqrt(max(diag(D)));
        else
            x(ind) = 1/sqrt(max(diag(D)));
        end
    end
end

r = ones(length(b)+1,1);
eta = 1;
while any(abs(r) > tol) && eta > tol
    
    eta = -lambda*(.5*x'*D*x - 1); %surrogate duality gap
    t = mu/eta; %barrier parameter

    r = [(A+lambda*D)*x + b; -lambda*(.5*x'*D*x-1) - 1/t]; %residual
    dr = [A + lambda*D, D*x; -lambda*x'*D, -(.5*x'*D*x-1)];
    delta = -dr\r; %descent direction
    dx = delta(1:end-1);
    dl = delta(end);

    %line search
    if dl < 0
        smax = min(1, -lambda/dl);
    else
        smax = 1;
    end
    s = .99*smax;
    xnew = x + s*dx;
    while 0.5*xnew'*D*xnew > 1
        s = beta*s;
        xnew = x + s*dx;
    end

    lnew = lambda + s*dl;
    rnew = [A*xnew + b + lnew*D*xnew; -lnew*(.5*xnew'*D*xnew-1) - 1/t];
    while sqrt(rnew'*rnew) > (1-alpha*s)*sqrt(r'*r)
        s = beta*s;
        lnew = lambda + s*dl;
        xnew = x + s*dx;
        rnew = [A*xnew + lnew*D*xnew; -lnew*(.5*xnew'*D*xnew-1) - 1/t];
    end

    x = xnew;
    lambda = lnew;

end

end
