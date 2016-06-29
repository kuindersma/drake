function [f,df] = smax(alpha,varargin)
  n = 0;
  d = 0;
  df = [];
  for i=1:length(varargin)
    x = varargin{i};
    n = n + x*exp(alpha*x);
    d = d + exp(alpha*x);
  end
  f = n/d;

  if nargout > 1
    for i=1:length(varargin)
      x = varargin{i};
      df = [df, exp(alpha*x)/d * (1 + alpha*(x-f))];
    end    
  end
end

