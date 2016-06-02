function [y,dydx] = smin(x,v,beta)
  
  y = -0.5*(sqrt((x-v).^2 + beta^2) - x-v);
  dydx = -0.5*(sqrt((x-v).^2 + beta^2)*(x-v) - 1);

end

