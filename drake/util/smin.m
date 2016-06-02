function [y,dydx] = smin(x,v,beta)
  
  sq = sqrt((x-v).^2 + beta^2);

  y = -0.5*(sq-x-v);
  dydx = -0.5*(1.0/sq*(x-v)-1);

end

