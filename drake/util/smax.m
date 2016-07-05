function [y,dydx] = smax(x,v,beta)
  
  [y,dydx] = smin(-x,-v,beta);
  y = -y;
  dydx = dydx;

end

