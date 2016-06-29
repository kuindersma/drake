function [f,df,d2f] = barrier(x,k)
  x = x;
  f = -k * log(x);
  df = -k * 1./x;
  d2f = k * 1./x.^2;
end