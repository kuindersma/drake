function [f,df,d2f] = smooth_barrier(x,k,q)

  x0 = k/q;

  i_ge = x>=x0;
  [fg,dfg,d2fg] = barrier(x(i_ge),k);
  
  i_le = ~i_ge;
  a0 = 0.5*k * 1./x0.^2;
  a1 = -2*k * 1./x0;
  a2 = -k * (log(x0) + 0.5 - 2);
 
  fl = a0*x(i_le).^2 + a1*x(i_le) + a2;
  dfl = 2*a0*x(i_le) + a1;
  d2fl = 2*a0;
  
  f = 0*x;
  f(i_ge) = fg;
  f(i_le) = fl;

  df = 0*x;
  df(i_ge) = dfg;
  df(i_le) = dfl;
  
  d2f = 0*x;
  d2f(i_ge) = d2fg;
  d2f(i_le) = d2fl;
end