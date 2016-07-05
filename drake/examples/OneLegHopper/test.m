% 
%   R_min = 0.1;
%   R_max = 20;
%   k = 10;
% 
%   y = randn();
%   
%   df = 1./(1+exp(-k*y)).^2 .* R_max*k .* exp(-k*y);
% 
%   
%   [f1,df1] = geval(@logistic,y,struct('grad_method','taylorvar'));
%   
%   valuecheck(df,df1);
%   function r = logistic(y)
%       r = R_min + R_max./(1+exp(-k*y));
%   end
% 

  k = 1;
  xmin = 0;
  xmax = 0.02;
  x = linspace(xmin-0.01,xmax+.1,200)';

  ymin = 0.001;
  ymax = 100;
  [y,dy] = sig(x,k,xmin,xmax,ymin,ymax);

  plot(x,-y)

  tmp = @(x) sig(x,k,xmin,xmax,ymin,ymax);
  [y1,dy1] = geval(tmp,x,struct('grad_method','numerical'));
  
%   valuecheck(dy,diag(dy1),1e-3);


  
  figure (2)
  
  [y,dydx] = smoothlin(x);
  [y1,dydx1] = geval(@smoothlin,x,struct('grad_method','taylorvar'));
  

  
  plot(x,y) 
  
  