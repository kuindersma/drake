function [ y,dydx ] = smoothlin(x)

 
  xmin = 0;
  xmax = 0.02;
  
  ymin = 0.01;
  ymax = 100;
  
  y1 = ymin + (ymax - ymin)*(x./xmax);
  dy1_dx = (ymax - ymin)*(1./xmax);
  
  beta = 1.0;
  [y2,dy2_dy1] = smax(y1,ymin,beta);
  
  [y,dy_dy2] = smin(y2,ymax,beta);
 
  dydx = dy_dy2.*dy2_dy1.*dy1_dx;
  
  
  

end

