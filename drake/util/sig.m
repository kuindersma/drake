    
    function [y,dydx] = sig(x,k,xmin,xmax,ymin,ymax)
      
      xscaled = 2*(x-xmin)./ (xmax-xmin) - 1; % scale to between -1 and 1
      y = ymin + ymax./(1+exp(-k*xscaled));
      dydx = 1./(1+exp(-k*xscaled)).^2 .* ymax*k .* exp(-k*xscaled) * 2./ (xmax-xmin);
    
    end
    
    
    
    
    
    
    
    