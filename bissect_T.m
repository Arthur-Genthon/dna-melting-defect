function T_star=bissect_T(a,b,xi,c,saltConc_molar,tolT)

    % ---- T_star is obtained from the ------------------------     
    % ---- bisection method, where the file h.m contains ------
    % ---- the function h(x) whose root is found --------------
    
    it=0;               % set iteration counter to zero 
    fa=h(a,xi,c,saltConc_molar); 
    fb=h(b,xi,c,saltConc_molar); 
    
    if sign(fa)==sign(fb), error('Root not in bracket'); end 
    
    while (b-a)>tolT     
       it=it+1; 
       d=a+0.5*(b-a); 
       fd=h(d,xi,c,saltConc_molar); 
       if sign(fa)~=sign(fd);   % root lies in [a,d] 
          b=d; 
          fb=fd; 
       else                     % root lies in [d,b] 
          a=d; 
          fa=fd; 
       end 
    end 
    T_star=0.5*(a+b);   % solution in Kelvin
    err=0.5*abs(b-a);   % error

end