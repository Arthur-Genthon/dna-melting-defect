function z0=bissect_z(a,b,bb,sb,c,tolp)

    % ---- The root, 0<=x0<=1, is obtained from the      ----
    % ---- bisection method, where the file f.m contains ----
    % ---- the function f(x) to which the root is found  ----
    
    it=0;               % set iteration counter to zero 
    fa=f(a,sb,bb,c,tolp); 
    fb=f(b,sb,bb,c,tolp); 

    if sign(fa)==sign(fb), error('Root not in bracket'); end 

    while (b-a)>tolp     
       it=it+1; 
       d=a+0.5*(b-a); 
       fd=f(d,sb,bb,c,tolp); 
       if sign(fa)~=sign(fd);   % root lies in [a,d] 
          b=d; 
          fb=fd; 
       else                     % root lies in [d,b] 
          a=d; 
          fa=fd; 
       end 
    end 
    z0=0.5*(a+b);       % solution 
    err=0.5*abs(b-a);   % error

end