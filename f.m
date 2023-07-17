function y = f(x,sigma0,betap,c,tolp)   

    % Function whose root defines z_0, found using the bissection
    % method in bissect_z.m

    temp=polylogT(c,x,tolp); 
    y=sigma0*x*temp-(betap-x);
   
