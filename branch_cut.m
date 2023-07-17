function I=branch_cut(delta,k,beta,sigma_0,c)

    % Branch cut contribution, integrand.m contains the function to
    % integrate

    int=@(x) integrand(x,delta,k,beta,sigma_0,c);
    I=sigma_0*quadgk(int,1,50)/(beta*gamma(c));
    
    % The integral is normally to infinity but numerical errors appear if
    % we increase the upper boundary. However, integrand decreases really
    % fast so good approximation

end