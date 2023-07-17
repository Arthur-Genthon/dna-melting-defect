function y=integrand(x,delta,k,beta,sigma_0,c)

    % Integrand of the integral defining the branch cut contribution I
    % in branch_cut.m
    
    y=x.^(-k).*(log(x)).^(c-1)./((1-x.*((1+sigma_0*real(polylog_prolong(c,x+delta*1i)))/beta)).^2+(pi*sigma_0*x.*(log(x)).^(c-1)/(beta*gamma(c))).^2);
end