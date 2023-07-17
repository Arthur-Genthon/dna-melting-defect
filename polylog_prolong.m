function y=polylog_prolong(c,z)  

    % Bose-Einstein representation of the polylogarithm is used for
    % the branch cut contribution
     
    if z==1                                                       
         
        y=zeta(c);
        
     elseif real(z)>1 & imag(z)==0
         
        error('Error. z is in the branch cut.')

     else
         
         be=@(t) bose_einstein(t,z,c);                         
         y=(1/gamma(c))*integral(be,0,Inf,'ArrayValued',true);   

        % Arrayvalued is necessary, see:
        % https://se.mathworks.com/matlabcentral/answers/103319-error-in-numerical-integration-matrix-dimensions-must-agree
    end
end
