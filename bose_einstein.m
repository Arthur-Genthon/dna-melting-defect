function y=bose_einstein(t,z,c)

    % Integrand of the Bose-Einstein integral used to compute the polylog
    % in polylog_prolong.m
 
    y=t.^(c-1)./(exp(t)./z-1);
    
end