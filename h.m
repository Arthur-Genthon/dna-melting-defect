function y = h(t_kel,xi,c,saltConc_molar)

% Compute T_star such that beta(T_star) = 1+sigma(T_star)*zeta(c)
% for a homo-DNA strand 'AAA...AAA'
% The root of h is computed in bissect_T.m with the bissection method
    
    [u_hbW_aTaS, u_hbS_aTaS, mat_u_st_ACGT_aTaS] = calculate_breaking_weights(t_kel, saltConc_molar);
    beta=u_hbW_aTaS*mat_u_st_ACGT_aTaS(1,1);
    delta=mat_u_st_ACGT_aTaS(1,1);

    y=beta-1-xi*delta*zeta(c);
   
