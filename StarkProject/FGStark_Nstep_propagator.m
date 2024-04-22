function [ X_hist, accuracy_hist, status_hist, acc_hist ] = FGStark_Nstep_propagator( X0, acc, dtau, N, order, ... 
        cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
        accuracyFlag, Xf, accuracyOut, statusFlag ) 

    % dealing with acceleration  
    acc_norm = norm(acc) ; 
    q = 10^(-3) ; 

    Xk   = X0 ; 
    Xkp1 = Xf ; 
    X_hist        = [ Xk ] ; 
    accuracy_hist = [ accuracyOut ] ; 
    status_hist   = [ statusFlag ] ; 
    acc_hist      = [  ] ; 

    for i = 1 : N 

        % rv to oe 
        oe_k = rv2oe( Xk(1:6), 1 ) ; 
        nu_k = oe_k(6) ; 
        r_norm = norm( Xk(1:3) ) ;

        % orient acc in same direction as velocity 
        v_vec = Xk(4:6) ; 
        acc   = q / r_norm * ( 1 + cos(nu_k) ) * v_vec ;  
        acc   = 0 * acc ; 

        % propagate one dtau step 
        [ Xkp1, accuracyOut, statusFlag ] = FGStark_oneStep_propagator( Xk, acc, dtau, order, ... 
            cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
            accuracyFlag, Xkp1, accuracyOut, statusFlag ) ;     

        % save hist 
        X_hist        = [ X_hist ; Xkp1 ] ; 
        accuracy_hist = [ accuracy_hist ; accuracyOut ] ; 
        status_hist   = [ status_hist ; statusFlag ] ; 
        acc_hist      = [ acc_hist ; norm(acc) ] ; 

        % ready for next iter 
        Xk = Xkp1 ; 

    end 

end 