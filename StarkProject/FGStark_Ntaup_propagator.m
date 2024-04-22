function [ X_hist, accuracy_hist, status_hist, acc_hist ] = FGStark_Ntaup_propagator( ... 
        X0, acc, dtau, N, order, mu, ... 
        cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
        accuracyFlag, Xf, accuracyOut, statusFlag ) 

    % dealing with acceleration  
    q = 10^(-3) ; 

    Xk   = X0 ; 
    Xkp1 = Xf ; 
    X_hist        = [ Xk ] ; 
    accuracy_hist = [ accuracyOut ] ; 
    status_hist   = [ statusFlag ] ; 
    acc_hist      = [  ] ; 

    for i = 1 : N 
        
        oe_k = rv2oe( Xk, mu ) ; 
        a = oe_k(1) ; e = oe_k(2) ; 

        % tau_p 
        [ tau_p, alphaSun ] = taup_caseSun( caseSun, a, e, mu, cSun ) ; 

        % propagate N steps step 
        N_step = 20 ; 
        dtau   = tau_p / N_step ; 
        [ X_Nstep_hist, accuracy_Nstep_hist, status_Nstep_hist, acc_Nstep_hist ] = ... 
            FGStark_Nstep_propagator( ... 
            Xk, acc, dtau, N_step, order, ... 
            cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
            accuracyFlag, Xf, accuracyOut, statusFlag ) ;   

        % save hist 
        X_hist        = [ X_hist ; X_Nstep_hist ] ; 
        accuracy_hist = [ accuracy_hist ; accuracy_Nstep_hist ] ; 
        status_hist   = [ status_hist ; status_Nstep_hist ] ; 
        acc_hist      = [ acc_hist ; acc_Nstep_hist ] ; 

        % ready for next iter 
        Xk = X_hist(end,:) ; 

    end 

end 