function [ X_hist, accuracy_hist, status_hist ] = FGStark_Nstep_propagator( X0, acc, dtau, N, order, ... 
        cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
        accuracyFlag, Xf, accuracyOut, statusFlag ) 

Xk   = X0 ; 
Xkp1 = Xf ; 
X_hist        = [ Xk' ] ; 
accuracy_hist = [ accuracyOut ] ; 
status_hist   = [ statusFlag ] ; 

for i = 1 : N 
    
    % propagate one dtau step 
    [ Xkp1, accuracyOut, statusFlag ] = FGStark_oneStep_propagator( Xk, acc, dtau, order, ... 
        cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
        accuracyFlag, Xkp1, accuracyOut, statusFlag ) ;     
    
    % save hist 
    X_hist        = [ X_hist ; Xkp1' ] ; 
    accuracy_hist = [ accuracy_hist ; accuracyOut ] ; 
    status_hist   = [ status_hist ; statusFlag ] ; 
    
    % ready for next iter 
    Xk = Xkp1 ; 
    
end 

end 