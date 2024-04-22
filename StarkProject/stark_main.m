%% routine for FGStark_oneStep_propagator in MATLAB 

% declare variables 

% X0: Initial Conditions, acc: Stark acceleration (3D), 
X0   = [ 0.2d0, 0.d0, 0.d0, 0.d0, 3.0d0, 0.d0, 0.d0 ]' ; 
acc  = 0.1 * [ 1.d-1, 1.d-1, 1.d-1 ]' ; 
% dtau = 1e-5 ; 

% get OEs 
mu  = 1 ; 
rv0 = X0(1:6) ; 
oe0 = rv2oe( rv0, mu ) ;  
a = oe0(1) ; e = oe0(2) ; 

% compute eccentric anomaly 

% Sundman transformation parameters (dt = cSun * r**alphaSun * dtau)
cSun     = 1 ; 
alphaSun = 2 ; 

% Select a case for the Sundman transformation: 
%   0: alpha = 0, tau = time
%   1: alpha = 1, tau = eccentric anomaly
%   2: alpha = 2, tau = true anomaly
%   3: alpha = 3/2, tau = intermediate anomaly
%   4: alpha = alphaSun, generic Sundman transformation: user provided alpha
caseSun = 1 ; 

% tau period 
[ tau_p, alphaSun ] = taup_caseSun( caseSun, a, e, mu, cSun ) ; 
tau_p = 10 * tau_p ; 

% dtau: delta tau for integration for N steps 
% dtau = 6.283185307179589D-1 ; 
N    = 20 ; 
dtau = tau_p / N ; 

% order: order of the TS integration 
order = 8 ; 

% Xf: Integrated vector 
Xf = zeros(7,1) ; 

% Accuracy measurements 
accuracyOut = zeros(1,4) ; 

% Indicates whether the user wants a Stark or Kepler propagation 
caseStark = true ; 

% Indicates whether the user wants to perform divergence checks (expensive)
divSafeguard = true ; 

% Indicates possible divergence of the series. 
% If statusFlag = -1, the series likely diverged
statusFlag = 0 ; 

% Indicates which accuracy measurements are to be computed 
% (1 = Hamiltonian, 2 = delR, 3 = delV, 4 = delT)
accuracyFlag = [ true, true, true, true ]' ; 


%% begin subroutine 

disp('*** Starting the F&G Stark series propagation ***')
disp('')

disp('X0: ') ; disp(X0)
disp('')

disp(['Case: ', num2str(caseSun), ' cSun: ', num2str(cSun), ' alphaSun: ', num2str(alphaSun), ' order: ', num2str(order)])
disp('')

% propagate N dtau steps 
% [ X_hist, accuracy_hist, status_hist, acc_hist ] = FGStark_Nstep_propagator( X0, acc, dtau, N, order, ... 
%         cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
%         accuracyFlag, Xf, accuracyOut, statusFlag ) ; 
[ X_hist, accuracy_hist, status_hist, acc_hist ] = FGStark_Ntaup_propagator( ... 
        X0, acc, dtau, N, order, mu, ... 
        cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
        accuracyFlag, Xf, accuracyOut, statusFlag ) ; 
    

%% display checks and plot 

oe_hist = [ ] ; E_hist = [ ] ; 
for i = 1 : size(X_hist, 1) 
    
    rv = X_hist( i, 1:6 ) ; 
    oe = rv2oe( rv, mu ) ;  
    
    e  = oe(2) ; 
    nu = oe(6) ; 
    E  = nu2E_fn( nu, e ) ; 
    
    oe_hist = [ oe_hist ; oe' ] ; 
    E_hist  = [ E_hist ; E ] ; 
    
end 

% make E_hist propagate forward 
for i = 1 : length(E_hist) - 1 
    
    if E_hist(i) > E_hist(i+1) 
        E_hist(i+1) = pi + ( pi - E_hist(i+1) ) ; 
    end 
    
end 

disp('Status: ' ) ; disp(statusFlag) ; 
disp('Xf: ') ; disp(Xf) 
disp('del Hamiltonian:') ; disp(accuracyOut(1)) ; 

figure() ; 
    subplot(3,1,1) ; 
        semilogy( E_hist, accuracy_hist(:,1), 'b' ) ; hold on ; grid on ; 
        scatter( E_hist, accuracy_hist(:,1), 'b' ) ; 
        xlabel( {'Eccentric anomaly E'; ''} ) ; 
        ylabel( 'delta H' ) ; 
        title( sprintf( 'alpha = %.2g, order = %d', alphaSun, order ) ) ; 
    subplot(3,1,2:3) ; 
        scatter3( 0,0,0, '+k' ) ; hold on ; grid on ; 
        plot3( X_hist(:,1), X_hist(:,2), X_hist(:,3), 'b' ) ; 
        scatter3( X_hist(:,1), X_hist(:,2), X_hist(:,3), 'b' )
        scatter3( X_hist(1,1), X_hist(1,2), X_hist(1,3), '^g', 'filled' ) ; 
        scatter3( X_hist(end,1), X_hist(end,2), X_hist(end,3), 'diamondr', 'filled' ) ; 
        view(0,90)
        xlabel('x') ; ylabel('y') ; zlabel('z') ; 
        title( sprintf( 'trajectory' ) ) ; 
        

%% control policy 

figure() ; 
    subplot(3,1,1) ; 
        semilogy( acc_hist, 'b' ) ; hold on ; grid on ; 
        scatter( [ 1 : length(acc_hist) ], acc_hist, 'b' ) ; 
        xlabel( {'Segment'; ''} ) ; ylabel('LU/TU^2') ; 
        title( 'Acc mag' ) ; 
    subplot(3,1,2:3) ; 
        scatter3( 0,0,0, '+k' ) ; hold on ; grid on ; 
        plot3( X_hist(:,1), X_hist(:,2), X_hist(:,3), 'b' ) ; 
        scatter3( X_hist(:,1), X_hist(:,2), X_hist(:,3), 'b' )
        scatter3( X_hist(1,1), X_hist(1,2), X_hist(1,3), '^g', 'filled' ) ; 
        scatter3( X_hist(end,1), X_hist(end,2), X_hist(end,3), 'diamondr', 'filled' ) ; 
        view(0,90)
        xlabel('x') ; ylabel('y') ; zlabel('z') ; 
        title( sprintf( 'trajectory' ) ) ; 
    sgtitle( sprintf( 'alpha = %.2g, order = %d', alphaSun, order ) ) ;         


%% subfunctions 

% "Convert true anomaly to eccentric anomaly"
function E = nu2E_fn( nu, e ) 

    % elliptic 
    if e < 1.0 
        E = 2 * atan( sqrt( (1-e)/(1+e) ) * tan(nu/2) ) ; 

    % hyperbolic 
    else 
        E = acosh( (e+cos(nu)) / (1+e*cos(nu)) ) ; 
        if nu < 0.0 || nu > pi 
            E = -E ; 
        end 
    end 

    if E < 0.0 
        E = 2*pi + E ; 
    end

%     return E        # eccentric anomaly [rad] 
end 
        






