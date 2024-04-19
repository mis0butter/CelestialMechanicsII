%% routine for FGStark_oneStep_propagator in MATLAB 

% declare variables 

% X0: Initial Conditions, acc: Stark acceleration (3D), 
X0   = [ 0.2d0, 0.d0, 0.d0, 0.d0, 3.0d0, 0.d0, 0.d0 ]' ; 
acc  = [ 1.d-1, 1.d-1, 1.d-1 ]' ; 
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
caseSun = 2 ; 

% tau period 
n = sqrt( mu / a^3 ) ; 

% 0: alpha = 0, tau = time 
if caseSun == 0 
    tau_p = 2 * pi / ( n * cSun ) ; 
    alphaSun = 0 ; 

% 1: alpha = 1, tau = eccentric anomaly 
elseif caseSun == 1 
    tau_p = 2 * pi / ( n * cSun * a ) ;  
    alphaSun = 1 ; 

% 2: alpha = 2, tau = true anomaly 
elseif caseSun == 2 
    tau_p = 2 * pi / ( n * cSun * sqrt( a * ( 1 - e^2 ) ) ) ; 
    alphaSun = 2 ; 

% 3: alpha = 3/2, tau = intermediate anomaly 
elseif caseSun == 3 
    M     = sqrt( 2 * e / ( 1 + e ) ) ; 
    tau_p = 4 * K( M ) / ( cSun * sqrt( mu * ( 1 + e ) ) ) ; 
    alphaSun = 3/2 ; 

% 4: alpha = alphaSun, generic Sundman transformation: user provided alpha  
else 
    tau_p = 2 * pi ; 
    
end 

% dtau: delta tau for integration for N steps 
% dtau = 6.283185307179589D-1 ; 
N    = 50 ; 
dtau = tau_p / N ; 

% order: order of the TS integration 
order = 8 ; 

% Xf: Integrated vector 
Xf = zeros(7,1) ; 

% Accuracy measurements 
accuracyOut = zeros(4,1) ; 

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

% propagate one dtau step 
[ Xf, accuracyOut, statusFlag ] = FGStark_oneStep_propagator( X0, acc, dtau, order, ... 
    cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
    accuracyFlag, Xf, accuracyOut, statusFlag ) ; 

% propagate N dtau steps 
Xk   = X0 ; 
Xkp1 = Xf ; 
X_hist = [ Xk' ] ; 
accuracy_hist = [ accuracyOut ] ; 
for i = 1 : N 
    
    % propagate one dtau step 
    [ Xkp1, accuracyOut, statusFlag ] = FGStark_oneStep_propagator( Xk, acc, dtau, order, ... 
        cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
        accuracyFlag, Xkp1, accuracyOut, statusFlag ) ;     
    
    % save hist 
    X_hist = [ X_hist ; Xkp1' ] ; 
    accuracy_hist = [ accuracy_hist ; accuracyOut ] ; 
    
    % ready for next iter 
    Xk = Xkp1 ; 
    
end 

%% display checks and plot 

disp('Status: ' ) ; disp(statusFlag) ; 
disp('Xf: ') ; disp(Xf) 
disp('del Hamiltonian:') ; disp(accuracyOut(1)) ; 

figure() ; hold on ; grid on ; 
    semilogy( accuracy_hist(:,1) ) ; 
    title( sprintf( ' alpha = %d', alphaSun ) ) ; 

        






