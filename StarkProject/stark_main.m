%% routine for FGStark_oneStep_propagator in MATLAB 

% declare variables 

% X0: Initial Conditions, acc: Stark acceleration (3D), dtau: delta tau for integration
X0   = zeros(7,1) ; 
acc  = zeros(3,1) ; 
dtau = 0 ; 

% Sundman transformation parameters (dt = cSun * r**alphaSun * dtau)
cSun     = 0 ; 
alphaSun = 0 ; 

% Select a case for the Sundman transformation: 0: alpha = 0, tau = time
%   1: alpha = 1, tau = eccentric anomaly
%   2: alpha = 2, tau = true anomaly
%   3: alpha = 3/2, tau = intermediate anomaly
%   4: alpha = alphaSun, generic Sundman transformation: user provided alpha
caseSun = 0 ; 

% order: order of the TS integration 
order = 0 ; 

% Xf: Integrated vector 
Xf = zeros(7,1) ; 

% Accuracy measurements 
accuracyOut = zeros(4,1) ; 

% Indicates whether the user wants a Stark or Kepler propagation 
caseStark = 0 ; 

% Indicates whether the user wants to perform divergence checks (expensive)
divSafeguard = 0 ; 

% Indicates possible divergence of the series. 
% If statusFlag = -1, the series likely diverged
statusFlag = 0 ; 

% Indicates which accuracy measurements are to be computed 
% (1 = Hamiltonian, 2 = delR, 3 = delV, 4 = delT)
accuracyFlag = logical(zeros(4,1)) ; 


%% begin subroutine 

X0   = [ 0.2d0, 0.d0, 0.d0, 0.d0, 3.0d0, 0.d0, 0.d0 ]' ; 

acc  = [ 1.d-1, 1.d-1, 1.d-1 ]' ; 
dtau = 6.283185307179589D-1 ; 

cSun     = 1.0d0 ; 
caseSun  = 0 ; 
alphaSun = 2.d0 ; 

caseStark    = true ; 
divSafeguard = true ; 

order = 12 ; 
accuracyFlag = [ true, true, true, true ]' ; 

disp('*** Starting the F&G Stark series propagation ***')
disp('')

disp('X0: ') ; disp(X0)
disp('')

disp(['Case: ', num2str(caseSun), ' cSun: ', num2str(cSun), ' alphaSun: ', num2str(alphaSun), ' order: ', num2str(order)])
disp('')

[ Xf, accuracyOut, statusFlag ] = FGStark_oneStep_propagator( X0, acc, dtau, order, ... 
    cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
    accuracyFlag, Xf, accuracyOut, statusFlag ) ; 

disp('Xf: ') ; disp(Xf) 

%% subfunctions 





