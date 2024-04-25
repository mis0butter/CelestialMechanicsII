%% routine for FGStark_oneStep_propagator in MATLAB 

% declare variables 

% X0: Initial Conditions, acc: Stark acceleration (3D), 
X0   = [ 0.2d0, 0.d0, 0.d0, 0.d0, 3.0d0, 0.d0, 0.d0 ] ; 
acc  = 0.1 * [ 1.d-1, 1.d-1, 1.d-1 ] ; 
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
% tau_p = 10 * tau_p ; 

% dtau: delta tau for integration for N steps 
% dtau = 6.283185307179589D-1 ; 
N    = 20 ; 
dtau = tau_p / N ; 

% order: order of the TS integration 
order = 8 ; 

% Xf: Integrated vector 
Xf = zeros(1,7) ; 

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
accuracyFlag = [ true, true, true, true ] ; 

%% acc hist in same direction as orbit velocity 

acc_hist = [ ... 
                         0                      0.03                         0
      -0.00881282482074017        0.0162890426326284                         0
       -0.0048435006051358       0.00401130451122418                         0
      -0.00173045139806898      0.000760741892234286                         0
     -0.000586534249761895      0.000117725374996398                         0
     -0.000202307340871176      1.81654153516559e-06                         0
     -6.93678326527898e-05     -1.25178774311193e-05                         0
     -2.20031301561146e-05     -9.02660936666008e-06                         0
     -5.54933306140563e-06     -4.27240674054442e-06                         0
      -7.2131430275336e-07      -1.1821146497067e-06                         0
     -9.61283771440908e-10     -1.51085172675286e-08                         0
      3.30057953060062e-07     -7.00575422748803e-07                         0
      3.63916385209905e-06     -3.27005164495005e-06                         0
      1.59963602396139e-05     -7.60845810961064e-06                         0
       5.2369996003557e-05      -1.1845759104201e-05                         0
      0.000154294458507784     -4.86414220223277e-06                         0
      0.000445849217276345      6.96315415458615e-05                         0
       0.00131059271017454      0.000495300179178469                         0
       0.00377856447509965       0.00269236428476427                         0
        0.0082199439820226        0.0120298333835178                         0 ... 
        ] ; 

%% begin subroutine 

disp('*** Starting the F&G Stark series propagation ***')
disp('')

disp('X0: ') ; disp(X0)
disp('')

disp(['Case: ', num2str(caseSun), ' cSun: ', num2str(cSun), ' alphaSun: ', num2str(alphaSun), ' order: ', num2str(order)])
disp('')

% propagate N dtau steps 
[ X_hist, accuracy_hist, status_hist, acc_hist ] = FGStark_Nstep_ctrl_prop( ... 
        X0, acc_hist, dtau, N, order, ... 
        cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
        accuracyFlag, Xf, accuracyOut, statusFlag ) ; 
% [ X_hist, accuracy_hist, status_hist, acc_hist ] = FGStark_Nstep_propagator( X0, acc, dtau, N, order, ... 
%         cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
%         accuracyFlag, Xf, accuracyOut, statusFlag ) ; 
% [ X_hist, accuracy_hist, status_hist, acc_hist ] = FGStark_Ntaup_propagator( ... 
%         X0, acc, dtau, N, order, mu, ... 
%         cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
%         accuracyFlag, Xf, accuracyOut, statusFlag ) ; 
    

%% display checks and plot 

disp('Status: ' ) ; disp(statusFlag) ; 
disp('Xf: ') ; disp(Xf) 
disp('del Hamiltonian:') ; disp(accuracyOut(1)) ; 

f = plot_dH_traj( X_hist, mu, accuracy_hist, alphaSun, order ) ; 


%% plot control policy 

f = plot_ctrl_traj( X_hist, acc_hist, alphaSun, order ) ; 

        
%% try fmincon 

fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2 ;  
x0  = [-1, 2] ; 

x = fmincon( fun, x0 ) 
        
        
%% try fmincon 

% final position 
rf = [ 0.5, 0, 0 ] ;  

delta_r = miss_distance( ... 
        X0, rf, acc_hist, dtau, N, order, ... 
        cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
        accuracyFlag, Xf, accuracyOut, statusFlag ) ; 

fun = @(acc_hist) miss_distance( ... 
        X0, rf, acc_hist, dtau, N, order, ... 
        cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
        accuracyFlag, Xf, accuracyOut, statusFlag ) ; 
x0  = acc_hist ; 

acc_soln = fmincon( fun, x0 ) 


%% try fmincon solution 

[ X_hist, accuracy_hist, status_hist, acc_soln ] = FGStark_Nstep_ctrl_prop( ... 
        X0, acc_soln, dtau, N, order, ... 
        cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
        accuracyFlag, Xf, accuracyOut, statusFlag ) ; 
        

%% plot control policy 

figure() ; 
    subplot(3,1,1) ; 
        semilogy( acc_hist_norm, 'b' ) ; hold on ; grid on ; 
        scatter( [ 1 : length(acc_hist) ], acc_hist_norm, 'b' ) ; 
        xlabel( {'Segment'; ''} ) ; 
        ylabel( 'Acc mag (LU/TU^2)' ) ; 
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


%% subfunctions 
        






