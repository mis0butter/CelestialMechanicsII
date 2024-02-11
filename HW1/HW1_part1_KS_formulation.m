%% PART 1: KS FORMULATION 

close all 
init 

%% PART 1.1: UNPERTURBED 


%% part a: Code the transformation subroutines for transformations in both 
% directions from { u, u' } <--> { r, v }. Give a screenshot of the 
% subroutines. (separate from providing src in the appendix)

a       = 1.1 ;             % LU 
e       = 0.9 ; 
i       = 40 * pi/180;      % deg
w       = 10e-15 * pi/180;  % deg 
Omega   = 10e-15 * pi/180;  % deg 
nu      = 10e-15 * pi/180;  % deg 
mu      = 1.25 ;            % LU^3 / TU^2 

% orbital elements and cartesian 
oe0     = [ a; e; i; w; Omega; nu ] ; 
rv0     = oe2rv( oe0, mu ) ; 

r0 = rv0(1:3) ; 
v0 = rv0(4:6) ; 

% { u, u' } <--> { r, v } 
KS0 = rv2KS( rv0 ) ; 
rv  = KS2rv( KS0 ) ;    % check 


%% part b: Code the analytic solution that takes the KS variables, energy 
% and tau as inputs and converts them to the KS variables and time at the 
% specified tau.

u0      = KS0(1:4) ; 
uprime0 = KS0(5:8) ; 
t0      = 0 ; 
eps0    = 2 * dot( uprime0, uprime0 ) / ( dot( u0, u0 ) ) - mu / ( dot( u0, u0 ) ) ; 

% dummy test 
tau     = 2 ; 
KS_time = KS_time_fn( u0, uprime0, tau, eps0 ) ;  


%% part c: What is the theoretical tau period? What is the initial 
% cartesian state? What is the initial KS state?

% period in time (s) 
T = 2*pi*sqrt( a^3 / mu ) ; 
    
% period in tau (tau units) 
k = 1 ; 
T_tau = T / ( k * a ) ; 

% check your answer 
t_hist   = [] ; 
tau_hist = [] ; 
for tau = 1 : 0.1 : 2*T 
    
    KS_time = KS_time_fn( u0, uprime0, tau, eps0 ) ; 
    t = KS_time(end) ; 
    
    t_hist   = [ t_hist ; t ] ; 
    tau_hist = [ tau_hist ; tau ] ; 
    
end 

clf
figure() ; hold on ; grid on ; 
    plot( t_hist, t_hist ) ; 
    plot( t_hist, T * ones(size(t_hist)), 'red' ) ; 
    plot( tau_hist, t_hist, 'green' ) ;
    legend( 't', 'T', 'tau' ) ; 
    ylabel('t (s)') 
    
    
%% part d: Analytically propagate the orbit for 60 equal steps in tau for 1 
% period. For each step, traverse from the initial conditions all the way 
% to the intermediate tau.

% tau steps 
N = 60 ; 

% analytically propagate unperturbed KS 
tau_hist   = [ 0 : T_tau / N : T_tau ] ; 
state_hist = [] ; 
for tau = tau_hist  
    u_uprime_t = KS_time_fn( u0, uprime0, tau, eps0 ) ; 
    state_hist = [ state_hist ; u_uprime_t' ] ; 
end 

% extract 
u_hist      = state_hist(:, 1:4) ; 
uprime_hist = state_hist(:, 5:8) ; 
t_hist      = state_hist(:, 9) ; 
% eps_hist    = state(:, 10) ; 

% get Cartesian hist 
rv_hist = [] ; 
for i = 1 : length(u_hist) 
    KS      = [ u_hist(i,:)' ; uprime_hist(i,:)' ] ; 
    rv      = KS2rv( KS )' ;    
    rv_hist = [ rv_hist ; rv ] ; 
end 

% ----------------------- % 
% i: Confirm the time period matches the expected value from the two-body
% formula.

t_period = t_hist(end) ; 
sprintf( "tau propagated t_period = %.10g s", t_period ) 
sprintf( "T = %.10g s", T ) 


% ----------------------- % 
% ii: Plot the Cartesian trajectory, with clear markers at the nodes and 
% lines connecting them. Use >>view(30,30) to see the 3D component. 

plot_rv(rv_hist, 'Analytical Cartesian Trajectory') ; 
    
    
% ----------------------- % 
% iii: Plot each of the output {u, u'} states as a function of tau on a 
% single subplot. 

% tau_hist = [ 0 : T_tau / 60 : T_tau ]' ; 
% tau_hist = [ tau_hist ; T_tau ] ; 

% okay, t_hist, u_hist and uprime_hist last 2 rows are almost same ...
% whatever 

figure(2) ; 
    plot( tau_hist, u_hist ) ; 
    hold on ; grid on ; 
    plot( tau_hist, uprime_hist, '--' ) ; 
    legend( 'u1', 'u2', 'u3', 'u4', 'uprime1', 'uprime2', 'uprime3', 'uprime4', 'location', 'eastoutside' ) ; 
    xlabel('\tau') ; ylabel('u, uprime') ; 
    title( 'Analytical \tau vs \{ u, uprime \}' ) ; 
    
    
% ----------------------- %  
% iv: Plot the output time as a function of tau on a single plot.

figure(3) ; hold on ; grid on ; 
    plot( tau_hist, t_hist ) ; 
    xlabel('\tau') ; ylabel('t (s)') 
    title( 'Analytical \tau vs t' ) ; 

  
% ----------------------- %      
% v: Compute the K numerical check at each step, and plot the results on a
% single plot 

[ K_A_hist, K_B_hist ] = K_check( state_hist, mu ) ; 

figure() ; hold on ; grid on ; 
    scatter( tau_hist, K_A_hist ) ; 
    scatter( tau_hist, K_B_hist, 'x' ) ; 
    legend( 'method 1', 'method 2', 'location', 'northwest' ) ; 
    xlabel('\tau') ; ylabel('K') ; 
    title( 'Analytical \tau vs K' ) ; 
    
    
% ----------------------- %      
% vi. Report | u0 - uf | and | u'0 - u'f | 

disp('[ u_hist(1,:) - u_hist(end,:) ] = ') 
disp([ u_hist(1,:) - u_hist(end,:) ]')

disp('[ uprime_hist(1,:) - uprime_hist(end,:) ] = ')
disp([ uprime_hist(1,:) - uprime_hist(end,:) ]')

disp('[ rv_hist(1,:) - rv_hist(end,:) ] = ') 
disp([ rv_hist(1,:) - rv_hist(end,:) ]') 
    
    
%% part e: Code a subroutine to numerically propagate, using the KS 
% formulation, 1 period of the same orbit for N equal steps in tau, using 
% the fixed step integrator. Repeat the following for N={3,5,10,20,40,80}

%tol should very from 1e-17 fine to 1e-4 course (or 0.0 for fixed step
tol = 0.0;         

% 0 perturbation 
a = zeros(4,1) ; 

% Solve ODE 
state0   = [ KS0 ; t0 ; eps0 ] ; 

% loop for tau steps 
for N = [ 3, 5, 10, 20, 40, 80 ] 
    
    % propagate 
    [ state_num, rv_num ] = prop_KS_num( state0, a, mu, T_tau, N, tol ) ; 
    [ state_ana, rv_ana ] = prop_KS_ana( state0, T_tau, N ) ; 

    % i: Plot the K numerical check as a function of tau 
    
    % ii: On the same plot as i, give the norms of the position, velocity and 
    % time differences (total of 4 curves on the one plot) between the 
    % analytical solution using the subroutine from 1b) and the numerically 
    % integrated solution with N equal steps. (don't use log scale since K can 
    % be negative, and put N in the title.) 
    [ K_A_num, K_B_num ]  = K_check( state_num, mu, 1, N, 'Numerical' ) ;   
%     [ K_A_ana, K_B_ana ]  = K_check( state_ana, mu, 1, N, 'Analytical' ) ;   

    % extract 
    r_num = rv_num(:, 1:3) ; v_num = rv_num(:, 4:6) ; 
    r_ana = rv_ana(:, 1:3) ; v_ana = rv_ana(:, 4:6) ; 
    t_num = state_num(:, 9) ; 
    t_ana = state_ana(:, 9) ; 

    % position, velocity, and time error norm 
    pos_err = [] ; vel_err = [] ; t_err = [] ; 
    for i = 1 : N + 1 
        pos_err(i,:) = norm( r_num(i,:) - r_ana(i,:) ) ; 
        vel_err(i,:) = norm( v_num(i,:) - v_ana(i,:) ) ; 
        t_err(i,:)   = norm( t_num(i,:) - t_ana(i,:) ) ; 
    end 
    
    tau_hist = state_num(:,11) ; 
    plot( tau_hist, pos_err ) ;
    plot( tau_hist, vel_err ) ; 
    plot( tau_hist, t_err ) ; 
    
    legend( 'K method 1', 'K method 2', '|pos err|', '|vel err|', '|t err|', ... 
        'location', 'northwest' ) ; 
    
end 


% save all figs 

% folder_name = 'outputs';   % Your destination folder
% save_all_figs(folder_name) 


%% PART 1.2: PERTURBED
% Add the perturbation acceleration a to your numerical simulation of KS.
% Include energy as the 10th state. You can use the examples in the notes 
% to validate your code. 


%% part a: Analytically, prove that that epsilon_prime = 0 when dot(a, 
% rdot) = 0. Validate with your numerical code. 

% set velocity to 0 
rv_v0       = rv0 ; 
rv_v0(4:6)  = zeros(3,1) ; 
KS_v0       = rv2KS( rv_v0 ) ; 

% get tau derivatives with nonzero perturbation 
a           = rand(4,1) ; 
state_prime = KS_EOM( tau, KS_v0, a, mu ) ; 
eps_prime   = state_prime(end)  


%% part b: Apply a thrust in the VTN frame, using a thrust magnitude of 
% 0.03 LU/TU2, longitude=60 deg, latitude=30 deg. Propagate for 4 periods 
% (using the IC’s assuming ballistic motion) and 120 total equal steps in 
% tau. Use the same fixed step integrator.

% a in VTN frame 
Tmag  = 0.005 ; 
lon   = 0 * pi/180 ; 
lat   = 90 * pi/180 ; 
a_VTN = latlon2VTN( Tmag, lat, lon ) ; 

a       = 1 ;             % LU 
e       = 0.7 ; 
i       = 40 * pi/180;      % deg
w       = 0 ;  % deg 
Omega   = 0 ;  % deg 
nu      = 0 ;  % deg 
mu      = 1 ;            % LU^3 / TU^2 

% orbital elements and cartesian 
oe0     = [ a; e; i; w; Omega; nu ] ; 
rv0     = oe2rv( oe0, mu ) ; 

r0 = rv0(1:3) ; 
v0 = rv0(4:6) ; 

% { u, u' } <--> { r, v } 
KS0 = rv2KS( rv0 ) ; 
rv  = KS2rv( KS0 ) ;    % check 

u0      = KS0(1:4) ; 
uprime0 = KS0(5:8) ; 
t0      = 0 ; 
eps0    = 2 * dot( uprime0, uprime0 ) / ( dot( u0, u0 ) ) - mu / ( dot( u0, u0 ) ) ; 

% period in time (s) 
T = 2*pi*sqrt( a^3 / mu ) ; 
    
% period in tau (tau units) 
k = 1 ; 
T_tau = T / ( k * a ) ; 

state0  = [ u0 ; uprime0 ; t0 ; eps0 ] ; 

[ state_num, rv_num ] = prop_KS_num( state0, a_VTN, mu, T_tau * 10, N * 10, tol ) ; 
plot_rv( rv_num ) ; 














