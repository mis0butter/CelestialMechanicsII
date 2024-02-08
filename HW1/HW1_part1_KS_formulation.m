%% PART 1: KS FORMULATION 

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
eps     = 2 * dot( uprime0, uprime0 ) / ( dot( u0, u0 ) ) - mu / ( dot( u0, u0 ) ) ; 

% dummy test 
tau     = 2 ; 
KS_time = KS_time_fn( u0, uprime0, tau, eps ) ;  


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
    
    KS_time = KS_time_fn( u0, uprime0, tau, eps ) ; 
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

close all 

tau_hist   = [ 0 : T_tau / 60 : T_tau ] ; 
state_hist = [] ; 
for tau = tau_hist  

    u_uprime_t = KS_time_fn( u0, uprime0, tau, eps ) ; 
    state_hist = [ state_hist ; u_uprime_t' ] ; 
    
end 

% extract 
u_hist      = state_hist(:, 1:4) ; 
uprime_hist = state_hist(:, 5:8) ; 
t_hist      = state_hist(:, 9) ; 
% eps_hist    = state(:, 10) ; 

% get Cartesian hist 
rv_hist = [] 
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

figure(1) ; 
    plot3( rv_hist(:,1), rv_hist(:,2), rv_hist(:,3) ) ; 
    hold on ; grid on ; 
    scatter3( rv_hist(:,1), rv_hist(:,2), rv_hist(:,3) ) ; 
    view(30, 30) 
    xlabel('x (LU)') ; ylabel('y (LU)') ; zlabel('z (LU)') ; 
    title('Cartesian trajectory ') ; 

    
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
    title( '\tau vs \{ u, uprime \}' ) ; 
    
    
% ----------------------- %  
% iv: Plot the output time as a function of tau on a single plot.

figure(3) ; hold on ; grid on ; 
    plot( tau_hist, t_hist ) ; 
    xlabel('\tau') ; ylabel('t (s)') 
    title( '\tau vs t' ) ; 

  
% ----------------------- %      
% v: Compute the K numerical check at each step, and plot the results on a
% single plot 

K_A_hist = [] ; 
K_B_hist = [] ; 
for i = 1 : length(tau_hist) 
    
%     eps         = eps_hist(i,:) ; 
    u           = u_hist(i,:) ; 
    uprime      = uprime_hist(i,:) ; 
    
    K_A         = - eps/2 * dot(u,u) + dot(uprime, uprime) - mu/2 ;     
    K_A_hist    = [ K_A_hist ; K_A ] ;     
    
    K2          = eps - 2 * dot( uprime,uprime ) / dot( u,u ) + mu / dot( u,u ) ; 
    K_B         = - dot( u,u ) / 2 * K2 ; 
    K_B_hist    = [ K_B_hist ; K_B ] ; 
    
end 

figure(4) ; hold on ; grid on ; 
    scatter( tau_hist, K_A_hist ) ; 
    scatter( tau_hist, K_B_hist, 'x' ) ; 
    legend( 'method 1', 'method 2', 'location', 'northwest' ) ; 
    xlabel('\tau') ; ylabel('K') ; 
    title( '\tau vs K' ) ; 
    
    
% ----------------------- %      
% vi. Report |u0- uf | and | u�0- u�f | 

disp('[ u_hist(1,:) - u_hist(end,:) ] = ') 
disp([ u_hist(1,:) - u_hist(end,:) ]')

disp('[ uprime_hist(1,:) - uprime_hist(end,:) ] = ')
disp([ uprime_hist(1,:) - uprime_hist(end,:) ]')

disp('[ rv_hist(1,:) - rv_hist(end,:) ] = ') 
disp([ rv_hist(1,:) - rv_hist(end,:) ]') 
    
    
%% part e: Code a subroutine to numerically propagate, using the KS 
% formulation, 1 period of the same orbit for N equal steps in tau, using 
% the fixed step integrator. Repeat the following for N={3,5,10,20,40,80}

clf 

% options 
tol     = 1e-14;         % 1e-14 accurate; 1e-6 coarse 
options = odeset('reltol', tol, 'abstol', tol ); 

% 0 perturbation 
a = zeros(4,1) ; 

% Solve ODE 
state0   = [ KS0 ; t0 ; eps ] ; 
tau_span = [ 0 T_tau ] ; 
% [tau_hist, state] = ode45(@(tau, KS_t_eps) KS_EOM( tau, KS_t_eps, a, mu), tau_span, state0, options) ; 
[tau_hist, state] = ode78rpr(@(tau, KS_t_eps) KS_EOM( tau, KS_t_eps, a, mu), 0, T_tau, T_tau / 60, state0, 0.0) ; 

% extract 
u_hist      = state(:, 1:4) ; 
uprime_hist = state(:, 5:8) ; 
t_hist      = state(:, 9) ; 
eps_hist    = state(:, 10) ; 

% get Cartesian hist 
rv_hist = [] 
for i = 1 : length(u_hist) 
    KS      = [ u_hist(i,:)' ; uprime_hist(i,:)' ] ; 
    rv      = KS2rv( KS )' ;    
    rv_hist = [ rv_hist ; rv ] ; 
end 
    




%% subfunctions 

function KS_zero = KS_zero_fn( u0, uprime0, tau, eps ) 

    KS_time = KS_time_fn( u0, uprime0, tau, eps ) ; 
    t = KS_time(end) ; 


end 





