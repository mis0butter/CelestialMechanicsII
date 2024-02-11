%% PART 2: SUNDMAN 

close all 
init 


%% part a: Code the transformation subroutines for transformations in both 
% directions from { u, u' } <--> { r, v }. Give a screenshot of the 
% subroutines. (separate from providing src in the appendix)

a       = 1.0 ;             % LU 
e       = 0.9 ; 
i       = 10e-15 * pi/180;  % deg 
w       = 10e-15 * pi/180;  % deg 
Omega   = 10e-15 * pi/180;  % deg 
nu      = 10e-15 * pi/180;  % deg 
mu      = 1.0 ;             % LU^3 / TU^2 

% orbital elements and cartesian 
oe0     = [ a; e; i; w; Omega; nu ] ; 
rv0     = oe2rv( oe0, mu ) ; 

r0 = rv0(1:3) ; 
v0 = rv0(4:6) ; 

% { u, u' } <--> { r, v } 
KS0 = rv2KS( rv0 ) ; 
rv  = KS2rv( KS0 ) ;    % check 

% value of k 
k = 1 ; 

% start with alpha 
alpha = 0 ; 

for alpha = [ 0, 1, 2, 3/2 ]

    % ----------------------- % 
    % part 2.a.1: Convert the equations of motion to a new independent 
    % variable τ using the general Sundman transformation 


    % ----------------------- % 
    % part 2.a.2: numeric value of τ for 1 period of the ellipse (τP)?

    [ T_time, T_tau ] = T_tau_fn( oe0, k, mu, alpha ) ; 
    sprintf( 'T_tau = %.6g for n = %.2d', T_tau, alpha ) 


    % ----------------------- % 
    % part 2.a.3: Using the general ODE solver provided on canvas to solve the 
    % equations of motion for τ=0.. τP . Use a fixed step integration in tau 
    % Δτ = τP/q for q=50. If you don’t use Matlab, you can use any ODE solver 
    % as long as it is very clearly a fixed step integrator.

    q    = 40 ; 
    dtau = T_tau / q ; 

    % initial state = r, v, t 
    state0 = [ rv0 ; 0 ] ; 

    % zero perturbation 
    a_pert = zeros(3,1) ; 

    [tau_hist, state] = ode78rpr( @(tau, rvt) rv_t_EOM( tau, rvt, a_pert, mu, k, alpha ), ... 
        0, T_tau, dtau, state0, tol) ; 


    % ----------------------- % 
    % part 2.a.4: Plot the XY view of the trajectory including nodes and 
    % connecting lines. 

    figure() 
        plot_rv( state, sprintf( 'Cartesian trajectory for n = %d', alpha ) ) ; 
%         view(0,90) 


    % ----------------------- % 
    % part 2.a.5: Plot the first six states vs. time and tau using subplots (2 
    % columns, 6 rows) in a single figure.

    titles = { 'r_x', 'r_y', 'r_z', 'v_x', 'v_y', 'v_z' } ; 
    ylabels = { 'LU', 'LU', 'LU', 'LU/TU', 'LU/TU', 'LU/TU' } ; 
    figure() 
        for i = 1 : 6 
            subplot(6,1,i) 
                plot( tau_hist, state(:,i) ) ; 
                title( titles{i} ) ; 
                ylabel( ylabels{i} ) ; 
        end 
        xlabel('tau') 
        sgtitle( sprintf( 'tau vs rv for n = %d', alpha ) ) ; 


    % ----------------------- % 
    % part 2.a.6: Report |r0-rf | |v0-vf | and |tf-2πsqrt(a^3/μ)|
    r0 = state(1, 1:3) ; rf = state(end, 1:3) ; 
    v0 = state(1, 4:6) ; vf = state(end, 4:6) ; 
    tf = state(end, end) ; tf_ana = 2*pi*sqrt( a^3 / mu ) ;  
    sprintf( '| r0 - rf | = %.6g ', norm( r0 - rf ) ) 
    sprintf( '| v0 - vf | = %.6g ', norm( v0 - vf ) ) 
    sprintf( '| tf-2πsqrt(a^3/μ) | = %.6g ', norm( tf - tf_ana ) ) 

end 

    
%% save all figs 












