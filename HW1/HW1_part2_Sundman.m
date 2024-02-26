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

%% loop for rv_t_EOM 

for alpha = [ 0, 1, 2, 3/2 ]

    % ----------------------- % 
    % part 2.a.1: Convert the equations of motion to a new independent 
    % variable τ using the general Sundman transformation 

    % ----------------------- % 
    % part 2.a.2: numeric value of τ for 1 period of the ellipse (τP)?

    [ T_time, T_tau ] = T_tau_fn( oe0, k, mu, alpha ) ; 
    sprintf( 'n = %g: T_tau = %.6g', alpha, T_tau ) 

    % ----------------------- % 
    % part 2.a.3: Using the general ODE solver provided on canvas to solve the 
    % equations of motion for τ=0.. τP . Use a fixed step integration in tau 
    % Δτ = τP/q for q=50. If you don’t use Matlab, you can use any ODE solver 
    % as long as it is very clearly a fixed step integrator.

    q    = 50 ; 
    dtau = T_tau / q ; 

    % initial state = r, v, t 
    state0 = [ rv0 ; 0 ] ; 

    % zero perturbation 
    a_pert = zeros(3,1) ; 

    [tau_hist, state] = ode78rpr( @(tau, rvt) rv_t_EOM( tau, rvt, a_pert, mu, k, alpha ), ... 
        0, T_tau, dtau, state0, tol) ; 
    t_hist = state(:,7) ; 

    % ----------------------- % 
    % part 2.a.4: Plot the XY view of the trajectory including nodes and 
    % connecting lines. 

    figure() 
        plot_rv( state, sprintf( 'Cartesian trajectory with time velocity for n = %g', alpha ) ) ; 
        view(0,90) 

    % ----------------------- % 
    % part 2.a.5: Plot the first six states vs. time and tau using subplots (2 
    % columns, 6 rows) in a single figure.

    titles = { 'r_x', 'r_y', 'r_z', 'v_x', 'v_y', 'v_z' } ; 
    ylabels = { 'LU', 'LU', 'LU', 'LU/TU', 'LU/TU', 'LU/TU' } ; 
    
    figure() 
        for i = 1 : 6 
            subplot(2,3,i) 
                plot( t_hist, state(:,i) ) ; 
                title( titles{i} ) ; 
                if i > 3 ; xlabel('t (s)') ; end 
                if i == 1 ; ylabel('LU') ; elseif i == 4 ; ylabel('LU/TU') ; end 
        end 
        sgtitle( sprintf( 'time vs rv with time velocity for n = %g', alpha ) ) ; 
    
    figure() 
        for i = 1 : 6 
            subplot(2,3,i) 
                plot( tau_hist, state(:,i) ) ; 
                title( titles{i} ) ; 
                if i > 3 ; xlabel('tau') ; end 
                if i == 1 ; ylabel('LU') ; elseif i == 4 ; ylabel('LU/TU') ; end 
        end 
        sgtitle( sprintf( 'tau vs rv with time velocity for n = %g', alpha ) ) ; 

    % ----------------------- % 
    % part 2.a.6: Report |r0-rf | |v0-vf | and |tf-2πsqrt(a^3/μ)|
    r0 = state(1, 1:3)' ; rf = state(end, 1:3)' ; 
    v0 = state(1, 4:6)' ; vf = state(end, 4:6)' ; 
    tf = state(end, end) ; tf_ana = 2*pi*sqrt( a^3 / mu ) ;  
    sprintf( 'n = %g: | r0 - rf | = %.6g ', alpha, norm( r0 - rf ) ) 
    sprintf( 'n = %g: | v0 - vf | = %.6g ', alpha, norm( v0 - vf ) ) 
    sprintf( 'n = %g: | tf-2*pi*sqrt(a^3/mu) | = %.6g ', alpha, norm( tf - tf_ana ) ) 

end 


%% loop for rv_tau_EOM 

% loop 
for alpha = [ 0, 1, 2, 3/2 ]
% for alpha = 3/2 

    % oh god at least I have the ICs saved here 
    r0 = rv0(1:3) ; 
    v0 = rv0(4:6) ; 
    
    % initial state = r, v_tau, t 
    v_tau0 = v0 * k * norm(r0)^(alpha) ; 
    state0 = [ r0; v_tau0 ; 0 ] ; 

    % ----------------------- % 
    % part 2.a.1: Convert the equations of motion to a new independent 
    % variable τ using the general Sundman transformation 

    % ----------------------- % 
    % part 2.a.2: numeric value of τ for 1 period of the ellipse (τP)?

    [ T_time, T_tau ] = T_tau_fn( oe0, k, mu, alpha ) ; 
    sprintf( 'n = %.2g: T_tau = %.6g', alpha, T_tau ) 

    % ----------------------- % 
    % part 2.a.3: Using the general ODE solver provided on canvas to solve the 
    % equations of motion for τ=0.. τP . Use a fixed step integration in tau 
    % Δτ = τP/q for q=50. If you don’t use Matlab, you can use any ODE solver 
    % as long as it is very clearly a fixed step integrator.

    q    = 50 ; 
    dtau = T_tau / q ; 

    % zero perturbation 
    a_pert = zeros(3,1) ; 

    [tau_hist, state] = ode78rpr( @(tau, rvt) rv_tau_EOM( tau, rvt, a_pert, mu, k, alpha ), ... 
        0, T_tau, dtau, state0, tol) ; 
    t_hist = state(:,7) ; 

    % ----------------------- % 
    % part 2.a.4: Plot the XY view of the trajectory including nodes and 
    % connecting lines. 

    figure() 
        plot_rv( state, sprintf( 'Cartesian trajectory with tau velocity for n = %g', alpha ) ) ; 
        view(0,90) 

    % ----------------------- % 
    % part 2.a.5: Plot the first six states vs. time and tau using subplots (2 
    % columns, 6 rows) in a single figure.

    titles = { 'r_x', 'r_y', 'r_z', 'v_x', 'v_y', 'v_z' } ; 
    ylabels = { 'LU', 'LU', 'LU', 'LU/TU', 'LU/TU', 'LU/TU' } ; 
    
    figure() 
        for i = 1 : 6 
            subplot(2,3,i) 
                plot( t_hist, state(:,i) ) ; 
                title( titles{i} ) ; 
                if i > 3 ; xlabel('t (s)') ; end 
                if i == 1 ; ylabel('LU') ; elseif i == 4 ; ylabel('LU/TU') ; end 
        end 
        sgtitle( sprintf( 'time vs rv with tau velocity for n = %g', alpha ) ) ; 
        
    figure() 
        for i = 1 : 6 
            subplot(2,3,i) 
                plot( tau_hist, state(:,i) ) ; 
                title( titles{i} ) ; 
                if i > 3 ; xlabel('tau') ; end 
                if i == 1 ; ylabel('LU') ; elseif i == 4 ; ylabel('LU/TU') ; end 
        end 
        sgtitle( sprintf( 'tau vs rv with tau velocity for n = %g', alpha ) ) ; 

    % ----------------------- % 
    % part 2.a.6: Report |r0-rf | |v0-vf | and |tf-2πsqrt(a^3/μ)|
    r0 = state(1, 1:3)' ; rf = state(end, 1:3)' ; 
    v0 = state(1, 4:6)' ; vf = state(end, 4:6)' ; 
    tf = state(end, end) ; tf_ana = 2*pi*sqrt( a^3 / mu ) ;  
    sprintf( 'n = %g: | r0 - rf | = %.6g ', alpha, norm( r0 - rf ) ) 
    sprintf( 'n = %g: | v0 - vf | = %.6g ', alpha, norm( v0 - vf ) ) 
    sprintf( 'n = %g: | tf-2*pi*sqrt(a^3/mu) | = %.6g ', alpha, norm( tf - tf_ana ) ) 

end 

    
%% subpart 2C: Using the simulation from subpart 2A, perform a trade study to examine the effect of accuracy due
% to changing q and eccentricity from 0 to 0.99. Generate 7 plots, each one dealing with a separate
% value of: q=20, 40, 80, 160, 320, 640, 1280.
% In each fixed-q plot, include 5 curves, one for each of the following {n=0, n=1, n=3/2, n=2, n=3}.
% The plot should have as its y axis: miss distance after one period = |r0-rf |, and its x axis: initial
% eccentricity. If a particular data point result is NaN, simply skip that point. Use semiology and
% make the ylim[1e-16 1]). Comment on the results.

% repeating stuff from subpart 2A 

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

% left-hand side terms 
n = sqrt( mu / a^3 ) ; 
p_pb = a * (1 - e^2) ; 
tau_terms = n * k * a^2 * p_pb * ( 1-e^2 )^(1/2) ; 

% right-hand side term (integrated) 
nu = 2*pi ; 
nu_terms = nu + e * sin(nu) ; 

alpha = 3 ; 

T_tau = nu_terms / tau_terms ; 
% [ T_time, T_tau ] = T_tau_fn( oe0, k, mu, alpha ) ; 

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


%% subpart 2C: Using the simulation from subpart 2A, perform a trade study 
% to examine the effect of accuracy due to changing q and eccentricity 
% from 0 to 0.99. Generate 7 plots, each one dealing with a separate
% value of: q=20, 40, 80, 160, 320, 640, 1280.
% In each fixed-q plot, include 5 curves, one for each of the following 
% {n=0, n=1, n=3/2, n=2, n=3}.
% The plot should have as its y axis: miss distance after one period = 
% |r0-rf |, and its x axis: initial eccentricity. If a particular data 
% point result is NaN, simply skip that point. Use semiology and
% make the ylim[1e-16 1]). Comment on the results.


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

% initial state = r, v, t 
state0 = [ rv0 ; 0 ] ; 

% value of k 
k = 1 ; 

% loop for rv_t_EOM 

for q = [ 20, 40, 80, 160, 320, 640, 1280 ] 
% for q = [ 1280 ] 
    
    figure() ; 
    for alpha = [ 0, 1, 3/2, 2, 3 ]
%     for alpha = 3 

        e_hist = [] ; pos_err_hist = [] ; 
        for e = [ 0 : 0.01 : 0.99 ] 

            % orbital elements and cartesian 
            oe0     = [ a; e; i; w; Omega; nu ] ; 
            rv0     = oe2rv( oe0, mu ) ; 

            % initial state = r, v, t 
            state0 = [ rv0 ; 0 ] ; 
            
            % numeric value of τ for 1 period of the ellipse (τP)?
            [ T_time, T_tau ] = T_tau_fn( oe0, k, mu, alpha ) ; 
%             sprintf( 'T_tau = %.6g for n = %.2d', T_tau, alpha ) 

            % ----------------------- % 

            dtau = T_tau / q ; 

            % zero perturbation 
            a_pert = zeros(3,1) ; 

            % propagate 
            [tau_hist, state] = ode78rpr( @(tau, rvt) rv_t_EOM( tau, rvt, a_pert, mu, k, alpha ), ... 
                0, T_tau, dtau, state0, tol) ; 

            % pos err 
            r0 = state(1, 1:3) ; rf = state(end, 1:3) ; 
            pos_err = norm( r0 - rf ) ; 
            
            e_hist = [ e_hist ; e ] ; 
            pos_err_hist = [ pos_err_hist ; pos_err ] ; 
            
        end 
        
        semilogy( e_hist, pos_err_hist ) ; hold on ; 
        ylim([1e-16 1]) ; 
        
    end 
    
    % finish plot 
    legend( 'n = 0', 'n = 1', 'n = 3/2', 'n = 2', 'n = 3', 'location', 'northwest' ) ; 
    plot_title = sprintf( 'Eccentricity vs miss distance q = %d ', ... 
        q ) ; 
    title(plot_title) 
    ylabel('Miss distance after one period | r0 - rf |')
    xlabel('e') 
        
end 


%% subpart 2D: Redo 2A-3 and 2A-4 for initial conditions r0=[1,0,0] LU, v0=[0, p,0] LU/TU, GM=1, and flight
% time of 10 TU in physical time for a) parabola (p=sqrt(2) ) and b) hyperbola (p=2). Below is a table
% for the tau of flight for each case that should correspond approximately 10 TU in time. You can
% hard code the stop on tau using this table, or use an approximate event finder. Comment on the
% results.

% parabola 
p_pb = sqrt(2) ; 
p_hb = 2 ; 

q     = 50 ; 

% zero perturbation 
a_pert = zeros(3,1) ; 

for alpha = [ 0, 1, 3/2, 2 ] 
    
    if alpha == 0 
        T_tau_pb = 10 ; 
        T_tau_hb = 10 ; 
    elseif alpha == 1 
        T_tau_pb = 3.407263066 ; 
        T_tau_hb = 2.150490755 ; 
    elseif alpha == 3/2 
        T_tau_pb = 2.281138482 ; 
        T_tau_hb = 1.297868738 ; 
    else 
        T_tau_pb = 1.665061454 ; 
        T_tau_hb = 0.909576946 ; 
    end 
    dtau_pb = T_tau_pb / q ; 
    dtau_hb = T_tau_hb / q ; 
    
    % PARABOLA 

    r0 = [1; 0; 0] ; v0 = [0; p_pb; 0] ; rv0 = [ r0 ; v0 ] ; 
    state0 = [ rv0 ; 0 ] ; 

    % part 2a-3: propagate 
    [tau_hist, state] = ode78rpr( @(tau, rvt) rv_t_EOM( tau, rvt, a_pert, mu, k, alpha ), ... 
        0, T_tau_pb, dtau_pb, state0, tol) ; 

    figure() 
        plot_rv( state, sprintf( 'Parabolic n = %d', alpha ) ) 
        view(0,90) 
    
    % HYPERBOLA

    r0 = [1; 0; 0] ; v0 = [0; p_hb; 0] ; rv0 = [ r0 ; v0 ] ; 
    state0 = [ rv0 ; 0 ] ; 

    % part 2a-3: propagate 
    [tau_hist, state] = ode78rpr( @(tau, rvt) rv_t_EOM( tau, rvt, a_pert, mu, k, alpha ), ... 
        0, T_tau_hb, dtau_hb, state0, tol) ; 

    figure() 
        plot_rv( state, sprintf( 'Hyperbolic n = %d', alpha ) ) 
        view(0,90) 
    
end 

















