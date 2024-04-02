
% get roots 

lambda = @(r) -2 * r^2 * a - R^2 * J2 * a + r * a^2 - r * a^2 * e^2 + r^3 ; 
Ra = fzero( lambda, ra ) ; 
Rp = fzero( lambda, rp ) ; 
R0 = fzero( lambda, 0 ) ; 


%% Problem 2 

% parameters 
R  = 1.1 ; 
J2 = 0.3 ; 
a  = 6 ; 
e  = 0.7 ; 

% Keplerian periapsis and apoapsis 
rp = a * ( 1 - e ) ; 
ra = a * ( 1 + e ) ; 

% initial conditions 
theta0 = 0 ; 
u0     = 1 / Rp  ; 
udot0  = 0 ; 

% ode options 
toler   = 1e-10 ;   % 1e-14 accurate; 1e-6 coarse 
options = odeset( 'reltol', toler, 'abstol', toler ) ; 

% ode loop for integration stopping 
u0_tmp    = u0 ; 
udot0_tmp = udot0 ; 
tspan     = [ 0 theta_full ] ; 
u_hist     = [] ; 
theta_hist = [] ; 
for i = 1 : 7 
    
    [theta, u] = ode45(@(theta, u) planar_J2(theta, u, R, J2, a, e), tspan, [u0_tmp udot0_tmp], options); 
    
    u0_tmp    = u(end,1) ; 
    udot0_tmp = u(end,2) ; 
    
    tspan = [ theta(end) theta(end) + theta_full ] ; 
    
    u_hist = [ u_hist ; u ] ; 
    theta_hist = [ theta_hist ; theta ] ; 
    
end 

% stopping integration 
options = odeset( 'reltol', toler, 'abstol', toler, 'Events', @stop_event ) ; 
% tspan 
[theta, u] = ode45(@(theta, u) planar_J2(theta, u, R, J2, a, e), tspan, [u0_tmp udot0_tmp], options); 

u_hist = [ u_hist ; u ] ; 
theta_hist = [ theta_hist ; theta ] ; 

r_hist_num = 1 ./ u_hist ; 

% plot 
polarplot( theta_hist, r_hist_num(:,1) ) ; 


%% analytic solution 

num = Rp * ( Ra - R0 ) ; 
den = Ra * Rp + R0 * Ra + R0 * Rp ; 
kappa = 1/2 * sqrt( num / den ) ; 

k2 = ( Ra - Rp ) / ( Ra - R0 ) * R0 / Rp ; 
k  = sqrt(k2) ; 

[sn, cn, dn] = ellipj( kappa * theta_hist, k^2 ) ; 

num = R0 * Rp * dn.^2 ; 
den = R0 - Rp * k^2 * sn.^2 ; 
r_hist_ana = num ./ den ;  

% half period 
Fk_tilde    = ellipke( k^2 ) ; 
theta_a     = Fk_tilde / kappa ; 
theta_full  = 2 * theta_a ; 

% perihelion advance 
domega = 2 * Fk_tilde / kappa - 2 * pi ; 


%% plot 

close all ; 

figure() 
polarplot( theta_hist, r_hist_num(:,1) ) ; hold on ; grid on ; 
polarplot( theta_hist, r_hist_ana, '--' ) 
legend('numerical', 'analytical') 

figure() 
error = r_hist_ana - r_hist_num(:,1) ; 
plot( theta_hist, error )
title( 'Analytical - Numerical Error' )
xlabel('Theta') 


%% subfunctions 

function [ value, isterminal, direction ] = stop_event( theta, u ) 

    r  = 1 / u(1) ; 
    Rp = 1.64525462898384 ; 

%     value = [ ( theta - 55 ) ; abs(r - Rp) - 1e-5 ] ; 
%     isterminal = [ 1 ; 1 ] ; 
%     direction  = [ 0 ; 0 ] ; 

    value      = ( r - Rp - 1e-5 ) ;
    isterminal = 1;   % Stop the integration
    direction  = 0;

end 

function output = state_dot_J2( t, state, R, J2, a, e )  

u = state(1) ; 
udot = state(2) ; 

u_der = udot ; 
udot_der = 1/(a*(1-e^2)) + (3*(R^2)*J2*(u^2))/(2*a*(1-e^2)) - u ; 

output = [ u_der; udot_der] ; 

end 

function r = analytic_J2( theta, kappa, k, R0, Rp )  

end 

function [ R0, Rp, Ra ] = lambda_fn( rp, ra, a, R, J2 )
    
    a = 1 ; 
    b = - ( ra + rp ) ; 
    c = ra * rp ; 
    d = a * R^2 * J2 ; 
    p = [ a b c d ] ; 

    R0_Rp_Ra = roots( p ) ; 
    
    Ra = max( R0_Rp_Ra ) ; 

    R0_Rp = R0_Rp_Ra( R0_Rp_Ra ~= Ra ) ; 
    Rp = max( R0_Rp ) ; 
    R0 = min( R0_Rp ) ; 
    
end 

function dudtheta = planar_J2( theta, u, R, J2, a, e ) 
    
    % central force 
    num = 1 + 3/2 * R^2 * J2 * u(1)^2 ; 
    den = a * (1-e^2) ; 
    
    dudtheta = zeros(2,1) ; 
    dudtheta(1) = u(2) ; 
    dudtheta(2) = num / den - u(1) ; 

end 

function drv = central_force( t, rv, k, c ) 
% Integrate equations of motion for Cartesian state and time velocity 

    r = rv(1:3) ; r_hat = r / norm(r) ; 
    v = rv(4:6) ; 
    
    % central force 
    P = k / norm(r)^2 - c / norm(r)^3 ; 

    % r, v, derivatives wrt time 
    dr = v ; 
    dv = -P * r_hat ; 
    
    % output  
    drv = [ dr ; dv ] ; 
    
end 

