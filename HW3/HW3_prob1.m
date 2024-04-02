
epsilon = 1/1000 ; 
A       = -3 ; 
B       = 2 ;  
omega   = 5 ; 
K       = 2 ; 
x0      = 1 ; 
xdot0   = 2 ;  
tspan   = [ 0, 10 ] ; 


%% symbolic stuff 

syms c1 c2 t 

g1 = sin(K*t) * ( -1/K * cos(omega*t) - A * ( -c1 * sin(K*t) + c2 * cos(K*t) ) - B/K * ( c1*cos(K*t) + c2*sin(K*t) )^3 ) ; 
g2 = cos(K*t) * ( 1/K * cos(omega*t) - A * ( -c1 * sin(K*t) + c2 * cos(K*t) ) - B/K * ( c1*cos(K*t) + c2*sin(K*t) )^3 ) ; 

% solve for 0'th order 
mat   = [ cos(K*t) sin(K*t) ; -K*sin(K*t) K*cos(K*t) ] ; 
mat   = subs(mat, t, 0) ; 
c_0th = inv(mat) * [ x0 ; xdot0 ] ; 

% solve for 1st order 
g1_eval = subs( g1, c1, c_0th(1) ) ; 
g1_eval = subs( g1_eval, c2, c_0th(2) ) ; 
g2_eval = subs( g2, c1, c_0th(1) ) ; 
g2_eval = subs( g2_eval, c2, c_0th(2) ) ; 

c_1st = [ int(g1_eval,t) ; int(g2_eval,t) ] ; 

% solve for 2nd order 
dg = [ diff(g1, c1), diff(g1, c2) ; diff(g2, c1), diff(g2, c2) ] ; 

dg_eval = subs(dg, c1, c_0th(1)) ; 
dg_eval = subs(dg_eval, c2, c_0th(2)) ; 

c_2nd = 2 * int( dg_eval * c_1st ) ; 

% integrate numerically 

% ode options 
toler   = 1e-8;         % 1e-14 accurate; 1e-6 coarse 
options = odeset('reltol', toler, 'abstol', toler ); 

% integrate  
[t_hist, x_hist] = ode45(@(t, x) duffing(t, x, K, epsilon, omega, A, B), tspan, [x0 xdot0], options); 
plot( t_hist, x_hist(:,1) ) 
title( 'x position' ) ; xlabel('time (s)')  

%% jesus christ ... true c 

c1_0 = xdot0 / K ; 
c2_0 = x0 ; 

[t_hist, c_true] = ode45(@(t, c) cdot_g_fn( t, c, K, omega, A, B, epsilon ), t_hist, [c1_0 c2_0], options); 



%% final solution of different orders 

% 1st order full solution 

c_0th_soln = c_0th ; 
c_1st_soln = c_0th + epsilon * c_1st ; 
c_2nd_soln = c_0th + epsilon * c_1st + epsilon^2 / factorial(2) * c_2nd ; 

c_0th_soln_fn = matlabFunction( c_0th ) ; 
c_1st_soln_fn = matlabFunction( c_1st_soln ) ; 
c_2nd_soln_fn = matlabFunction( c_2nd_soln ) ; 

c_0th_fn = matlabFunction(c_0th) ; 
c1_1st_fn = matlabFunction(c_1st(1)) ; 
c2_1st_fn = matlabFunction(c_1st(2)) ; 
c1_2nd_fn = matlabFunction(c_2nd(1)) ; 
c2_2nd_fn = matlabFunction(c_2nd(2)) ; 

x1_0th_approx_hist = [] ; 
x1_1st_approx_hist = [] ; 
x1_2nd_approx_hist = [] ; 
x1_true_hist       = [] ; 
    
for i = 1 : length(t_hist) 
    
    t = t_hist(i) ; 
    
    % APPROX 
    x1_0th_approx = c_0th_soln_fn()' * [ sin(K*t) ; cos(K*t) ] ; 
    x1_1st_approx = c_1st_soln_fn(t)' * [ sin(K*t) ; cos(K*t) ] ; 
    x1_2nd_approx = c_2nd_soln_fn(t)' * [ sin(K*t) ; cos(K*t) ] ;  
    
    % TRUE 
    x1_true = c_true(i,:) * [ sin(K*t) ; cos(K*t) ] ; 
    
    x1_0th_approx_hist = [ x1_0th_approx_hist ; x1_0th_approx ] ; 
    x1_1st_approx_hist = [ x1_1st_approx_hist ; x1_1st_approx ] ; 
    x1_2nd_approx_hist = [ x1_2nd_approx_hist ; x1_2nd_approx ] ; 
    x1_true_hist = [ x1_true_hist ; x1_true ] ; 
    
end 




%% 12 plots 

figure() 

% c1 
subplot(3,4,1) 
    plot( t_hist, c_0th(1) * ones(size(t_hist)) ) ; 
    title( { '( \epsilon / n! ) c_1^n ' ; 'n = 0 '} ) 
subplot(3,4,5) 
    plot( t_hist, epsilon / 1 * c1_1st_fn( t_hist ) ) ; 
    title( 'n = 1' ) 
subplot(3,4,9) 
    plot( t_hist, epsilon^2 / 2 * c1_2nd_fn( t_hist ) ) ; 
    title( 'n = 2' ) 
    
% c2 
subplot(3,4,2) 
    plot( t_hist, c_0th(2) * ones(size(t_hist)) ) 
    title( { '( \epsilon / n! ) c_2^n ' ; 'n = 0 '} ) 
subplot(3,4,6) 
    plot( t_hist, epsilon / 1 * c2_1st_fn( t_hist ) ) ; 
    title( 'n = 1' ) 
subplot(3,4,10) 
    plot( t_hist, epsilon^2 / 2 * c2_2nd_fn( t_hist ) ) ; 
    title( 'n = 2' ) 

% x1 
subplot(3,4,3) 
    plot( t_hist, x1_0th_approx_hist ) 
    title( { 'x_1^n(t) ' ; 'n = 0 '} ) 
subplot(3,4,7) 
    plot( t_hist, x1_1st_approx_hist ) 
    title( 'n = 1' ) 
subplot(3,4,11) 
    plot( t_hist, x1_2nd_approx_hist ) 
    title( 'n = 2' ) 

% error 
subplot(3,4,4) 
    plot( t_hist, x1_true_hist - x1_0th_approx_hist ) 
    title( { 'error' ; 'n = 0' } ) 
subplot(3,4,8) 
    plot( t_hist, x1_true_hist - x1_1st_approx_hist ) 
    title( 'n = 1' ) 
subplot(3,4,12) 
    plot( t_hist, x1_true_hist - x1_2nd_approx_hist ) ; 
    title( 'n = 2' ) 


%% subfunctions 

function cdot = cdot_g_fn( t, c, K, omega, A, B, epsilon ) 

    c1 = c(1) ; 
    c2 = c(2) ; 

%     g1 = -1/K * sin(K*t) * cos(omega*t) + A * ( -c1 * sin(K*t) + c2 * cos(K*t) ) ; 
%     g2 = 1/K * cos(K*t) * cos(omega*t) - A * ( -c1 * sin(K*t) + c2 * cos(K*t) ) ; 
    
    g1 = sin(K*t) * ( -1/K * cos(omega*t) - A * ( -c1 * sin(K*t) + c2 * cos(K*t) ) - B/K * ( c1*cos(K*t) + c2*sin(K*t) )^3 ) ; 
    g2 = cos(K*t) * ( 1/K * cos(omega*t) - A * ( -c1 * sin(K*t) + c2 * cos(K*t) ) - B/K * ( c1*cos(K*t) + c2*sin(K*t) )^3 ) ; 

    
    cdot = epsilon * [ g1 ; g2 ] ; 

end 


function dxdt = duffing( t, x, K, epsilon, omega, A, B ) 

    dxdt = zeros(2,1) ; 
    dxdt(1) = x(2) ; 
    dxdt(2) = -K^2 * x(1) + epsilon * ( cos(omega*t) - A*x(2) - B*x(1)^3 ) ; 

end 




