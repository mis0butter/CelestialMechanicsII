% HW3_prob2_Noah

R = 1.1 ; 
J2 = 0.3 ; 
a = 6 ; 
e = 0.7 ; 

fun = @(r) r*(r^2 - 2*a*r + a*a*(1-e^2))-a*R^2*J2 ; 
rp = fzero(fun, a*(1-e)) ; 
ra = fzero(fun, a*(1+e)) ; 
r0 = fzero(fun, 0) ; 

state0 = [ 1 / rp ; 0 ] ; 

opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10 ) ; 

times = [] ; 
states = [] ; 
t_tmp = 0 ; 
state_tmp = state0 ; 

% for i = 1:16 
%     
%     [times_loop, states_loop] = ode45( @(t,y) state_dot_J2( t, y, R, J2, a, e ), [ t_tmp, t_tmp + 10 ], state_tmp, opts ) ; 
%     times = [ times ; times_loop ] ; 
%     states = [ states ; states_loop ] ; 
%     t_tmp = times(end) ; 
%     state_tmp = states(end,:) ; 
%     
% end 
[times, states] = ode45(@(theta, u) state_dot_J2(theta, u, R, J2, a, e), [ t_tmp, t_tmp + 10 ], state_tmp, options); 


k = sqrt( (ra-rp)*r0/(( ra-r0 )*rp ) ) ; 
kappa = 0.5 * sqrt( (rp*(ra-r0))/(ra*rp+r0*ra+r0*rp) ) ; 
[SN, CN, DN] = ellipj(kappa*times, k^2) ; 
r_analytic = r0 * rp * DN.^2./( r0 - rp*k*k*SN.^2 ) ; 

close all ; 
polarplot( times, 1./states(:,1) ) ; hold on ; grid on ; 
polarplot( times, r_analytic, '--' ) 



%% subfunctions 


function output = state_dot_J2( t, state, R, J2, a, e )  

u = state(1) ; 
udot = state(2) ; 

u_der = udot ; 
udot_der = 1/(a*(1-e^2)) + (3*(R^2)*J2*(u^2))/(2*a*(1-e^2)) - u ; 

output = [ u_der; udot_der] ; 

end 

