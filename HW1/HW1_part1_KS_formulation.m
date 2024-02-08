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
% tau = 2 

% KS_time = KS_time_fn( u0, uprime0, tau, eps ) ;  

% period 
T = 2*pi*sqrt( a^3 / mu ) ; 

t_hist   = [] ; 
tau_hist = [] ; 
for tau = 1 : 100 
    
    KS_time = KS_time_fn( u0, uprime0, tau, eps ) ; 
    t = KS_time(end) ; 
    
    t_hist   = [ t_hist ; t ] ; 
    tau_hist = [ tau_hist ; tau ] ; 
    
end 

figure(1) ; hold on ; grid on ; 
    plot( t_hist ) ; 
    yline( T ) ; 
    
%% part d: 

toler   = 1e-8;         % 1e-14 accurate; 1e-6 coarse 
options = odeset('reltol', toler, 'abstol', toler ); 
a = zeros(4,1) ; 

% Solve ODE 
KS_t_eps_0 = [ KS0 ; t0 ; eps0 ] ; 
tau_span = [ 0 6 ] ; 
[tau, x] = ode45(@(tau, KS_t_eps) KS_EOM( tau, KS_t_eps, a, mu), tau_span, KS_t_eps_0, options) ; 

    
%% 

test_fn = @(x) KS_time_fn( u0, uprime0, x, eps) ; 
X = fzero( test_fn, 6 ) 


%% subfunctions 

function KS_t_eps_prime = KS_EOM( tau, KS_t_eps, a, mu ) 

    KS_t_eps_prime = zeros(10, 1) ; 
    
    u        = KS_t_eps(1:4) ; 
    u_prime  = KS_t_eps(5:8) ; 
    
    KS_t_eps_prime(1:4) = u_prime ; 
    
    r = dot( u, u ) ; 
    eps = 2 * dot( u_prime, u_prime ) / r - mu / r ; 
    
    u_prime_prime = eps/2 * u + 1/2 * r * Lu_fn(u)' * a ; 
    
    t_prime = r ; 
    
    eps_prime = 2 * dot(Lu_fn(u)' * a, u_prime ) ; 
    
    KS_t_eps_prime = [ u_prime ; u_prime_prime ; t_prime ; eps_prime ] ; 

end 

function KS_zero = KS_zero_fn( u0, uprime0, tau, eps ) 

    KS_time = KS_time_fn( u0, uprime0, tau, eps ) ; 
    t = KS_time(end) ; 


end 





