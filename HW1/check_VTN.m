%% check that plot matches Dr. Russell's slides 

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
eps0    = 2 * dot( uprime0, uprime0 ) / dot( u0, u0 ) - mu / dot( u0, u0 ) ; 

% period in time (s) 
T = 2*pi*sqrt( a^3 / mu ) ; 
    
% period in tau (tau units) 
k = 1 ; 
T_tau = T / ( k * a ) ; 

state0  = [ u0 ; uprime0 ; t0 ; eps0 ] ; 

[ state_num, rv_num ] = prop_KS_num( state0, a_VTN, mu, T_tau * 10, N * 10, tol ) ; 
plot_rv( rv_num ) ; 