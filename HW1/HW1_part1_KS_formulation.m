%% PART 1: KS FORMULATION 

init 

%% PART 1.1: UNPERTURBED 

% part a: Code the transformation subroutines for transformations in both 
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

%% { u, u' } <--> { r, v } 

% going from u to r (straightforward) 

KS0 = rv2KS( rv0 ) ; 
rv  = KS2rv( KS0 ) ; 





