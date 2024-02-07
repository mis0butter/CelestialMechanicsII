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

%% subfunctions 

function KS = rv2KS( rv ) 

    % get position and velocity 
    r = rv(1:3) ;   x = r(1) ; y = r(2) ; z = r(3) ; 
    v = rv(4:6) ; 
    
    % deal with redundant variable 
    if r(1) >= 0 
        u1 = sqrt( 1/2 * ( x + norm(r) ) ) ; 
        u2 = y / ( 2 * u1 ) ; 
        u3 = z / ( 2 * u1 ) ; 
        u4 = 0 ;         
    else 
        u2 = sqrt( 1/2 * ( norm(r) - x ) ) ; 
        u1 = y / ( 2 * u2 ) ; 
        u3 = 0 ; 
        u4 = z / ( 2 * u2 ) ; 
    end 
    
    % KS parameters 
    u       = [ u1 ; u2 ; u3 ; u4 ] ; 
    uprime  = 1/2 * Lu_fn(u)' * [ v ; 0 ] ; 
    KS      = [ u ; uprime ] ; 

end 

function rv = KS2rv( KS )

    % KS in u basis 
    u       = KS(1:4) ;
    uprime  = KS(5:8) ; 
    
    % recover position 
    r = Lu_fn(u) * u ; 
    r = r(1:3) ; 

    % recover velocity 
    v = 2/norm(r) * Lu_fn(u) * uprime ; 
    v = v(1:3) ; 
    
    % Cartesian state 
    rv = [ r ; v ] ; 

end 

function Lu = Lu_fn(u) 

    Lu = [  u(1)    -u(2)   -u(3)   u(4) ; 
            u(2)    u(1)    -u(4)   -u(3) ; 
            u(3)    u(4)    u(1)    u(2) ; 
            u(4)    -u(3)   u(2)    -u(1) ] ; 

end 


