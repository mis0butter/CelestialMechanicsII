
z = pi/4 

Cz = C_stumpff(z) 
Sz = S_stumpff(z) 

%% Newton Raphson: find universal variable x 

mu = 1 ; 

% initial params 
r0 = [ 1; 0; 0 ] ;          r0_norm = norm(r0) ; 
v0 = [ 0; 0.01 ; 0 ] ;      v0_norm = norm(v0) ; 
t_given = 1 ; 

% rv to oe 
oe0 = rv2oe( [ r0 ; v0 ], mu ) ; 
a   = oe0(1) ; 

% first guess 
x_n   = sqrt(mu) * abs( 1/a ) * t_given ; 
f     = t_fn( x_n, mu, r0, v0 ) - t_given ; 
dfdx  = dtdx_fn( x_n, mu, r0, v0 ) ; 
x_np1 = x_n - f / dfdx ; 

% iterate 
tol = 1e-5 ; 
k   = 0 ; 
while abs(x_np1 - x_n) > tol 
    
    k = k + 1 ; 
    if k > 10000 
        display( 'too much' ) ; 
        break 
    end 

    x_n   = x_np1 ; 
    f     = t_fn( x_n, mu, r0, v0 ) - t_given ; 
    dfdx  = dtdx_fn( x_n, mu, r0, v0 ) ; 
    x_np1 = x_n - f / dfdx ; 
    
end 

% set values 
x     = x_np1 ; 
alpha = 2 / r0_norm - v0_norm^2 / mu ; 
z     = x^2 * alpha ; 

S = S_stumpff(z) ; 
C = C_stumpff(z) ; 


%% find f and g 

f = 1 - x^2 / r0_norm * C ; 
g = t_given - x^3 / sqrt(mu) * S ; 

% r_norm = sqrt(mu) * dtdx_fn( x, mu, r0, v0 ) ; 
r_norm = norm( f * r0 + g * v0 ) ; 

fdot = sqrt(mu) / ( r0_norm * r_norm ) * x * ( z * S - 1 ) ; 
gdot = 1 - x^2 / r_norm * C ; 

r = f * r0 + g * v0 ; 
v = fdot * r0 + gdot * v0 ; 


%%  subfunctions 

function t_out = t_fn( x, mu, r0, v0 ) 

    r0_norm = norm(r0) ; 
    v0_norm = norm(v0) ; 
    
    alpha = 2 / r0_norm - v0_norm^2 / mu ; 
    z     = x^2 * alpha ; 
    
    S = S_stumpff(z) ; 
    C = C_stumpff(z) ; 
    
    t_out = x^3 / sqrt(mu) * S + ... 
            dot(r0, v0) / mu * x^2 * C + ... 
            r0_norm * x / sqrt(mu) * ( 1 - z * S ) ; 
    
end


function dtdx_out = dtdx_fn( x, mu, r0, v0 ) 

    r0_norm = norm(r0) ; 
    v0_norm = norm(v0) ; 
    
    alpha = 2 / r0_norm - v0_norm^2 / mu ; 
    z     = x^2 * alpha ; 
    
    S = S_stumpff(z) ; 
    C = C_stumpff(z) ; 
    
    dtdx_out = x^2 / sqrt(mu) * C + ... 
               dot(r0, v0) / mu * x * ( 1 - z * S ) + ... 
               r0_norm / sqrt(mu) * ( 1 - z * C ) ; 

end 


function Cz = C_stumpff( z )
    
    if z > 0 
        
        num = 1 - cos( sqrt(z) ) ;         
        den = z ; 
        
    elseif z < 0 

        num = 1 - cosh( sqrt(-z) ) ; 
        den = z ; 
        
    else 
        
        num = 1 ; 
        den = 2 ; 
        
    end 
    
    Cz = num / den ; 

end 

function Sz = S_stumpff( z ) 

    if z > 0 
        
        num = sqrt( z ) - sin( sqrt(z) ) ; 
        den = sqrt( z^3 ) ; 
    
    elseif z < 0 
        
        num = sinh( sqrt(-z) ) - sqrt(-z) ; 
        den = sqrt( (-z)^3 ) ; 
        
    else 
        
        num = 1 ; 
        den = 6 ; 
        
    end 
    
    Sz = num / den ; 

end 

