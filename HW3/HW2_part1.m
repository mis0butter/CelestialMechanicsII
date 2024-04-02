% HW 2 part 1 


%% 1.a: 

mu = 1.0 ; 


%% 1.b: Assume e < 1, discuss the character of the orbit for the c=0 and c>0 cases. Verify your
% conclusions using numerical integration in a Cartesian frame 

k = 5 ;  
h = 1 ;  
e = 0.2 ; 
theta = 1 ; 

% ode options 
toler   = 1e-8;         % 1e-14 accurate; 1e-6 coarse 
options = odeset('reltol', toler, 'abstol', toler ); 

% Circular orbit ICs 

r0 = [ 1*r_c ; 0 ; 0 ] ; 
v0 = [ 0 ; 1.1*v_c ; 0 ] ; 

% r 
alpha = sqrt( 1 + c / h^2 ) ; 
num   = h^2 / k ; 
den   = 1 / alpha^2 + e * cos( alpha * theta ) ;  
r     = num / den ; 

% c > 0 case 
c = 0.2 ; 

% integrate 
tspan = [ 0 T_c*250 ] ; 
[t,x] = ode45(@(t, rv) central_force(t, rv, k, c), tspan, [r0 v0], options); 

% plot 
figure 
    plot3( x(:,1), x(:,2), x(:,3) ) ; 
    hold on ; grid on ; 
    
% c = 0 case 
c = 0 ; 
tspan = [ 0 T_c*10  ] ; 
[t,x] = ode45(@(t, rv) central_force(t, rv, k, c), tspan, [r0 v0], options); 
    plot3( x(:,1), x(:,2), x(:,3) ) ; 
    scatter3( 0,0,0, 'filled', 'k' ) ; 
    legend( 'c > 0', 'c = 0' ) ; 
    title( 'Character of orbit for c = 0, c > 0 cases' ) 
    xlabel('x') ; ylabel('y') ; zlabel('z') ; 
    


%% 1.c: Assuming h=1, k=1, c=e=0.2, solve (via quadrature) for the time required from theta=0..1 

% treat r as a definite integral 
syms x
expr = x*log(1+x);
F = int(expr,[0 1])

syms theta 
c = 0.2 ; 
k = 1 ;  
h = 1 ;  
e = 0.2 ; 

alpha = sqrt( 1 + c / h^2 ) ; 

num   = h^2 / k ; 
den   = 1 / alpha^2 + e * cos( alpha * theta ) ;  

r     = num / den ; 
F     = int( r^2,[0 1] ) ; 
double(F) 


%% 1.d, e, f: circular orbit radius, velocity, period 

h = 1 ; 
k = 1 ;  
c = 0.2 ; 

r_c = ( h^2 + c ) / k  
v_c = sqrt( k/r_c - c/r_c^2 )  
T_c = 2 * pi * r_c^2 / h 


%% 1.f: stability 

u_c = 1 / r_c 
P_u = k * u_c^2 - c * u_c^3 ; 
dP_u = 2 * k * u_c - 3 * c * u_c^2 

term = 3 - u_c * dP_u / P_u 



%% subfunctions 

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


