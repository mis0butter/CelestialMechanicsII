function [tout, yout] = ode78rpr(FUNCTION, t0, tfinal, dtmin, y0, tol)

%1-30-2013, RPR (ryan.russell@utexas.edu) modified to allow for fixed step from existing rk78 found on web, and
%other minor changes...

%variable step rk7(8) integrator stepping from t0 to tf
%where dtmin is the required first step as input
%if you want fixed step set tol=0.0 and then it will step according to
%dtmin.  gives back the time and state output at each time step.
%the ydot function should be of the form  f = FUNCTION(t,y)
%tol should very from 1e-17 fine to 1e-4 course (or 0.0 for fixed step
%it works in backwards time also if tf<t0

% The Fehlberg coefficients:
% From Matlab website, 1996
global numcalls

%rpr, made these global so you didnt have to keep recomputing them, renamed
%them to avoid interfering with other variables in the calling routines
global alphaGL betaGL chiGL psiGL powGL firstcall 


if(isempty(firstcall))
firstcall=1;
alphaGL = [ 2./27. 1/9 1/6 5/12 .5 5/6 1/6 2/3 1/3 1 0 1 ]';
betaGL = [ [  2/27  0  0   0   0  0  0  0  0  0  0   0  0  ]
    [  1/36 1/12  0  0  0  0  0  0   0  0  0  0  0  ]
    [  1/24  0  1/8  0  0  0  0  0  0  0  0  0  0 ]
    [  5/12  0  -25/16  25/16  0  0  0  0  0  0   0  0  0  ]
    [ .05   0  0  .25  .2  0  0  0  0  0  0  0  0 ]
    [ -25/108  0  0  125/108  -65/27  125/54  0  0  0  0  0  0   0  ]
    [ 31/300  0  0  0  61/225  -2/9  13/900  0  0  0   0  0  0  ]
    [ 2  0  0  -53/6  704/45  -107/9  67/90  3  0  0  0  0  0  ]
    [ -91/108  0  0  23/108  -976/135  311/54  -19/60  17/6  -1/12  0  0  0  0 ]
    [2383/4100 0 0 -341/164 4496/1025 -301/82 2133/4100 45/82 45/164 18/41 0 0 0]
    [ 3/205  0   0  0   0    -6/41  -3/205   -3/41     3/41   6/41   0   0  0 ]
    [-1777/4100 0 0 -341/164 4496/1025 -289/82 2193/4100 ...
    51/82 33/164 12/41 0 1 0]...
    ]';
chiGL = [ 0 0 0 0 0 34/105 9/35 9/35 9/280 9/280 0 41/840 41/840]';
psiGL = [1  0  0  0  0  0  0  0  0  0  1 -1  -1 ]';
powGL = 1/8;
%pause
end

if(tol==0)
    varstep=false;
else
    varstep=true;
end
last=false;

% Initialization
t = t0;
hmax = (tfinal - t)/1.0;  %/2.5
hmin = dtmin;

hmin = abs(hmin); %rpr
hmax = abs(hmax); %rpr

h = sign(tfinal - t)*dtmin;
y = y0(:);
f = y*zeros(1,13);
tout = t;
yout = y.';
tau = tol * max(norm(y, 'inf'), 1);


forward=1;  %rpr
if (h<0), forward=-1; end    %rpr

count=0;
% The main loop
while (t*forward < tfinal*forward) %& (abs(h) >= hmin)  %rprmod
    count=count+1;
       %[ith t tfinal]
    
    if (t + h)*forward > tfinal*forward 
        h = tfinal - t;
        last=true;
    end
%h;
    % Compute the slopes
    f(:,1) = feval(FUNCTION,t,y);
    for j = 1: 12
        
        %MODE RYAN 7-20-09  this is faster way to calculate not wasting
        %zero multiplications
        %yint=y+h*f*betaGL(:,j);
        yint=y+h*f(:,1:j)*betaGL(1:j,j);
        
        f(:,j+1) = feval(FUNCTION, t+alphaGL(j)*h, yint);
    end

    % Truncation error term
    if(varstep)
        gamma1 = h*41/840*f*psiGL;

        % Estimate the error and the acceptable error
        delta = norm(gamma1,'inf');
        tau = tol*max(norm(y,'inf'),1.0);
    end
    % Update the solution only if the error is acceptable
    if varstep==false   |  delta <= tau  | abs(h)==hmin     
        t = t + h; 
        y = y + h*f*chiGL;
        tout = [tout; t];
        yout = [yout; y.'];
    end

    % Update the step size
    if varstep  &  delta ~= 0.0
        h = forward*min(hmax, abs(0.8*h*(tau/delta)^powGL)); 
        if last==false
            h = forward*max(abs(h),hmin);
        end

    end
end;
count;
if (t*forward < tfinal*forward)
    disp('SINGULARITY LIKELY.')
    t
    tfinal
    h
    hmin
end

