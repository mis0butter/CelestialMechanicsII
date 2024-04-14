%% routine for FGStark_oneStep_propagator in MATLAB 

X0   = [ 0.2d0, 0.d0, 0.d0, 0.d0, 3.0d0, 0.d0, 0.d0 ] ; 

acc  = [ 1.d-1, 1.d-1, 1.d-1 ] ; 
dtau = 6.283185307179589D-1 ; 

cSun     = 1.0d0 ; 
caseSun  = 0 ; 
alphaSun = 2.d0 ; 

caseStark = true ; 
divSafeguard = true ; 

order = 12 ; 
accuracyFlag = [ true, true, true, true ] ; 

% initialize outputs ; 
Xf          = X0 ; 
accuracyOut = zeros(4,1) ; 
statusFlag  = 0 ; 

[ Xf, accuracyOut, statusFlag ] = FGStark_oneStep_propagator( X0, acc, dtau, order, ... 
    cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
    accuracyFlag, Xf, accuracyOut, statusFlag ) ; 


%% subfunctions 

function [ Xf, accuracyOut, statusFlag ] = FGStark_oneStep_propagator( X0, acc, dtau, order, ... 
    cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
    accuracyFlag, Xf, accuracyOut, statusFlag )

MAXORD = 30 ; 

%         !-- Declaration of variables
%         ! I/O
%         real*8, intent(in)  :: X0(7)          ! X0: Initial Conditions. Dimension 7x1, containing r0 (3D), v0 (3D), t0 
%         real*8, intent(in)  :: acc(3), dtau   ! acc: Stark acceleration (3D), dtau: delta tau for integration
%         real*8, intent(in)  :: cSun           ! Sundman transformation constant (dt = cSun * r**alphaSun * dtau). Typically cSun = 1
%         integer, intent(in) :: caseSun        ! Select a case for the Sundman transformation: 0: alpha = 0, tau = time
%                                               !     1: alpha = 1, tau = eccentric anomaly ; This is the case resulting in best efficiency of the series
%                                               !     2: alpha = 2, tau = true anomaly
%                                               !     3: alpha = 3/2, tau = intermediate anomaly
%                                               !     4: alpha = alphaSun, generic Sundman transformation: user provided alpha
%         integer, intent(in) :: order          ! order: order of the desired Taylor Series integration
%         real*8, intent(in)  :: alphaSun       ! Value of alpha when case 4 is selected. For other cases, alphaSun is not used
%         logical, intent(in) :: caseStark      ! Indicates whether a Stark or Kepler propagation should be computed
%         logical, intent(in) :: divSafeguard   ! Indicates whether divergence tests should be realized
%         logical, intent(in) :: accuracyFlag(4)! Indicates which accuracy measurements are to be computed (1 = Hamiltonian, 2 = delR, 3 = delV, 4 = delT)
% 
%         real*8, intent(out) :: Xf(7)          ! Xf: Output vector, dimension 7x1
%         real*8, intent(out) :: accuracyOut(4) ! Accuracy measurements: change in Hamiltonian, last term error estimates for r, v, p
%         integer, intent(out) :: statusFlag    ! Indicates possible overflow/divergence of the series. If statusFlag = -1, the routine exited before finishing the
%                                               ! computation. The step size has to be reduced.


%         ! Internal variables
%         integer, parameter  :: rng = RANGE(1.d0)
%         real*8, parameter  :: frst = 10.d0**(rng-10)
% 
%         real*8  :: dtaun, dtaunprime, normR, v2
%         real*8  :: hamil, hamilPlus1
%         real*8  :: delRvec(3), delVvec(3)
%         real*8  :: Ftot, Gtot, Htot, Ttot, dtaunprimej 
%         real*8  :: Ftotprime, Gtotprime, Htotprime
%         real*8  :: first(4), new(4), newprime(3), numerator
%         integer :: j, i, k, set(4)
%         real*8  :: FGHT(4*MAXORD), coeff(4)

% internal variables 
rng  = 307 ; 
frst = 10.d0^(rng-10) ; 

dtaun = 0 ; dtaunprime = 0 ; normR = 0 ; v2 = 0 ; 
hamil = 0 ; hamilPlus1 = 0 ; 
delRvec = zeros(3,1) ; delVvec = zeros(3,1) ; 
Ftot  = 0 ; Gtot = 0 ; Htot = 0 ; Ttot = 0 ; dtaunprimej = 0 ; 
Ftotprime = 0 ; Gtotprime = 0 ; Htotprime = 0 ; 
first = zeros(4,1) ; new = zeros(4,1) ; newprime = zeros(3,1) ; numerator = 0 ; 
j = 0 ; i = 0 ; k = 0 ; set = zeros(4,1) ; 
FGHT = zeros( 4 * MAXORD, 1 ) ; coeff = zeros(4,1) ; 

%% Begin subroutine

% Variables needed for overflow safeguard
% first = frst ; 
% set   = 0 ; 
statusFlag  = 0 ; 
accuracyOut = 0.d0 ; 

% Initialize the sums
Ftot = 0.d0 ; 
Gtot = 0.d0 ; 
Htot = 0.d0 ; 
Ttot = 0.d0 ; 
Ftotprime = 0.d0 ; 
Gtotprime = 0.d0 ; 
Htotprime = 0.d0 ; 

dtaun = 1.d0 ; 
dtaunprime = 1.d0 ; 

% Compute the Hamiltonian
normR = sqrt( X0(1)*X0(1) + X0(2)*X0(2) + X0(3)*X0(3) ) ; 
v2    = X0(4)*X0(4) + X0(5)*X0(5) + X0(6)*X0(6) ; 
hamil = 0.5d0 * v2 - 1.d0/normR - dot(X0(1:3), acc) ; 
        
% Call the Maple generated file to get the values of
% the F&G Stark series coefficients. Call depending on the
% caseSun and caseStark

switch caseSun

    case (0) % The independent variable is proportional to time, dt = dtau 
        FGHT = starkCoeffs_time(X0, acc, order, cSun, FGHT) ; 

    case (1) % The independent variable is proportional to the eccentric anomaly
        FGHT = starkCoeffs_ecc(X0, acc, order, cSun, FGHT) ; 

    case (2) % The independent variable is proportional to the true anomaly
        FGHT = starkCoeffs_true(X0, acc, order, cSun, FGHT) ; 

    case (3) % The independent variable is proportional to the intermediate anomaly
        FGHT = starkCoeffs_inter(X0, acc, order, cSun, FGHT) ; 

    case (4) % Generic Sundman transformation: dt = cSun * r**alphaSun * dtau
        FGHT = starkCoeffs_generic(X0, acc, order, cSun, alphaSun, FGHT) ; 

    case default
        display('ERROR: Please select a Sundman case between 0 and 4') ; 

end  
    

% Form the Taylor Series from the coefficients
k = 1 ; 
for j = 1 : order

    dtaun = dtaun*dtau ; 

    % Compute position TS
    coeff(1) = FGHT(k) ; 
    coeff(2) = FGHT(k+1) ; 
    coeff(3) = FGHT(k+2) ; 
    coeff(4) = FGHT(k+3) ; 
    k = k + 4 ; 

    new(1) = coeff(1)*dtaun ; 
    new(2) = coeff(2)*dtaun ; 
    new(3) = coeff(3)*dtaun ; 
    new(4) = coeff(4)*dtaun ; 


    if (divSafeguard) 
        
        % Heuristically check convergence of the series. If the new
        % term in the series is bigger than 100* the very first term, we
        % exit the integrator
        for i = 1 : 4

          % Set the "first" value (first non zero value taken by each
          % series)
          if ( (set(i) == 0) && ( abs(new(i)) > 0.d0 ) ) 
            first(i) = 100.d0*(1.d0 + abs(new(i))) ; 
            set(i) = 1 ; 
          end 

          % Compare the new value to "first"
          if (abs(new(i)) > first(i)) 
              statusFlag = -1 ; 
              break 
          end 
          
        end 
    end 

    Ftot = Ftot + new(1) ; 
    Gtot = Gtot + new(2) ; 
    Htot = Htot + new(3) ; 
    Ttot = Ttot + new(4) ; 

    % Compute velocity TS: derivative of the position.
    dtaunprimej = dtaunprime*j ; 

    newprime(1) = coeff(1)*dtaunprimej ; 
    newprime(2) = coeff(2)*dtaunprimej ; 
    newprime(3) = coeff(3)*dtaunprimej ; 

    Ftotprime = Ftotprime + newprime(1) ; 
    Gtotprime = Gtotprime + newprime(2) ; 
    Htotprime = Htotprime + newprime(3) ; 

    dtaunprime = dtaunprime*dtau ; 

end 
% 
% 
%         ! Finish the TS integration: multiply by the basis vectors
% 
%         ! r = F*r0 + G*v0 + H*thrust
%         Xf(1:3) = X0(1:3) + Ftot*X0(1:3) + Gtot*X0(4:6) +
%      *              Htot * acc
% 
%         ! v = Fdot*r0 + Gdot*v0 + Hdot*thrust
%         Xf(4:6) = Ftotprime * X0(1:3) + Gtotprime * X0(4:6) +
%      *              Htotprime * acc
% 
%         ! Careful: to go from r to v, differentiation wrt 
%         ! time. However, the coefficients have been
%         ! differentiated with respect to tau. v = dr/dt =
%         ! dr/dtau*dtau/dt = dr/dtau * 1/(cSun * r**alphaSun)
%         normR = dsqrt( Xf(1)*Xf(1) + Xf(2)*Xf(2) + Xf(3)*Xf(3) )
%         select case(caseSun)
%           case(0)
%             numerator = (cSun)
%           case(1)
%             numerator = (cSun*normR)
%           case(2)
%             numerator = (cSun*normR**2)
%           case(3)
%             numerator = (cSun*sqrt(normR**3))
%           case(4)
%             numerator = (cSun*normR**alphaSun)
%         end select
%         
%         Xf(4:6) = Xf(4:6) / numerator
% 
%         ! time = time0 + dtime
%         Xf(7) = X0(7) + Ttot
%         
%         ! Accuracy computations:
%         if (accuracyFlag(1)) then
%           ! Compute new Hamiltonian
%           v2 = Xf(4)*Xf(4) + Xf(5)*Xf(5) + Xf(6)*Xf(6)
%           hamilPlus1 = 0.5d0 * v2 - 1.d0 / normR -dot_product(Xf(1:3),
%      *      acc)
%           accuracyOut(1) = dabs(hamilPlus1 - hamil)
%         endif
%         if (accuracyFlag(2)) then
%           delRvec = new(1) * X0(1:3) + new(2) * X0(4:6) + new(3) * acc
%           accuracyOut(2) = dsqrt( delRvec(1)*delRvec(1) + 
%      *                    delRvec(2)*delRvec(2) + delRvec(3)*delRvec(3))
%         endif
%         if (accuracyFlag(3)) then 
%           delVvec = ( newprime(1) * X0(1:3) + newprime(2) * X0(4:6) + 
%      *              newprime(3) * acc ) / numerator
%           accuracyOut(3) = dsqrt( delVvec(1)*delVvec(1) + 
%      *                   delVvec(2)*delVvec(2) + delVvec(3)*delVvec(3) )
%         endif
%         if (accuracyFlag(4)) then
%           accuracyOut(4) = dabs(new(4))
%         endif


end 

function FGHT = starkCoeffs_time(X, acc, order, cSun, FGHT) 

%         use maxOrder
%         implicit double precision (t)
% 
%         integer, intent(in) :: order
%         real*8, intent(in)  :: X(7), acc(3), cSun
% 
%         real*8  :: FGHT(MAXORD*4)
%         real*8  :: r, RdotV, PdotR, PdotV, v2, p2, k

    r  = sqrt( X(1)*X(1) + X(2)*X(2) + X(3)*X(3) ) ; 
    v2 = X(4)*X(4) + X(5)*X(5) + X(6)*X(6) ; 
    
    RdotV = dot(X(1:3),X(4:6)) ; 
    PdotR = dot(X(1:3),acc) ; 
    PdotV = dot(X(4:6),acc) ; 
    
    p2 = acc(1)*acc(1) + acc(2)*acc(2) + acc(3)*acc(3) ; 
    k  = cSun ; 

    switch (order)

        case (0)

        case (1)
            
            FGHT(1) = 0.0D0 ; 
            FGHT(2) = k ; 
            FGHT(3) = 0.0D0 ; 
            FGHT(4) = k ; 
            
        case (2) 

            t1 = k^2 ; 
            t2 = r^2 ; 
            FGHT(1) = 0.0D0 ; 
            FGHT(2) = k ; 
            FGHT(3) = 0.0D0 ; 
            FGHT(4) = k ; 
            FGHT(5) = -0.5000000000000000D0 * t1 / t2 / r ; 
            FGHT(6) = 0.0D0 ; 
            FGHT(7) = 0.5000000000000000D0 * t1 ; 
            FGHT(8) = 0.0D0 ; 
      
    end 

end 


