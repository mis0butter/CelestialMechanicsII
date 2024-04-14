
!     File: main_fg.f90
!     Description: PROGRAM for testing purposes
!     Date: Feb 10, 2014
!     Version: 1.10
!
! ***************************************************************************************
!     This file is electronic supplementary material to the paper:
!     Pellegrini, E., Russell, R.P., Vittaldev, V., "F and G Taylor Series
!     Solutions to the Stark and Kepler Problems with Sundman Transformations", 
!     Celestial Mechanics and Dynamical Astronomy, DOI: 10.1007/s10569-014-9538-7
!     Contact: Etienne Pellegrini, etienne.pellegrini@utexas.edu
!
! ***************************************************************************************
!
!     Copyright (C) 2014  Etienne Pellegrini and Ryan P. Russell
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! ***************************************************************************************
!
!     IMPORTANT: Please cite the scholarly article from where this code was obtained.
!       Pellegrini, E., Russell, R.P., Vittaldev, V., "F and G Taylor Series
!       Solutions to the Stark and Kepler Problems with Sundman Transformations", 
!       Celestial Mechanics and Dynamical Astronomy, DOI: 10.1007/s10569-014-9538-7
!
!
!     An up-to-date version of the code ( with any big fixes etc.) is available at:
!       http://russell.ae.utexas.edu/index_files/fgstark.htm
!
!     NOTE: Minor changes have been made to the code since producing the results
!     presented in the paper. Efficiency has been improved at low
!     orders, and divergence safeguards have been added. Moreover, the
!     code used for the timings was more specialized for each case (in
!     particular each case of the Sundman transformation).
!
! ***************************************************************************************
!
!     Email : Etienne Pellegrini: etienne.pellegrini@utexas.edu
!             Ryan P. Russell: ryan.russell@utexas.edu
!
! ***************************************************************************************


      PROGRAM mainFG


        implicit none

        !-- Declaration of variables
        ! I/O
        real*8  :: X0(7), acc(3), dtau    ! X0: Initial Conditions, acc: Stark acceleration (3D), dtau: delta tau for integration
        real*8  :: cSun, alphaSun         ! Sundman transformation parameters (dt = cSun * r**alphaSun * dtau)
        integer :: caseSun                ! Select a case for the Sundman transformation: 0: alpha = 0, tau = time
                                          !     1: alpha = 1, tau = eccentric anomaly
                                          !     2: alpha = 2, tau = true anomaly
                                          !     3: alpha = 3/2, tau = intermediate anomaly
                                          !     4: alpha = alphaSun, generic Sundman transformation: user provided alpha
        integer :: order                  ! order: order of the TS integration
        real*8  :: Xf(7)                  ! Xf: Integrated vector
        real*8  :: accuracyOut(4)         ! Accuracy measurements
        logical :: caseStark              ! Indicates whether the user wants a Stark or Kepler propagation
        logical :: divSafeguard           ! Indicates whether the user wants to perform divergence checks (expensive)
        integer :: statusFlag             ! Indicates possible divergence of the series. 
                                          ! If statusFlag = -1, the series likely diverged
        logical :: accuracyFlag(4)        ! Indicates which accuracy measurements are to be computed 
                                          ! (1 = Hamiltonian, 2 = delR, 3 = delV, 4 = delT)

        !-- Begin subroutine


        X0 = (/0.2d0, 0.d0, 0.d0, 0.d0, 3.0d0, 0.d0, 0.d0/)
        acc = (/1.d-1, 1.d-1, 1.d-1/)
        dtau = 6.283185307179589D-1

        cSun = 1.0d0
        caseSun = 2
        alphaSun = 2.d0

        caseStark = .TRUE.
        divSafeguard = .FALSE.

        order = 12
        accuracyFlag = (/ .TRUE., .TRUE., .TRUE., .TRUE./)

        print*, '*** Starting the F&G Stark series propagation ***'
        print*, ''

        print*, 'X0: ', X0
        print*, ''

        print*, 'Case: ', caseSun, 'cSun: ', cSun, 'alphaSun: ',&
          alphaSun, 'order: ', order
        print*, ''

        if (.NOT.caseStark) then
          acc = 0
        endif

        CALL FGStark_oneStep_propagator(X0, acc, dtau, order, &
        cSun, caseSun, alphaSun, caseStark, divSafeguard, accuracyFlag,&
        Xf, accuracyOut, statusFlag)
        
        print*, '================================================='
        print*, ''
       
        print*, 'Status: ', statusFlag
        print*, 'Xf: ', Xf
        print*,''
        print*, 'del Hamiltonian: ', accuracyOut(1)
        print*, 'del norm R:      ', accuracyOut(2)
        print*, 'del norm V:      ', accuracyOut(3)
        print*, 'del time:        ', accuracyOut(4)
        print*, '*** Propagation ended ***'

      end PROGRAM
