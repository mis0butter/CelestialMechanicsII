
     File: Readme.txt
     Description: Contains information about the F&G Stark and Kepler propagator
     Date: FEB 10, 2014
     Version: 2.0

 ***************************************************************************************
     This file is electronic supplementary material to the paper:
     Pellegrini, E., Russell, R.P., Vittaldev, V., "F and G Taylor Series
     Solutions to the Stark and Kepler Problems with Sundman Transformations", 
     Celestial Mechanics and Dynamical Astronomy, DOI: 10.1007/s10569-014-9538-7
     Contact: Etienne Pellegrini, etienne.pellegrini@utexas.edu

 ***************************************************************************************

     Copyright (C) 2014  Etienne Pellegrini and Ryan P. Russell

     This program is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with this program.  If not, see <http://www.gnu.org/licenses/>.

 ***************************************************************************************

     IMPORTANT: Please cite the scholarly article from where this code was obtained.
       Pellegrini, E., Russell, R.P., Vittaldev, V., "F and G Taylor Series
       Solutions to the Stark and Kepler Problems with Sundman Transformations", 
       Celestial Mechanics and Dynamical Astronomy, DOI: 10.1007/s10569-014-9538-7


     An up-to-date version of the code ( with any big fixes etc.) is available at:
       http://russell.ae.utexas.edu/index_files/fgstark.htm

     NOTE: Minor changes have been made to the code since producing the results
     presented in the paper. Efficiency has been improved at low
     orders, and divergence safeguards have been added. Moreover, the
     code used for the timings was more specialized for each case (in
     particular each case of the Sundman transformation).

 ***************************************************************************************

     Email : Etienne Pellegrini: etienne.pellegrini@utexas.edu
             Ryan P. Russell: ryan.russell@utexas.edu

 ***************************************************************************************
     *** The F&G Stark Series propagator is organised as follows:

       - main_fg.f90 is an example PROGRAM that can be run to test the propagator
      
       - FGStark_propagator.for is the file containing the driver subroutine as
         well as the subroutines computing the Series coefficients.

     *** Driver subroutine FGStark_oneStep_propagator:

       - Used to integrate position, velocity and time over one tau-step (has to
         be wrapped in another routine in order to propagate for longer periods).

       - Also computes different accuracy measurements: the delta Hamiltonian as
         well as the norm of the delta in position, velocity and time resulting from
         the last order taken by the series. The Hamiltonian can not be
         lower than machine precision (1e-16), whereas the other three
         scalars can be lower. The velocity error will usually be one or two
         orders of magnitude larger than the position and time errors.
         The accuracy measurements are controled with accuracyFlag, an
         array or booleans. When an element of accuracyFlag is .TRUE.
         the corresponding accuracy indicator is computed. The cost (in
         % of the compute time) of each of those computations is given
         for order 12.
           - 1: Hamiltonian  (2%)
           - 2: delR         (1.5%)
           - 3: deV          (3%)
           - 4: delT         (0%)
         As can be observed, the time error estimate is essentially
         free, and its use as a rough estimate of the accuracy is suggested.

       - Depending on a flag set by the user, chooses to use the
         Stark or the Kepler coefficients

       - If activated with divSafeguard = .TRUE., the program tests the series for divergence.
         If the statusFlag returns a value of -1, the series most likely diverged, and the user
         should reduce the tau-step. The divergence safeguard is an expensive computation and 
         should only be turned on if divergence is possible

       - Inputs:
         - Initial state
         - Perturbation
         - Delta tau
         - Desired order
         - Sundman parameters (see below)
         - caseStark: indicates if a Stark propagation is to be computed
         - divSafeguard: boolean indicating if the user wants to check
         for divergence
         - Accuracy flag, an array of booleans indicating which
         accuracy measurements are to be computed.
      
       - Outputs:
         - Final state
         - Accuracy measurements
         - statusFlag: -1 indicates divergence

     *** The Sundman transformation is described by: dt = cSun * r**alphaSun * dtau 
       In order to choose which Sundman transformation to use, the user selects a
       case ('caseSun' variable): 
       - caseSun = 0 => alphaSun = 0, independent variable is proportional to
         time
       - caseSun = 1 => alphaSun = 1, independent variable is proportional to
         eccentric anomaly
       - caseSun = 2 => alphaSun = 2, independent variable is proportional to
         true anomaly
       - caseSun = 3 => alphaSun = 3/2, independent variable is proportional to
         intermediate anomaly
       - caseSun = 4 => alphaSun = user provided. It can be any double. CAREFUL:
         when used for a value of alpha cited above (0, 1, 2, 3/2), the generic
         case is less efficient than the corresponding specialized case, and the
         specialized cases should be preferred.

     *** Units: Position (LU), velocity (LU/TU), acceleration (LU/TU^2). 
       The outputs' units are consistent. The standard gravitational
       parameter mu is taken to be 1 LU^3/TU^2


 ***************************************************************************************
  Note on the computational efficiency:

     The code has been changed since producing the results
     presented in the paper. Efficiency has been improved at low
     orders, and divergence safeguards have been added. Moreover, the
     code used for the timings was more specialized for each case (in
     particular each case of the Sundman transformation).


     In the case of a particular propagation (particular Sundman
     transformation, only Kepler, etc...) the user is strongly
     encouraged in stripping the unused parts of the code in order to
     improve computational efficiency.

 ***************************************************************************************
  *** WARNING: The large size of the file can cause problems with the
    debug mode of Microsoft Visual Studio. The user should be able to
    step over the FG Stark propagation, but stepping inside this file
    might crash Visual Studio. If the user wants to step in those
    routines, he/she might have to strip the file and make it smaller.

  *** WARNING: The FandG_KeplerStarkSundman file is large and might take
    several dozens of minutes to compile. Compiling it as a library
    might help, since only one compilation would be necessary

 ***************************************************************************************
