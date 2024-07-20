function                                                                     ...
[                                                                            ...
  actualSimulationTimes,                                                     ...
  scalarFirstUnitQuaternions,                                                ...
  angularVelocitiesBodyCoords                                                ...
] =                                                                          ...
generateSolutionsForFreeRotationRigidBodyKinematics                          ...
        (                                                                    ...
          varargin                                                           ...
        )
%===============================================================================
%|
%|  FUNCTION:
%|
%|    generateSolutionsForFreeRotationRigidBodyKinematics
%|
%|-----------------------------------------------------------------------------
%|
%|  PURPOSE:
%|
%|    Generate solutions at specified times for torque free rotation of
%|    rigid bodies.
%|
%|-----------------------------------------------------------------------------
%|
%|  INPUTS:
%|
%|    J
%|      A three by three matrix of double precision floating point Inertia
%|      Tensor values.
%|      UNIT(s):  [kilograms / meter squared]
%|
%|    initialScalarFirstUnitQuaternion
%|      A one dimensional vector of four elements containing the initial
%|      scalar first unit quaternion.
%|      UNIT(s):  [nondimensional]
%|
%|    intitalAngularVelocityBodyCoords
%|      A one dimensional vector of length three containing the initial
%|      angular velocity vector in non-inertial body coordinates.
%|      UNIT(s):  [radians/second]
%|
%|    numberTimeSamples
%|      Number of time sample values generated for the results.
%|      UNIT(s):  [nondimensional]
%|
%|    simulationEndTime
%|      Maximum time value for the simulation.
%|      UNIT(s):  [seconds[
%|
%|-----------------------------------------------------------------------------
%|
%|  OUTPUTS:
%|
%|    actualSimulationTimes
%|      A one dimensional vector the ith element of which is the actual ith
%|      simulation time value.
%|      UNIT(s):  [seconds]
%|
%|    scalarFirstUnitQuaternions
%|      A two dimensional matrix having four rows and as many columns as
%|      there are time samples.
%|      The jth column of this matrix is the scalar first unit quaternion
%|      for the jth sample time.
%|      UNIT(s):  [nondimensional]
%|
%|    angularVelocitiesBodyCoords
%|      A two dimensional vector the ith element of which is the
%|      angular velocity in non-inertial body coordinates at the
%|      ith sample time.
%|      UNIT(s):  [radians/second]
%|
%|-----------------------------------------------------------------------------
%|
%|  NOTE(s):
%|
%|    [ 1 ]  Rotating object kinematics is "scalar first unit quaternion"
%|           based.
%|
%|    [ 2 ]  Scalar First Unit Quaternion Definition:
%|           Use Equation (22) on Page 264 of Reference [ 1 ]
%|           EXCEPT:
%|             [ 1 ]  Equation (22) is scalar LAST.
%|             [ 2 ]  Here scalar FIRST is used.
%|
%|    [ 3 ]  Ordinary Differential Equations (ODEs) are solved using the
%|           Matlab "ode45" function to apply the Dorman-Prince 4(5)
%|           explicit embedded variable time step Runge-Kutta method.
%|
%|    [ 4 ] State Vector Components:
%|
%|          [ 4.1 ]  State Vector components 1, 2, 3 and 4 are the
%|                   "scalar first unit quaternion".
%|
%|          [ 4.2 ]  State Vector components 5, 6 and 7 are the
%|                   angular velocity vector.
%|
%|    [ 5 ] Scalar First Unit Quaternion
%|
%|          Use Equation (22) on Page 264 of Reference [ 1 ]
%|          EXCEPT:
%|            [ 5.1 ]  Equation (22) is scalar LAST.
%|            [ 5.2 ]  Here scalar FIRST is used.
%|
%|
%|    [ 6 ] Time rate of change of the "scalar first unit quaternion":
%|
%|          Use the combination of the following two items:
%|
%|          [ 6.1 ] The FIRST PART of Equation (46) on Page 266
%|                  of Reference [ 1 ].
%|                  This states that the time rate of change of the
%|                  attitude quaterion is equal to the cross product
%|                  of two quaternions"
%|                  [ 6.1.1 ] The pure quaternion container of the
%|                            angular velocity 3D vector.
%|                  [ 6.1.2 ] The attitude quaternion.
%|
%|          [ 6.2 ] Equation (305) on Page 482 of Reference [ 2 ].
%|                  This equation shows how to use the quaternion
%|                  product when the quaternion cross product is
%|                  specified.
%|                  Here the angular velocity is contained in a
%|                  "pure" scalar first quaternion in the sense of
%|                  Equation (182) on Page 466 of Reference [ 2 ]
%|                  (except that Equation (182) is scalar last).
%|                  Also note Equation (30) on Page 265 of
%|                  Reference [ 1 ]
%|
%|-----------------------------------------------------------------------------
%|
%|  REFERENCE(s):
%|
%|    [ 1 ]  "The Kinematic Equation for the Rotation Vector",
%|           Malcolm D. Shuster,
%|           IEEE Transactions on Aerospace and Electronic Systems,
%|           Vol. 29, No. 1, pp. 263 - 267,
%|           January 1993
%|           malcolmdshuster.com/Pub_1993c_J_RotVec_IEEE.pdf
%|
%|    [ 2 ]  "A Survey of Attitude Representations",
%|           Malcolm D. Shuster,
%|           The Journal of the Astronautical Sciences,
%|           Vol. 41, No. 4, pp. 430 - 517,
%|           October-December, 1993
%|           malcolmdshuster.com/Pub_1993h_J_Repsurv_scan.pdf
%|
%|------------------------------------------------------------------------------
%|
%|  USAGE:
%|
%|    [                                                         ...
%|      actualSimulationTimes,                                  ...
%|      scalarFirstUnitQuaternions,                             ...
%|      angularVelocitiesBodyCoords                             ...
%|    ] = generateSolutionsForFreeRotationRigidBodyKinematics   ,,,
%|                (                                             ...
%|                  J,                                          ...
%|                  initialScalarFirstUnitQuaternion,           ...
%|                  intitalAngularVelocityBodyCoords            ...
%|                  numberTimeSamples,                          ...
%|                  simulationEndTime                           ...
%|                );
%|
%===============================================================================

%===============================================================================
%        1         2         3         4         5         6         7         8
%2345678901234567890123456789012345678901234567890123456789012345678901234567890
%===============================================================================

%{------------------------------------------------------------------------------
   expectedNumberInputArguments   =  5;
   expectedNumberOutputArguments  =  3;
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     actualNumberInputArguments   =  nargin;
     actualNumberOutputArguments  =  nargout;
%-------------------------------------------------------------------------------
   if( actualNumberInputArguments == expectedNumberInputArguments )
    %{--------------------------------------------------------------------------
    %  Have encountered expected number of function input arguments.
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %  Extract input arguments from variable input argument cell array.
    %---------------------------------------------------------------------------
       J                                = varargin{ 1 };
       initialScalarFirstUnitQuaternion = varargin{ 2 };
       initialAngularVelocityBodyCoords = varargin{ 3 };
       numberTimeSamples                = varargin{ 4 };
       simulationEndTime                = varargin{ 5 };
    %---------------------------------------------------------------------------
    %  Check number of function output arguments.
    %---------------------------------------------------------------------------
       if( actualNumberOutputArguments == expectedNumberOutputArguments )
        %{----------------------------------------------------------------------
        %  Have encountered expected number of function output arguments.
        %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        %  Continue processing.
        %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        %  Define a six component state vector consisting of:
        %     A three element rotation vector.
        %     A three element angular velocity vector.
        %-----------------------------------------------------------------------
           initialStateVector = [                                            ...
                                  initialScalarFirstUnitQuaternion( : );     ...
                                  initialAngularVelocityBodyCoords( : )      ...
                                ];
        %-----------------------------------------------------------------------
        %  Set the error tolerance option values for the ODE solver.
        %-----------------------------------------------------------------------
           ordinaryDifferentialEquationsSolverOptions =                      ...
                   odeset                                                    ...
                     (                                                       ...
                       'RelTol',   1.0E-13,                                  ...
                       'AbsTol',   1.0E-13                                   ...
                     );
        %-----------------------------------------------------------------------
        %  Obtain the principal Moment of Inertia values from the specified
        %  Moment of Inertia matrix.
        %-----------------------------------------------------------------------
           principalMomentInertiaVector = 1.0 ./ linsolve( J, ones( 3, 1 ) );
        %-----------------------------------------------------------------------
        %  Generate a vector of specified simulation time sample values.
        %-----------------------------------------------------------------------
           specifiedSimulationTimes     = linspace                           ...
                                             (                               ...
                                               0.0,                          ...
                                               simulationEndTime,            ...
                                               numberTimeSamples             ...
                                             );
        %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           stateDerivativesFunctionHandle =                                  ...
                  @(                                                         ...
                     currentTime,                                            ...
                     currentStateVector                                      ...
                   )                                                         ...
                  determineStateDerivatives                                  ...
                           (                                                 ...
                             currentTime,                                    ...
                             currentStateVector,                             ...
                             principalMomentInertiaVector                    ...
                           );
        %-----------------------------------------------------------------------
        %  Use the Matlab implementation of the Dorman-Prince 45 explicit
        %  variable time step Runge-Kutta method to solve the equations of
        %  motion.
        %-----------------------------------------------------------------------
           [                                                                 ...
             actualSimulationTimes,                                          ...
             simulationStateVectors                                          ...
           ] = ode45                                                         ...
                (                                                            ...
                  stateDerivativesFunctionHandle,                            ...
                  specifiedSimulationTimes,                                  ...
                  initialStateVector,                                        ...
                  ordinaryDifferentialEquationsSolverOptions                 ...
                );
        %-----------------------------------------------------------------------
        % 
        % NOTE(s):
        % 
        %   [ 1 ] State Vector Components:
        % 
        %         [ 1.1 ]  State Vector components 1, 2, 3 and 4 are the
        %                  scalar first unit quaternions.
        % 
        %         [ 1.2 ]  State Vector components 4, 5 and 6 are the
        %                  angular velocity vector.
        % 
        %-----------------------------------------------------------------------
           scalarFirstUnitQuaternions  = simulationStateVectors              ...
                                                       ( :, 1 : 4 )';
           angularVelocitiesBodyCoords = simulationStateVectors              ...
                                                       ( :, 5 : 7 )';
        %}----------------------------------------------------------------------
       else
        %{----------------------------------------------------------------------
        %  Have not encountered expected number of function output arguments.
        %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        %  Generate a purpose message.
        %-----------------------------------------------------------------------
           generateKinematicsSolutionsPurposeMessage(  );
        %-----------------------------------------------------------------------
        %  Generate a usage message.
        %-----------------------------------------------------------------------
           generateKinematicsSolutionsUsageMessage(  );
        %-----------------------------------------------------------------------
        %  Generate an error message.
        %-----------------------------------------------------------------------
           STDOUT = 1;
        %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           fprintf                                                           ...
           (                                                                 ...
             STDOUT,                                                         ...
             [                                                               ...
               '\n\n\n'                                                      ...
               '%s\n%s\n%s\n%s\n%s\n%s\n%s\n'                                ...
               '%s%d\n'                                                      ...
               '%s%d\n'                                                      ...
               '%s\n%s\n%s\n%s\n'                                            ...
               '\n\n\n'                                                      ...
             ],                                                              ...
             '============================================================', ...
             '|',                                                            ...
             '|  ERROR:',                                                    ...
             '|',                                                            ...
             '|    Have encountered unexpected number function output',      ...
             '|    arguments.',                                              ...
             '|',                                                            ...
             '|    Expected number function output arguments:-->',           ...
             expectedNumberOutputArguments,                                  ...
             '|      Actual number function output arguments:-->',           ...
               actualNumberOutputArguments,                                  ...
             '|',                                                            ...
             '|    This is an error.',                                       ...
             '|',                                                            ...
             '============================================================'  ...
           );
        %-----------------------------------------------------------------------
           scalarFirstUnitQuaternions  = [ ];
           angularVelocitiesBodyCoords = [ ];
        %-----------------------------------------------------------------------
        %  Terminate the program.
        %-----------------------------------------------------------------------
           error                                                             ...
             (                                                               ...
               [                                                             ...
                 'Have encountered unexpected number of function '           ...
                 'output arguments.'                                         ...
               ]                                                             ...
             );
        %}----------------------------------------------------------------------
       end;
    %}--------------------------------------------------------------------------
   else
    %{--------------------------------------------------------------------------
    %  Have not encountered expected number of function input arguments.
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %  Generate a purpose message.
    %---------------------------------------------------------------------------
       generateKinematicsSolutionsPurposeMessage(  );
    %---------------------------------------------------------------------------
    %  Generate a usage message.
    %---------------------------------------------------------------------------
       generateKinematicsSolutionsUsageMessage(  );
    %---------------------------------------------------------------------------
    %  Generate an error message.
    %---------------------------------------------------------------------------
       STDOUT = 1;
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       fprintf                                                               ...
       (                                                                     ...
         STDOUT,                                                             ...
         [                                                                   ...
           '\n\n\n'                                                          ...
           '%s\n%s\n%s\n%s\n%s\n%s\n%s\n'                                    ...
           '%s%d\n'                                                          ...
           '%s%d\n'                                                          ...
           '%s\n%s\n%s\n%s\n'                                                ...
           '\n\n\n'                                                          ...
         ],                                                                  ...
         '============================================================',     ...
         '|',                                                                ...
         '|  ERROR:',                                                        ...
         '|',                                                                ...
         '|    Have encountered unexpected number function input',           ...
         '|    arguments.',                                                  ...
         '|',                                                                ...
         '|    Expected number function input arguments:-->',                ...
         expectedNumberInputArguments,                                       ...
         '|      Actual number function input arguments:-->',                ...
           actualNumberInputArguments,                                       ...
         '|',                                                                ...
         '|    This is an error.',                                           ...
         '|',                                                                ...
         '============================================================'      ...
       );
    %---------------------------------------------------------------------------
       scalarFirstUnitQuaternions  = [ ];
       angularVelocitiesBodyCoords = [ ];
    %---------------------------------------------------------------------------
    %  Terminate the program.
    %---------------------------------------------------------------------------
       error                                                                 ...
        ( 'Have encountered unexpected number of function input arguments.' );
    %}--------------------------------------------------------------------------
   end;
%------------------------------------------------------------------------------- 
   return;
%}------------------------------------------------------------------------------


%===============================================================================
function                                                                     ...
[                                                                            ...
  currentStateDerivatives                                                    ...
] =                                                                          ...
determineStateDerivatives                                                    ...
         (                                                                   ...
           currentTime,                                                      ...
           currentStateVector,                                               ...
           principalInertiaMomentVector                                      ...
         )
%{------------------------------------------------------------------------------
%
%  NOTE(s):
%
%    [ 1 ]  Rotating object kinematics is "scalar first unit quaternion"
%           based.
%
%    [ 2 ]  Scalar First Unit Quaternion Definition:
%|
%           Use Equation (22) on Page 264 of Reference [ 1 ]
%           EXCEPT:
%             [ 2.1 ]  Equation (22) is scalar LAST.
%             [ 2.2 ]  Here scalar FIRST is used.
%
%    [ 3 ] State Vector Components:
%
%          [ 3.1 ]  State Vector components 1, 2, 3 and 4 are the
%                   "scalar first unit quaternion".
%
%          [ 3.2 ]  State Vector components 5, 6 and 7 are the
%                   angular velocity vector.
%
%    [ 4 ] Time rate of change of the "scalar first unit quaternion":
%
%          Use the combination of the following two items:
%
%          [ 4.1 ] The FIRST PART of Equation (46) on Page 266
%                  of Reference [ 1 ].
%                  This states that the time rate of change of the
%                  attitude quaterion is equal to the cross product
%                  of two quaternions"
%                  [ 5.1.1 ] The pure quaternion container of the
%                            angular velocity 3D vector.
%                  [ 5.1.2 ] The attitude quaternion.
%
%          [ 4.2 ] Equation (305) on Page 482 of Reference [ 2 ].
%                  This equation shows how to use the quaternion
%                  product when the quaternion cross product is
%                  specified.
%                  Here the angular velocity is contained in a
%                  "pure" scalar first quaternion in the sense of
%                  Equation (182) on Page 466 of Reference [ 2 ]
%                  (except that Equation (182) is scalar last).
%
%------------------------------------------------------------------------------
%
%
%  REFERENCE(s):
%
%    [ 1 ]  "The Kinematic Equation for the Rotation Vector",
%           Malcolm D. Shuster,
%           IEEE Transactions on Aerospace and Electronic Systems,
%           Vol. 29, No. 1, pp. 263 - 267,
%           January 1993
%           malcolmdshuster.com/Pub_1993c_J_RotVec_IEEE.pdf
%
%    [ 2 ]  "A Survey of Attitude Representations",
%           Malcolm D. Shuster,
%           The Journal of the Astronautical Sciences,
%           Vol. 41, No. 4, pp. 430 - 517,
%           October-December, 1993
%           malcolmdshuster.com/Pub_1993h_J_Repsurv_scan.pdf
%
%-------------------------------------------------------------------------------
   X_INDEX                            = 1;
   Y_INDEX                            = 2;
   Z_INDEX                            = 3;
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   STATE_QUATERNION_S_INDEX           = 1;
   STATE_QUATERNION_X_INDEX           = 2;
   STATE_QUATERNION_Y_INDEX           = 3;
   STATE_QUATERNION_Z_INDEX           = 4;
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   STATE_ANG_VEL_X_INDEX              = 5;
   STATE_ANG_VEL_Y_INDEX              = 6;
   STATE_ANG_VEL_Z_INDEX              = 7;
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   STATE_DERIV_QUATERNAION_S_INDEX    = 1;
   STATE_DERIV_QUATERNAION_X_INDEX    = 2;
   STATE_DERIV_QUATERNAION_Y_INDEX    = 3;
   STATE_DERIV_QUATERNAION_Z_INDEX    = 4;
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   STATE_DERIV_ANG_ACCEL_X_INDEX      = 5;
   STATE_DERIV_ANG_ACCEL_Y_INDEX      = 6;
   STATE_DERIV_ANG_ACCEL_Z_INDEX      = 7;
%-------------------------------------------------------------------------------
   currentStateDerivatives            = zeros( size( currentStateVector ) );
%-------------------------------------------------------------------------------
   currentScalarFirstUnitQuaternion   = currentStateVector                   ...
                                               (                             ...
                                                 [                           ...
                                                   STATE_QUATERNION_S_INDEX; ...
                                                   STATE_QUATERNION_X_INDEX; ...
                                                   STATE_QUATERNION_Y_INDEX; ...
                                                   STATE_QUATERNION_Z_INDEX  ...
                                                 ]                           ...
                                               );
%-------------------------------------------------------------------------------
   currentAngularVelocityVectorBodyCoords                                    ...
                                      = currentStateVector                   ...
                                               (                             ...
                                                 [                           ...
                                                   STATE_ANG_VEL_X_INDEX;    ...
                                                   STATE_ANG_VEL_Y_INDEX;    ...
                                                   STATE_ANG_VEL_Z_INDEX     ...
                                                 ]                           ...
                                               );
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   pureQuaternionContainerForAngVelVectorBodyCoords =                        ...
                             [                                               ...
                               0.0;                                          ...
                               currentAngularVelocityVectorBodyCoords( : )   ...
                             ];
%-------------------------------------------------------------------------------
%  Compute product of two scalar first (not necessarily unit) quaternions.
%-------------------------------------------------------------------------------
   [                                                                         ...
     angVelQuaternionTimesCurrentQuaternion                                  ...
   ] = scalarFirstQuaternionMultiplication                                   ...
             (                                                               ...
               currentScalarFirstUnitQuaternion(                 : ),        ...
               pureQuaternionContainerForAngVelVectorBodyCoords( : )         ...
             );
%-------------------------------------------------------------------------------
%
%  Compute quaternion time rate of change.
%
%-------------------------------------------------------------------------------
   quaternionDerivative       = 0.5 *                                        ...
                                angVelQuaternionTimesCurrentQuaternion;
%-------------------------------------------------------------------------------
   currentStateDerivatives                                                   ...
         ( STATE_DERIV_QUATERNAION_S_INDEX ) = quaternionDerivative          ...
                                                ( STATE_QUATERNION_S_INDEX );
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   currentStateDerivatives                                                   ...
         ( STATE_DERIV_QUATERNAION_X_INDEX ) = quaternionDerivative          ...
                                                ( STATE_QUATERNION_X_INDEX );
   currentStateDerivatives                                                   ...
         ( STATE_DERIV_QUATERNAION_Y_INDEX ) = quaternionDerivative          ...
                                                ( STATE_QUATERNION_Y_INDEX );
   currentStateDerivatives                                                   ...
         ( STATE_DERIV_QUATERNAION_Z_INDEX ) = quaternionDerivative          ...
                                                ( STATE_QUATERNION_Z_INDEX );
%-------------------------------------------------------------------------------
%
%  Use Euler's equations to determine the angular acceleration in terms of
%  moments of inertia (J) and angular velocities (w).
%
%-------------------------------------------------------------------------------
   Jx                               = principalInertiaMomentVector( X_INDEX );
   Jy                               = principalInertiaMomentVector( Y_INDEX );
   Jz                               = principalInertiaMomentVector( Z_INDEX );
%-------------------------------------------------------------------------------
   wx                               = currentAngularVelocityVectorBodyCoords ...
                                                                  ( X_INDEX );
   wy                               = currentAngularVelocityVectorBodyCoords ...
                                                                  ( Y_INDEX );
   wz                               = currentAngularVelocityVectorBodyCoords ...
                                                                  ( Z_INDEX );
%-------------------------------------------------------------------------------
   xAngularAccelerationBodyCoords   = (  ( Jy - Jz ) / Jx ) * wy * wz;
   yAngularAccelerationBodyCoords   = ( -( Jx - Jz ) / Jy ) * wx * wz;
   zAngularAccelerationBodyCoords   = (  ( Jx - Jy ) / Jz ) * wx * wy;
%-------------------------------------------------------------------------------
   currentStateDerivatives                                                   ...
        ( STATE_DERIV_ANG_ACCEL_X_INDEX ) = xAngularAccelerationBodyCoords;
   currentStateDerivatives                                                   ...
        ( STATE_DERIV_ANG_ACCEL_Y_INDEX ) = yAngularAccelerationBodyCoords;
   currentStateDerivatives                                                   ...
        ( STATE_DERIV_ANG_ACCEL_Z_INDEX ) = zAngularAccelerationBodyCoords;
%-------------------------------------------------------------------------------
   return;
%}------------------------------------------------------------------------------
