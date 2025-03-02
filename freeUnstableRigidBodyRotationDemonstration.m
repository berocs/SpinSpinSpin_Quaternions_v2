function                                                                     ...
freeUnstableRigidBodyRotationDemonstration(  )
%===============================================================================
%|
%|  FUNCTION:
%|
%|    freeUnstableRigidBodyRotationDemonstration
%|
%|-----------------------------------------------------------------------------
%|
%|  PURPOSE:
%|
%|    Perform a deomonstration of free unstable rotation of a rigid body.
%|
%|-----------------------------------------------------------------------------
%|
%|  INPUT(s):
%|
%|    None
%|
%|-----------------------------------------------------------------------------
%|
%|  OUTPUT(s):
%|
%|    None
%|
%|-----------------------------------------------------------------------------
%|
%|  RETURN VALUE(s):
%|
%|    None
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
%|    [ 5 ] Time rate of change of the "scalar first unit quaternion":
%|
%|          Use the combination of the following two items:
%|
%|          [ 5.1 ] The FIRST PART of Equation (46) on Page 266
%|                  of Reference [ 1 ].
%|                  This states that the time rate of change of the
%|                  attitude quaterion is equal to the cross product
%|                  of two quaternions"
%|                  [ 5.1.1 ] The pure quaternion container of the
%|                            angular velocity 3D vector.
%|                  [ 5.1.2 ] The attitude quaternion.
%|
%|          [ 5.2 ] Equation (305) on Page 482 of Reference [ 2 ].
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
%|    [ 3 ]  "The Bizarre Behavior of Rotating Bodies -
%|            The Dzhanibekov Effect",
%|            Derek Alexander Muller
%|            YouTube channel Veritasium
%|            Spinning objects have strange instabilities known as
%|            "The Dzhanibekov Effect" or "Tennis Racket Theorem" -
%|            this video offers an intuitive explanation.
%|            https://www.youtube.com/watch?v=1VPfZ_XzisU
%|
%|    [ 4 ]  Elements of the following were used to help create this
%|           project:
%|
%|           [ 4.1 ]  https://github.com/alaricgregoire/EulerFreeBody
%|                    www.mathworks.com/matlabcentral/fileexchange/
%|                      75265-euler-free-body-motion
%|
%|           [ 4.2 ]  github.com/phymhan/matlab-axis-label-alignment
%|                    www.mathworks.com/matlabcentral/fileexchange/
%|                      49542-phymhan-matlab-axis-label-alignment
%|
%|           [ 4.3 ]  www.mathworks.com/matlabcentral/fileexchange/
%|                      25372-marrow3-m-easy-to-use-3d-arrow
%|
%|-----------------------------------------------------------------------------
%|
%|  USAGE:
%|
%|    freeUnstableRigidBodyRotationDemonstration(  )
%|
%===============================================================================

%===============================================================================
%        1         2         3         4         5         6         7         8
%2345678901234567890123456789012345678901234567890123456789012345678901234567890
%===============================================================================

%{------------------------------------------------------------------------------
   close   all;
   clear   all;
%-------------------------------------------------------------------------------
   STDOUT = 1;
%-------------------------------------------------------------------------------
   [                                                                         ...
     operatingSystemArchitectureName                                         ...
   ] = computer(  );
%-------------------------------------------------------------------------------
   WINDOW_OS_INDICATOR = 'PCWIN64';
   LINUX_OS_INDICATOR  = 'GLNXA64';
%-------------------------------------------------------------------------------
   switch( operatingSystemArchitectureName )
        %{----------------------------------------------------------------------
           case( 'PCWIN64' )
              %{----------------------------------------------------------------
              %  Operating system is MS Windows 64 bit.
              %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              %  Set directory separator character.
              %-----------------------------------------------------------------
                 directorySeparator = '\';
              %}----------------------------------------------------------------
        %-----------------------------------------------------------------------
           case( 'GLNXA64' )
              %{----------------------------------------------------------------
              %  Operating system is Linux 64 bit.
              %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              %  Set directory separator character.
              %-----------------------------------------------------------------
                 directorySeparator = '/';
              %}----------------------------------------------------------------
        %-----------------------------------------------------------------------
           otherwise
              %{----------------------------------------------------------------
              %  Default operating system.
              %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              %  There is no default operating system.
              %  Generate an error message.
              %-----------------------------------------------------------------
                 [                                                           ...
                   errorMessageString                                        ...
                 ] =                                                         ...
                   sprintf                                                   ...
                     (                                                       ...
                       [                                                     ...
                         '\n\n\n'                                            ...
                         '%s\n%s\n%s\n%s\n%s\n'                              ...
                         '%s"%s"\n'                                          ...
                         '%s\n%s\n%s\n%s\n'                                  ...
                         '\n\n\n'                                            ...
                       ],                                                    ...
                       '--------------------------------------------------', ...
                       '|',                                                  ...
                       '| ERROR:',                                           ...
                       '|',                                                  ...
                       '|   Encountered unknown operating system.',          ...
                       '|   Operating System Indicator was: ',               ...
                       operatingSystemArchitectureName,                      ...
                       '|',                                                  ...
                       '|   Program will terminate.',                        ...
                       '|',                                                  ...
                       '--------------------------------------------------'  ...
                     );
              %-----------------------------------------------------------------
              %  Terminate the program.
              %-----------------------------------------------------------------
                 error( errorMessageString );
              %}----------------------------------------------------------------
        %}----------------------------------------------------------------------
   end;
%===============================================================================
   currentWorkingDirectoryPath = cd;
%-------------------------------------------------------------------------------
   sourceCodeDirectory1 = [                                                  ...
                            currentWorkingDirectoryPath,                     ...
                            directorySeparator,                              ...
                            'Utilities'                
                          ];
%-------------------------------------------------------------------------------
   addpath( genpath( sourceCodeDirectory1 ) );
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%  Generate purpose message.
%-------------------------------------------------------------------------------
   generateFreeUnstableRigidBodyRotationDemoPurposeMessage(   );
%-------------------------------------------------------------------------------
%  Generate usage message.
%-------------------------------------------------------------------------------
   generateFreeUnstableRigidBodyRotationDemoUsageMessage(   );
%-------------------------------------------------------------------------------
%  Define the unstable free rotation of a rigid body demonstration.
%-------------------------------------------------------------------------------
   fprintf                                                                   ...
    (                                                                        ...
      STDOUT,                                                                ...
      [                                                                      ...
        '\n\n\n'                                                             ...
        '%s\n%s\n%s\n%s\n%s\n'                                               ...
        '\n\n\n'                                                             ...
      ],                                                                     ...
      '------------------------------------------------------------------',  ...
      '|',                                                                   ...
      '|  UNSTABLE FREE ROTATION RIGID BODY DEMONSTRATION:',                 ...
      '|',                                                                   ...
      '------------------------------------------------------------------'   ...
    );
%-------------------------------------------------------------------------------
%
%  Specify unstable free rotation of a rigid 3D rectangular block by
%  assigning dominant angular velocity about body axis with intermediate
%  inertia value (body axis length proportional to inertia value).
%
%-------------------------------------------------------------------------------
   simulationEndTime                = 1.0E3;
   numberSampleTimes                = 1.0E6;
          %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   momentOfInertiaMatrix            = diag( [         4.0,   1.0,      9.0 ] );
   initialScalarFirstUnitQuaternion =       [                                ...
                                               cos( 0.5 );                   ...
                                               sin( 0.5 );                   ...
                                                    0.0;                     ...
                                                    0.0                      ...
                                            ];
   initialAngularVelocityBodyCoord  =       [         0.10;  0.00001;  0.0 ];
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   [                                                                         ...
     demonstrationTitle                                                      ...
   ] = sprintf                                                               ...
        (                                                                    ...
          '%s',                                                              ...
          'UNSTABLE FREE ROTATION RIGID BODY DEMONSTRATION'                  ...
        );
%-------------------------------------------------------------------------------
   colorDefinitions;
%-------------------------------------------------------------------------------
%  profile on;
%-------------------------------------------------------------------------------
   fprintf                                                                   ...
    (                                                                        ...
      STDOUT,                                                                ...
      [                                                                      ...
        '\n\n\n'                                                             ...
        '%s\n%s\n%s\n%s\n%s\n'                                               ...
        '\n\n\n'                                                             ...
      ],                                                                     ...
      '------------------------------------------------------------------',  ...
      '|',                                                                   ...
      '|  Generating body kinematics at each time sample:',                  ...
      '|',                                                                   ...
      '------------------------------------------------------------------'   ...
    );
%-------------------------------------------------------------------------------
%  Generate solution for specified torque free rotation of specified rigid
%  body at each sample time.
%-------------------------------------------------------------------------------
   [                                                                         ...
     timeSampleValues,                                                       ...
     scalarFirstUnitQuaternions,                                             ...
     angularVelocitiesBodyCoords                                             ...
   ] = generateSolutionsForFreeRotationRigidBodyKinematics                   ...
               (                                                             ...
                 momentOfInertiaMatrix,                                      ...
                 initialScalarFirstUnitQuaternion,                           ...
                 initialAngularVelocityBodyCoord,                            ...
                 numberSampleTimes,                                          ...
                 simulationEndTime                                           ...
               );
%-------------------------------------------------------------------------------
%  Angular momentum values is calculated as the product of the moment
%  of inertia matrix and the angular velocity vectors.
%-------------------------------------------------------------------------------
   angularMomentumsBodyCoords = momentOfInertiaMatrix *                      ...
                                angularVelocitiesBodyCoords;
%-------------------------------------------------------------------------------
   numberTimeSamples = length( timeSampleValues );
%-------------------------------------------------------------------------------
   [                                                                         ...
     angularVelocitiesInertialCoords,                                        ...
     angularMomentumsInertialCoords,                                         ...
     rotationalKineticEnergies                                               ...
   ] = generateInertialResults                                               ...
               (                                                             ...
                 scalarFirstUnitQuaternions,                                 ...
                 angularVelocitiesBodyCoords,                                ...
                 angularMomentumsBodyCoords                                  ...
               );
%-------------------------------------------------------------------------------
%  profsave;
%-------------------------------------------------------------------------------
   fprintf                                                                   ...
    (                                                                        ...
      STDOUT,                                                                ...
      [                                                                      ...
        '\n\n\n'                                                             ...
        '%s\n%s\n%s\n%s\n%s\n'                                               ...
        '\n\n\n'                                                             ...
      ],                                                                     ...
      '------------------------------------------------------------------',  ...
      '|',                                                                   ...
      '|  Finished generating body kinematics at each time sample:',         ...
      '|',                                                                   ...
      '------------------------------------------------------------------'   ...
    );
%-------------------------------------------------------------------------------
   input( 'Press Enter to Continue with 2D plots:-->' );
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   generatePlots2D                                                           ...
           (                                                                 ...
              colors,                                                        ...
              timeSampleValues,                                              ...
              scalarFirstUnitQuaternions,                                    ...
              angularVelocitiesBodyCoords,                                   ...
              angularMomentumsBodyCoords,                                    ...
              angularVelocitiesInertialCoords,                               ...
              angularMomentumsInertialCoords,                                ...
              rotationalKineticEnergies,                                     ...
              demonstrationTitle                                             ...
           );
%-------------------------------------------------------------------------------
   input( 'Press Enter to Continue with 3D plots:-->' );
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   generatePlots3D                                                           ...
           (                                                                 ...
              colors,                                                        ...
              angularVelocitiesBodyCoords,                                   ...
              angularMomentumsBodyCoords,                                    ...
              rotationalKineticEnergies,                                     ...
              momentOfInertiaMatrix,                                         ...
              demonstrationTitle                                             ...
           );
%-------------------------------------------------------------------------------
   input( 'Press Enter to Continue with 3D Animation:-->' );
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   skipRateFactor = 200;
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   referencePosition = [                                                     ...
                          0.0;                                               ...
                          0.0;                                               ...
                          0.0                                                ...
                       ];
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   performNonRealTimeAnimationPlayBack3D                                     ...
          (                                                                  ...
            colors,                                                          ...
            timeSampleValues,                                                ...
            referencePosition,                                               ...
            scalarFirstUnitQuaternions,                                      ...
            angularVelocitiesInertialCoords,                                 ...
            angularMomentumsInertialCoords,                                  ...
            momentOfInertiaMatrix,                                           ...
            skipRateFactor,                                                  ...
            demonstrationTitle                                               ...
          );
%-------------------------------------------------------------------------------
   return;
%}------------------------------------------------------------------------------

