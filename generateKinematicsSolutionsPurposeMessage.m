function                                                                     ...
generateKinematicsSolutionsPurposeMessage(  )
%{------------------------------------------------------------------------------
   STDOUT = 1;
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   fprintf                                                                   ...
   (                                                                         ...
     STDOUT,                                                                 ...
     [                                                                       ...
       '\n\n\n'                                                              ...
       '%s%s\n'                                                              ...
       '%s\n%s\n%s\n%s\n%s\n'                                                ...
       '%s%s\n'                                                              ...
       '%s\n%s\n%s\n%s\n%s\n%s\n'                                            ...
       '%s%s\n'                                                              ...
       '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n'                            ...
       '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n'                            ...
       '%s\n%s\n%s\n%s\n%s\n%s\n%s\n'                                        ...
       '%s%s\n'                                                              ...
       '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n'                            ...
       '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n'                            ...
       '%s\n'                                                                ...
       '%s%s\n'                                                              ...
       '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n'                            ...
       '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n'                            ...
       '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n'                            ...
       '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n'                            ...
       '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n'                            ...
       '%s%s\n'                                                              ...
       '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n'                            ...
       '%s\n%s\n%s\n%s\n%s\n%s\n%s\n'                                        ...
       '%s%s\n'                                                              ...
       '\n\n\n'                                                              ...
     ],                                                                      ...
     '====================================================================', ...
     '===========',                                                          ...
     '|',                                                                    ...
     '|  FUNCTION:',                                                         ...
     '|',                                                                    ...
     '|    generateSolutionsForFreeRotationRigidBodyKinematics',             ...
     '|',                                                                    ...
     '|-------------------------------------------------------------------', ...
     '----------',                                                           ...
     '|',                                                                    ...
     '|  PURPOSE:',                                                          ...
     '|',                                                                    ...
     '|    Generate solutions at specified times for torque free rotation',  ...
     '|    of rigid bodies.',                                                ...
     '|',                                                                    ...
     '|-------------------------------------------------------------------', ...
     '----------',                                                           ...
     '|',                                                                    ...
     '|  INPUTS:',                                                           ...
     '|',                                                                    ...
     '|    J',                                                               ...
     '|      A three by three matrix of double precision floating point',    ...
     '|      Inertia Tensor values.',                                        ...
     '|      UNIT(s):  [kilograms / meter squared]',                         ...
     '|',                                                                    ...
     '|    initialScalarFirstUnitQuaternion',                                ...
     '|      A one dimensional vector of four elements containing the',      ...
     '|      initial scalar first unit quaternion.',                         ...
     '|      UNIT(s):  [nondimensional]',                                    ...
     '|',                                                                    ...
     '|    intitalAngularVelocityBodyCoords',                                ...
     '|      A one dimensional vector of length three containing the',       ...
     '|      initial angular velocity vector in non-inertial body',          ...
     '|      coordinates.',                                                  ...
     '|      UNIT(s):  [radians/second]',                                    ...
     '|',                                                                    ...
     '|    numberTimeSamples',                                               ...
     '|      Number of time sample values generated for the results.',       ...
     '|      UNIT(s):  [nondimensional]',                                    ...
     '|',                                                                    ...
     '|    simulationEndTime',                                               ...
     '|      Maximum time value for the simulation.',                        ...
     '|      UNIT(s):  [seconds[',                                           ...
     '|',                                                                    ...
     '|-------------------------------------------------------------------', ...
     '----------',                                                           ...
     '|',                                                                    ...
     '|  OUTPUTS:',                                                          ...
     '|',                                                                    ...
     '|    actualSimulationTimes',                                           ...
     '|      A one dimensional vector the ith element of which is the',      ...
     '|      actual ith simulation time value.',                             ...
     '|      UNIT(s):  [seconds]',                                           ...
     '|',                                                                    ...
     '|    scalarFirstUnitQuaternions',                                      ...
     '|      A two dimensional matrix having four rows and as many columns', ...
     '|      as there are time samples.',                                    ...
     '|      The jth column of this matrix is the scalar first unit',        ...
     '|      quaternion for the jth sample time.',                           ...
     '|      UNIT(s):  [nondimensional]',                                    ...
     '|',                                                                    ...
     '|    angularVelocitiesBodyCoords',                                     ...
     '|      A two dimensional vector the ith element of which is the',      ...
     '|      angular velocity in non-inertial body coordinates at the',      ...
     '|      ith sample time.',                                              ...
     '|      UNIT(s):  [radians/second]',                                    ...
     '|',                                                                    ...
     '|-------------------------------------------------------------------', ...
     '----------',                                                           ...
     '|',                                                                    ...
     '|  NOTE(s):',                                                          ...
     '|',                                                                    ...
     '|    [ 1 ]  Rotating object kinematics is "scalar first unit',         ...
     '|           quaternion" based.',                                       ...
     '|',                                                                    ...
     '|    [ 2 ]  Scalar First Unit Quaternion Definition:',                 ...
     '|           Use Equation (22) on Page 264 of Reference [ 1 ]',         ...
     '|           EXCEPT:',                                                  ...
     '|             [ 1 ]  Equation (22) is scalar LAST.',                   ...
     '|             [ 2 ]  Here scalar FIRST is used.',                      ...
     '|',                                                                    ...
     '|    [ 3 ]  Ordinary Differential Equations (ODEs) are solved using',  ...
     '|           the Matlab "ode45" function to apply the Dorman-Prince',   ...
     '|           4(5) explicit embedded variable time step Runge-Kutta',    ...
     '|           method.',                                                  ...
     '|',                                                                    ...
     '|    [ 4 ] State Vector Components:',                                  ...
     '|',                                                                    ...
     '|          [ 4.1 ]  State Vector components 1, 2, 3 and 4 are the',    ...
     '|                   "scalar first unit quaternion".',                  ...
     '|',                                                                    ...
     '|          [ 4.2 ]  State Vector components 5, 6 and 7 are the'        ...
     '|                   angular velocity vector.',                         ...
     '|',                                                                    ...
     '|    [ 5 ] Time rate of change of the "scalar first unit',             ...
     '|          quaternion":',                                              ...
     '|',                                                                    ...
     '|          Use the combination of the following two items:',           ...
     '|',                                                                    ...
     '|          [ 5.1 ] The FIRST PART of Equation (46) on Page 266',       ...
     '|                  of Reference [ 1 ].',                               ...
     '|                  This states that the time rate of change of the',   ...
     '|                  attitude quaterion is equal to the cross product',  ...
     '|                  of two quaternions"',                               ...
     '|                  [ 5.1.1 ] The pure quaternion container of the',    ...
     '|                            angular velocity 3D vector.',             ...
     '|                  [ 5.1.2 ] The attitude quaternion.',                ...
     '|',                                                                    ...
     '|          [ 5.2 ] Equation (305) on Page 482 of Reference [ 2 ].'     ...
     '|                  This equation shows how to use the quaternion',     ...
     '|                  product when the quaternion cross product is',      ...
     '|                  specified.',                                        ...
     '|                  Here the angular velocity is contained in a',       ...
     '|                  "pure" scalar first quaternion in the sense of',    ...
     '|                  Equation (182) on Page 466 of Reference [ 2 ]',     ...
     '|                  (except that Equation (182) is scalar last).',      ...
     '|                  Also note Equation (30) on Page 265 of',            ...
     '|                  Reference [ 1 ]',                                   ...
     '|',                                                                    ...
     '|-------------------------------------------------------------------', ...
     '----------',                                                           ...
     '|',                                                                    ...
     '|  REFERENCE(s):',                                                     ...
     '|',                                                                    ...
     '|    [ 1 ]  "The Kinematic Equation for the Rotation Vector",',        ...
     '|           Malcolm D. Shuster,',                                      ...
     '|           IEEE Transactions on Aerospace and Electronic Systems,',   ...
     '|           Vol. 29, No. 1, pp. 263 - 267,',                           ...
     '|           January 1993',                                             ...
     '|           malcolmdshuster.com/Pub_1993c_J_RotVec_IEEE.pdf',          ...
     '|',                                                                    ...
     '|    [ 2 ]  "A Survey of Attitude Representations",',                  ...
     '|           Malcolm D. Shuster,',                                      ...
     '|           The Journal of the Astronautical Sciences,',               ...
     '|           Vol. 41, No. 4, pp. 430 - 517,',                           ...
     '|           October-December, 1993',                                   ...
     '|           malcolmdshuster.com/Pub_1993h_J_Repsurv_scan.pdf',         ...
     '|',                                                                    ...
     '====================================================================', ...
     '==========='                                                           ...
   );
%------------------------------------------------------------------------------- 
   return;
%}------------------------------------------------------------------------------
