# SpinSpinSpin_Quaternions_v2

A Simple MATLAB Simulation for Unstable Free Rotation of a Rigid Body.

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=berocs/SpinSpinSpin_Quaternions_v2)

NOTE(s):
- Rotating object kinematics is "**scalar first unit quaternion**" based.
- The project "**SpinSpinSpin_Quaternion_v1**" is similar to this project.
- The project "**SpinSpinSpin_RotationVector**" is similar to this project but "**RotationVector**" based.
- Scalar First Unit Quaternion Definition:
  - A scalar first unit quaternion is a one dimensional row or column vector of length four.
    - **Use Equation (22) on Page 264 of Reference [ 1 ]** 
      - **EXCEPT:**
       - **[ 1 ]  Equation (22) is scalar LAST.**
       - **[ 2 ]  Here scalar FIRST is used.**
- Ordinary Differential Equations (ODEs) are solved using the  Matlab "ode45" function to apply the Dorman-Prince 4(5) explicit embedded adaptive time step Runge-Kutta method.
  - Error limits in Matlab "ode45" function set to 10<sup>-13</sup>.
- State Vector Components:
  - State Vector components 1, 2, 3 and 4 are the "Scalar First Unit Quaternion".
  - State Vector components 5, 6 and 7 are the angular velocity vector.
  - Initial Conditions:
    - Scalar First Unit Quaternion:  [ cos( 0.5 ); sin( 0.5 ); 0.0; 0.0 ].
      - This will result in a 1 [radian] (57.3 [degrees]) rotation about the body X axis.
    - Angular Velocity: [  0.10;  0.00001;  0.0 ] [radians/second]
      - Spin about body X axis with a small spin about the body Y axis.
    - Principal components of moments of inertia:  [  4.0;   1.0;      9.0 ]
      - Rectangular block dimensions determined by principal components of moments of inertia.
      - Spin will be about the rectangular block axis of intermediate length.
- Time rate of change of the Scalar First Unit Quaternion:
  - To determine the time rate of change of the Scalar First Unit Quaternion, **use a combination of the following two items**:
     - **The FIRST PART of Equation (46) on Page 266 of Reference [ 1 ]**.
       - **This states that the time rate of change of the attitude quaterion is equal to the cross product of two quaternions:**
         -  **The pure quaternion container of the angular velocity 3D vector.**
         -  **The attitude quaternion.**

    - **Equation (305) on Page 482 of Reference [ 2 ]**.
      - **This equation shows how to use the quaternion product when the quaternion cross product is specified.**
      - **Here the angular velocity is contained in a "pure" scalar first quaternion in the sense of Equation (182) on Page 466 of Reference [ 2 ]**
      - **(except that Equation (182) is scalar last).**
      - **Also note Equation (30) on Page 265 of Reference [ 1 ]**

- Note that the 2D plots of inertial angular momentum quantities reveal the constant nature of the inertial angular momentum over time as expected for free rotation.
  - Change is within 10<sup>-13</sup>.
- Note that the 2D plot of rotational kinetic energy reveals the constant nature of this energy over time as expected for free rotation.
  - Change is within 10<sup>-13</sup>.
- Note that the 3D animation video reveals the constant nature of the angular momentum vector over time as expected for free rotation.
  - Simulation generates a rotation vector solution every millisecond for 1000 seconds generating a total of one million rotation vector solutions.
  - Simulation creates a sequence of 3D animation frames, one 3D animation frame for every 200 th rotation vector solution.  3D simulation animation frames are at 5 Hertz rate.

CONTENTS:

- freeUnstableRigidBodyRotationDemonstration.m

  MATLAB source code file for top level simulation function.
- generateSolutionsForFreeRotationRigidBodyKinematics.m

  MATLAB source code file for non-real time solution of Ordinary Differential Equations (ODEs).
- generateInertialResults.m

  MATLAB source code file to convert state vector solutions at each time sample from body coordinates to inertial coordinates.
- generatePlots2D.m

  MATLAB source code file to generate 2D plots.
- generatePlots3D.m

  MATLAB source code file to generate 3D plots.
- performNonRealTimeAnimationPlayBack3D.m

  MATLAB source code file to generate 3D animation non-real time playback of solution applied to a rigid rectangular 3D block.
- generateFreeUnstableRigidBodyRotationDemoPurposeMessage.m
- generateFreeUnstableRigidBodyRotationDemoUsageMessage.m
- generateKinematicsSolutionsPurposeMessage.m
- generateKinematicsSolutionsUsageMessage.m

  Auxiliary MATLAB source code files.
- Utilities.

  Directory of utility MATLAB source code files.

EXAMPLES:
- nonRealTimeAnimationPartialPlayBackMovie.avi

  A non-real time 3D animation of a partial results play back AVI movie file.

- TwoDimensionalPlots

  Directory containing the two dimensional plots generated by this simulation.

REFERENCE(s):
    
    - [ 1 ]  "The Kinematic Equation for the Rotation Vector",      
              Malcolm D. Shuster,      
              IEEE Transactions on Aerospace and Electronic Systems,      
              Vol. 29, No. 1, pp. 263 - 267,      
              January 1993      
              https://malcolmdshuster.com/Pub_1993c_J_RotVec_IEEE.pdf
              

    - [ 2 ] "A Survey of Attitude Representations",
             Malcolm D. Shuster,      
             The Journal of the Astronautical Sciences,      
             Vol. 41, No. 4, pp. 430 - 517,      
             October-December, 1993      
             https://malcolmdshuster.com/Pub_1993h_J_Repsurv_scan.pdf
             

    - [ 3 ] "The Bizarre Behavior of Rotating Bodies -    - 
            The Dzhanibekov Effect",      
            Derek Alexander Muller      
            YouTube channel Veritasium      
            Spinning objects have strange instabilities known as      
            "The Dzhanibekov Effect" or "Tennis Racket Theorem" -      
            this video offers an intuitive explanation.      
            https://www.youtube.com/watch?v=1VPfZ_XzisU
            
      
    - [ 4 ] Elements of the following were used to help create this project:
      
        [ 4.1 ]  Gatech AE (2024). Euler Free Body Motion (https://github.com/alaricgregoire/EulerFreeBody/releases/tag/v1.0.0), GitHub. Retrieved July 15, 2024
        
        [ 4.2 ]   Ligong Han (2024). phymhan/matlab-axis-label-alignment (https://github.com/phymhan/matlab-axis-label-alignment), GitHub. Retrieved July 15, 2024. 
                 
        [ 4.3 ]   Georg Stillfried (2024). mArrow3.m - easy-to-use 3D arrow (https://www.mathworks.com/matlabcentral/fileexchange/25372-marrow3-m-easy-to-use-3d-arrow), MATLAB Central File Exchange. Retrieved July 15, 2024. 
      


    

