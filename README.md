# OpenSim-BoB-IMU-based-Comparison
Collection of model &amp; input files, scripts and methods used to perform a comparison between the OpenSim and Biomechanics of Bodies (BoB) simulation softwares.

Below, a breakdown is given of the steps required to replicate the performed experiments. __Pay attention to the following:__
- Be sure to add the 'Geometry' files required for the MoBL_ARMS_41 model to the same directory as the .osim model file
- Be sure to check file locations in order to succesfully run the provided scripts

# Biomechanics of Bodies (BoB) Workflow
To replicate the experiments inside of BoB, simply import the following files into BoB:
- XSens data file 'original_data.mvnx', containing IMU measurements
- Skeleton file 'skeleton_scaled.txt', containing the scaled skeleton
- Force file, 'bob_force.txt', containing the external force applied at the hand
- Change the model's height and mass to 1.86 m and 90 kg respectively

# OpenSim Workflow
To replicate the experiments inside of OpenSim, the following steps should be followed:
1. Segment &amp; mass scaling
Starting from the 'MoBL_ARMS_41_base.osim' model, use the scaling tool in OpenSim, combined with the results of the 'manual_scaling.m' script, to manually scale the segment sizes. Additionally, change the new model mass to the value computed in 'mass_scaling.m'.

Output: 'MoBL_ARMS_41_scaled.osim'

2. IMU calibration &amp; tracking
- use the 'save_imu_orientations.m' script to create IMU orientations file
- use the 'OpenSense_CalibrateModel.m' script to place the IMU's on the scaled model
- use the 'OpenSense_OrientationTracking.m' script to track the IMU orientations (IK)

Output: 'imu_orientations_sensors.sto', 'MoBL_ARMS_41_IMU_scaled.osim', 'IK_imu_scaled_segments.mot'

3. Inverse dynamics
- use the 'OpenSim_ExternalLoads' script to create the required force file
- use the post-processing 'MoBL_Arms_IKfilter.m' script to filter the found IK
- use the OpenSim inverse dynamics tool to compute joint torques. Use the filtered inverse kinematics, check the filter kinematics at 6 Hz, add the external forces.

Output: 'opensim_force.txt', 'opensim_force.xml', 'filteredIK_imu_scaled_segments.sto', 'ID_imu_scaled_segments.sto'

4. Static optimization
- use the static optimization tool in OpenSim to compute muscle forces and activations
- add the residual actuators from the 'Actuators_3dof.xml' file
- add the external loads from the 'opensim_force.xml' file
- filter the kinematics at 6 Hz
- change the precision to 20

Output: 'SO_imu_muscle_force_scaled_segments.sto', 'SO_imu_muscle_activations_scaled_segments.sto', 'SO_imu_settings_scaled_segments.xml'

# Post-processing Workflow
As the definition of joint angles differs between OpenSim and BoB, the OpenSim shoulder angles must be transformed into equivalent XYZ Euler angles. This is achieved by using the 'OpenSim_orientations.m' script.

Before applying the static optimization tool inside OpenSim, it is best to filter the results of the inverse kinematics tool. This is achieved using the 'MoBL_Arms_IKfilter.m' script, provided by the MoBL Arms project.

Plots of the results and computation of the correlation coefficients and errors can be repeated using the 'plotter.m' script.

# Acknowledgements
- The OpenSim model and IK filtering code were taken from the [MoBL Arms project](https://simtk.org/projects/upexdyn)
- The IMU calibration and tracking codes were taken from the [OpenSense project](https://simtk.org/projects/opensense)
