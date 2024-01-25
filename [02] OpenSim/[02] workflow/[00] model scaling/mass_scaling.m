clear all
clc
close all

%% Redefine mass properties in OpenSim based on scaling in BoB
% Massess < 1e-3 kg in OpenSim are maintained
% Total subject mass computed based on mass hand (direct correspondence in BoB)
% New mass based on scaling in BoB

%% Subject Specifications
TOT_MASS = 90;

%% Define Variables BOB

% BoB Fractions
bob_upper_arm_fraction = 2.5/100;
bob_forearm_fraction = 1.7/100;
bob_hand_fraction = 0.6/100;
bob_clavicle_fraction = 1e-4;
bob_scapula_fraction = 1e-3;
bob_humerus_fraction = 0.028;
bob_radius_fraction = 0.016;

%% Define Variables OpenSim

% OpenSim Masses
os_generic_tot_mass_model = 4.77882; % kg
os_generic_clavicle_mass = 0.156; % kg
os_generic_scapula_mass = 0.70396000000000003; %kg
os_generic_humerus_mass = 1.9975700000000001; % kg
os_generic_ulna_mass = 1.105299999999999; % kg
os_generic_radius_mass = 0.23358999999999999; % kg
os_generic_hand_mass = 0.58189999999999997; % kg

% OpenSim Inertias [Ixx Iyy Izz Ixy Ixz Iyz]
os_generic_clavicle_inertia = [0.00024258999999999999 0.00025525999999999999 4.4419999999999998e-05 -1.8980000000000001e-05 -6.9939999999999998e-05 5.3709999999999999e-05];
os_generic_scapula_inertia = [0.0012428999999999999 0.0011504 0.0013651 0.00044939999999999997 0.00040922000000000002 0.00024110000000000001];
os_generic_humerus_inertia = [0.0122776 0.0025513300000000001 0.0125789 -0.00034740999999999998 -0.00023250000000000001 0.0012293];
os_generic_radius_inertia = [0.00043855000000000001 8.8590000000000004e-05 0.00040257999999999998 3.0139999999999999e-05 -4.2400000000000001e-06 6.4179999999999999e-05];
os_generic_ulna_inertia = [0.0054130899999999997 0.00115318 0.0049436100000000002 0.00031686000000000003 -7.6149999999999994e-05 0.00109169];
os_generic_hand_inertia = [0.00011 6.0000000000000002e-05 0.00014999999999999999 8.9999999999999996e-07 -1.9999999999999999e-07 1.2e-05];

%% Compute Fractions Ulna/Radius from Single BoB Radius Fraction

% OpenSim ratio of radius/ulna mass
os_generic_radius_ulna_fraction = os_generic_radius_mass/os_generic_ulna_mass;

% Ulna & Radius seperate fractions
ulna_fraction = bob_radius_fraction/(1+os_generic_radius_ulna_fraction);
radius_fraction = ulna_fraction*os_generic_radius_ulna_fraction;

%% Compute OpenSim Total Subject Mass based on Hand Scaling

% total mass = ( hand mass ) / ( fraction BoB hand )
os_generic_tot_mass_subject = os_generic_hand_mass/bob_hand_fraction;

%% Compute OpenSim Bodies Mass based on BoB Fractions

% new mass = ( fraction BoB ) * ( total subject mass )
os_new_clavicle_mass = bob_clavicle_fraction*os_generic_tot_mass_subject;
os_new_scapula_mass = bob_scapula_fraction*os_generic_tot_mass_subject;
os_new_humerus_mass = bob_humerus_fraction*os_generic_tot_mass_subject;
os_new_hand_mass = bob_hand_fraction*os_generic_tot_mass_subject;
os_new_ulna_mass = ulna_fraction*os_generic_tot_mass_subject;
os_new_radius_mass = radius_fraction*os_generic_tot_mass_subject;

save("MOBL_ARMS41_adjusted_mass","os_new_clavicle_mass","os_new_scapula_mass",...
    "os_new_humerus_mass", "os_new_hand_mass","os_new_ulna_mass","os_new_radius_mass");

os_new_tot_mass_model = os_new_clavicle_mass + os_new_scapula_mass + ...
    os_new_humerus_mass + os_new_hand_mass + os_new_ulna_mass + os_new_radius_mass;

% fraction of subject mass relegated to OpenSim model
os_model_mass_fraction = os_new_tot_mass_model/os_generic_tot_mass_subject;

% required mass of OpenSim model when scaling
os_scaled_model_mass = os_model_mass_fraction*TOT_MASS;

%% Compute OpenSim Bodies Inertias based on New Computed Massess

% new inertia = ( new mass ) / ( old mass ) * ( old inertia )
os_new_clavicle_inertia = os_new_clavicle_mass/os_generic_clavicle_mass*os_generic_clavicle_inertia;
os_new_scapula_inertia = os_new_scapula_mass/os_generic_scapula_mass*os_generic_scapula_inertia;
os_new_humerus_inertia = os_new_humerus_mass/os_generic_humerus_mass*os_generic_humerus_inertia;
os_new_radius_inertia = os_new_radius_mass/os_generic_radius_mass*os_generic_radius_inertia;
os_new_ulna_inertia = os_new_ulna_mass/os_generic_ulna_mass*os_generic_ulna_inertia;
os_new_hand_inertia = os_new_hand_mass/os_generic_hand_mass*os_generic_hand_inertia;