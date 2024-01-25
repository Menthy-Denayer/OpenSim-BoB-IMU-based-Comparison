clear all
clc 
close all

%% Define XSENS Measurements
arm_span = 1.86; % m
shoulder_width = 0.44; % m
wrist_span = 1.46; % m
elbow_span = 0.96; % m

%% Define OpenSim Marker Locations
% Marker locations in ground frame

% Position of acromion in N-pose
right_acromion_marker = [-0.197467 0.025406 -0.0241164]; % m
left_acromion_marker = right_acromion_marker; left_acromion_marker(1)=left_acromion_marker(1)*-1;

% Position of elbow in T-pose
right_elbow_marker = [-0.445283 0.00185013 -0.0460539]; % m
left_elbow_marker = right_elbow_marker; left_elbow_marker(1)=left_elbow_marker(1)*-1; % m

% Position of wrist in T-pose
right_wrist_marker = [-0.698265 0.0461035 0.00837696]; % m
left_wrist_marker = right_wrist_marker; left_wrist_marker(1)=left_wrist_marker(1)*-1; % m

% Position of knuckles low in T-pose
hand1_marker = [-0.783496 0.0639017 0.0205733]; % m

% Position of knuckles up in T-pose
hand2_marker = [-0.804502 0.0279963 0.0334607]; % m

% Position of finger ends in T-pose
hand3_marker = [-0.787179 -0.00614436 0.0454131]; % m

%% Compute Scaling Factors

os_shoulder_width = norm(right_acromion_marker-left_acromion_marker);
os_elbow_span = norm(right_elbow_marker-left_elbow_marker);
os_wrist_span = norm(right_wrist_marker-left_wrist_marker);
os_hand_length = norm(hand1_marker-right_wrist_marker) + ...
    norm(hand2_marker-hand1_marker) + norm(hand3_marker-hand2_marker);

shoulder_scaling = (shoulder_width)/os_shoulder_width;
humerus_scaling = (elbow_span)/os_elbow_span;
forearm_scaling = (wrist_span)/os_wrist_span;
hand_scaling = ((arm_span-wrist_span)/2)/os_hand_length;

%% Print Results
fprintf("Scale scapula, clavicle and thorax by " + shoulder_scaling +"\n")
fprintf("Scale humerus by " + humerus_scaling +"\n")
fprintf("Scale radius and ulna " + forearm_scaling +"\n")
fprintf("Scale hand by " + hand_scaling +"\n")