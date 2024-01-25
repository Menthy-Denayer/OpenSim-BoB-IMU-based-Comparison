clear all
clc
close all

%% Define Variables

% General
OUT_FILE_NAME = "imu_orientations_sensors";
FRAME_RATE = 60;
N_IMU_MARKERS = 4;

% Inputs
EXAMPLE_FILE_LOC = "IMU_orientations_example.sto";
XSENS_FILE_LOC = ".../data.xlsx";

%% Load Data
example_file = importdata(EXAMPLE_FILE_LOC);
xsens_file = readtable(XSENS_FILE_LOC,"sheet","Sensor Orientation - Quat");
NFRAMES = size(xsens_file,1);

%% Create Dataset
data = zeros(NFRAMES,N_IMU_MARKERS*4+1);
k = 1;

% Time
data(:,k) = (0:1/60:(NFRAMES-1)/60)';
k=k+1;

% Thorax Sensor
data(:,k) = xsens_file.T8Q0;
k=k+1;
data(:,k) = xsens_file.T8Q1;
k=k+1;
data(:,k) = xsens_file.T8Q2;
k=k+1;
data(:,k) = xsens_file.T8Q3;
k=k+1;

% Right Shoulder
data(:,k) = xsens_file.RightShoulderQ0;
k=k+1;
data(:,k) = xsens_file.RightShoulderQ1;
k=k+1;
data(:,k) = xsens_file.RightShoulderQ2;
k=k+1;
data(:,k) = xsens_file.RightShoulderQ3;
k=k+1;

% Right Upper Arm
data(:,k) = xsens_file.RightUpperArmQ0;
k=k+1;
data(:,k) = xsens_file.RightUpperArmQ1;
k=k+1;
data(:,k) = xsens_file.RightUpperArmQ2;
k=k+1;
data(:,k) = xsens_file.RightUpperArmQ3;
k=k+1;

% Right Forearm
data(:,k) = xsens_file.RightForearmQ0;
k=k+1;
data(:,k) = xsens_file.RightForearmQ1;
k=k+1;
data(:,k) = xsens_file.RightForearmQ2;
k=k+1;
data(:,k) = xsens_file.RightForearmQ3;
k=k+1;

% Right Hand
% data(:,k) = xsens_file.RightHandQ0;
% k=k+1;
% data(:,k) = xsens_file.RightHandQ1;
% k=k+1;
% data(:,k) = xsens_file.RightHandQ2;
% k=k+1;
% data(:,k) = xsens_file.RightHandQ3;
% k=k+1;
% data(:,k) = xsens_file.RightForearmQ0;
% k=k+1;
% data(:,k) = xsens_file.RightForearmQ1;
% k=k+1;
% data(:,k) = xsens_file.RightForearmQ2;
% k=k+1;
% data(:,k) = xsens_file.RightForearmQ3;
% k=k+1;

%% Save New .sto File
fileID = fopen(strcat(OUT_FILE_NAME,'.sto'),'w');

% Print header info
for i = 1:5
    fprintf(fileID,string(example_file(i)) + "\n");
end

% fprintf(fileID,['time	thorax_imu	scapula_imu	humerus_imu	ulna_imu	hand_imu\n']);
fprintf(fileID,['time	thorax_imu	scapula_imu	humerus_imu	ulna_imu\n']);

% Print data
fprintf(fileID,['%5.4f\t' repmat('%5.4f,',1,3) '%5.4f\t' repmat('%5.4f,',1,3) '%5.4f\t' repmat('%5.4f,',1,3) '%5.4f\t' repmat('%5.4f,',1,3) '%5.4f','\n'],data');
fclose(fileID);