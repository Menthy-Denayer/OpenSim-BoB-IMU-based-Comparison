clear all
clc 
close all

%% OpenSim External Loads

% General
GRAVITY = 9.81; % m/s²
OUT_FILE_NAME = "opensim_force";
TEMPLATE_FILE_NAME = "OpenSim_ForceTemplate.mot";

% Experiment
WEIGHT = 3; % kg
TIME_INT = [421,541,901,1021]; % indices between which to apply force

% OpenSim results (IK results)
OPENSIM_RESULTS_FOLDER = "...";
OPENSIM_RESULTS_FILE = OPENSIM_RESULTS_FOLDER + "/" + "IK_imu_scaled_segments.mot";

%% Read Time Vector
time = read_kinematic_data("time",OPENSIM_RESULTS_FILE);

%% Define Force Vectors (absolute axes)

F_mass = [0;-WEIGHT*GRAVITY;0];
NDATA = length(time);

F_mass_rel = zeros(NDATA,3);
for i = 1:NDATA
    if( (TIME_INT(1) <= i && i <= TIME_INT(2)) || (TIME_INT(3) <= i && i <= TIME_INT(4)))
         F_mass_rel(i,:) = [0,-WEIGHT*GRAVITY,0];
    else
        F_mass_rel(i,:) = zeros(1,3);
    end
end

%% Load Template GRF File
template_file = importdata(TEMPLATE_FILE_NAME);

%% Create New Forces
NDATA = size(time,1);
k = 1;
new_motion_file.data(:,k) = time;        % Time
k=k+1;

% Mass
new_motion_file.data(:,k) = F_mass_rel(:,1);             % Force_x
k=k+1;
new_motion_file.data(:,k) = F_mass_rel(:,2);             % Force_y
k=k+1;
new_motion_file.data(:,k) = F_mass_rel(:,3);             % Force_z
k=k+1;
new_motion_file.data(:,k) = zeros(NDATA,1);             % Position_x
k=k+1;
new_motion_file.data(:,k) = zeros(NDATA,1);             % Position_y
k=k+1;
new_motion_file.data(:,k) = zeros(NDATA,1);             % Position_z
k=k+1;

NCOLUMNS = size(new_motion_file.data,2);

%% Create Header
% Copy template header
new_motion_file.colheaders = template_file.colheaders(:,1:7);
new_motion_file.textdata = template_file.textdata(:,1:7);

% Adjust headers
new_motion_file.textdata(4,1) = {strcat('nColumns=',string(NCOLUMNS))};    
new_motion_file.textdata(1,1) = {strcat(OUT_FILE_NAME,'.mot')};  

%% Create Force File (.mot)
fileID = fopen(strcat(OUT_FILE_NAME,'.mot'),'w');

% Print header info
fprintf(fileID,string(new_motion_file.textdata(1,1)) + "\n");
fprintf(fileID,string(new_motion_file.textdata(2,1)) + "\n");
fprintf(fileID,string(new_motion_file.textdata(3,1)) + "\n");
fprintf(fileID,string(new_motion_file.textdata(4,1)) + "\n");
fprintf(fileID,string(new_motion_file.textdata(5,1)) + "\n");
fprintf(fileID,string(new_motion_file.textdata(6,1)) + "\n");
fprintf(fileID,['time  hand_force_vx  hand_force_vy  hand_force_vz hand_force_px  hand_force_py  hand_force_pz\n']);

% Print data
data = new_motion_file.data';
fprintf(fileID,[repmat('%5.4f\t',1,size(data,1)),'\n'],data);
fclose(fileID);

%% Functions

function [data_list,time_list,range] = read_kinematic_data(column_name,file_loc)
% sto, mot, csv (physio angles BoB)
    results = importdata(file_loc);
    data_index = find(results.colheaders == column_name);
    if(isempty(data_index))
        disp(results.colheaders)
        error("Column name should be inside this list!")
    end
    data_list = results.data(:,data_index);
    time_list = results.data(:,1);
    range = max(data_list)-min(data_list);
end