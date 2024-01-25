clear all
clc
close all

%% Define Variables

% General
GRAVITY = 9.81; % m/s²
BOB_TIMESTEP = 0.1; % s
OUT_FILE_NAME = "bob_force.txt";

% Experiment
WEIGTH = 3; % kg
TRIAL_END = 35.4; % s
BOB_SEGMENT_NAMES = ["right_hand"];

% Create data vectors
time = 0:BOB_TIMESTEP:TRIAL_END;
NDATA = length(time);
force = zeros(NDATA,1);
force(71:91) = -WEIGTH*GRAVITY;
force(151:171) = -WEIGTH*GRAVITY;

% Save to .txt file
save_txt_force(OUT_FILE_NAME,time',force,BOB_SEGMENT_NAMES);

%% Functions

% Saving Data 
function save_txt_force(file_name, time, force, segment_names)
% Save data into .txt file for BoB
% Parameters
% ----------
% file_name         : string            : name of output file
% data              : Nx2 table         : time & torque to save 
% segment_name      : string            : BoB segment name to apply torque

    % Create new file
    output_file = fopen(file_name,'w');
    
    nmb_lines_tot = size(time,2)+3;    % Total number of lines
    newLines = cell(nmb_lines_tot,1);    % Store new created lines
    
    % Right hand
    k=1;
    newLines{k} = char("%% Lumped mass");
    k=k+1;
    newLines{k} = create_line(segment_names(1) + ".force.time",time);
    k=k+1;
    newLines{k} = create_line(segment_names(1) + ".force.type",3);
    k=k+1;
    newLines{k} = create_line(segment_names(1) + ".force.fz",force);
   
    fwrite(output_file,strjoin(newLines, '\n'));
    fclose(output_file);

end

function[new_line]  = create_line(variable_name,data_array)
% Create data line in BoB format Parameters ---------- variable_name :
% string data_array    : Nx1 array
%
% return
% new char line : 1x1 char

    str = mat2str(data_array');
    new_line = char(variable_name + ' = ' + string(str) + ';');
end