clear all
clc 
close all

%% Define Variables

% OpenSim IK results
OPENSIM_RESULTS_FOLDER = "...";
OPENSIM_RESULTS_FILE = OPENSIM_RESULTS_FOLDER + "/" + "IK_markers_non_fixed_thorax.mot";
os_txtdata = "colheaders";
delimiter = "	";

% Define Output Location
output_file = OPENSIM_RESULTS_FOLDER + "/" + "IK_os_markers_non_fixed_thorax.csv";

% Rotation axis elevation angle
rot_axis1 = [0.0048; 0.999089; 0.0424];
% rot_axis1 = [0;1;0];

% Rotation axis shoulder elevation
rot_axis2 = [-0.998261; 0.0023; 0.058898];
% rot_axis2 = [-1;0;0];

% Rotation axis shoulder1_r2 (= opposite to elevation angle)
rot_axis3 = [0.0048; 0.999089; 0.0424];
% rot_axis3 = [0;1;0];

% Rotation axis shoulder rotation
rot_axis4 = [0.0048; 0.999089; 0.0424];
% rot_axis4 = [0;1;0];


%% Define Angles OpenSim

% Read OpenSim angles
os_time = read_kinematic_data("time",OPENSIM_RESULTS_FILE,os_txtdata,delimiter);
NDATA = length(os_time);

os_a1 = read_kinematic_data("elv_angle",OPENSIM_RESULTS_FILE,os_txtdata,delimiter);
% os_a1 = 90*ones(NDATA,1);

os_a2 = read_kinematic_data("shoulder_elv",OPENSIM_RESULTS_FILE,os_txtdata,delimiter);
% os_a2 = 0:90/(NDATA-1):90;

os_a3 = read_kinematic_data("shoulder1_r2",OPENSIM_RESULTS_FILE,os_txtdata,delimiter);
% os_a3 = -os_a1;

os_a4 = read_kinematic_data("shoulder_rot",OPENSIM_RESULTS_FILE,os_txtdata,delimiter);
% os_a4 = 0*ones(NDATA,1);

%% Compute Total Rotation Matrix
% Opensim Reference Axes
os_x_vector = [1;0;0];
os_y_vector = [0;1;0];
os_z_vector = [0;0;1];

% Euler XYZ in BoB
new_ax_list = zeros(NDATA,1);
new_ay_list = zeros(NDATA,1);
new_az_list = zeros(NDATA,1);

for t = 1:NDATA

    R = opensim_shoulder_rot(os_a1(t),os_a2(t),os_a3(t),os_a4(t),rot_axis1,rot_axis2,rot_axis3,rot_axis4);
    
    % Humerus Axes (OpenSim)
    os_x_new = R*os_x_vector;
    os_y_new = R*os_y_vector;
    os_z_new = R*os_z_vector;

    % Opensim --> BoB Humerus Axes
    Rcorr1 = rotation_matrix(os_x_new,pi/2);
    Rcorr2 = rotation_matrix(Rcorr1*os_z_new,pi);
    Rcorr = Rcorr2*Rcorr1;

    % BoB reference --> BoB Humerus Transformation
    Rtot = Rcorr*R*(rotation_matrix(os_z_vector,pi)');

    % compute z-y-x euler angles

    % WEIRD?? Changing definition all the time?
    [euler_angles, euler_alt] = rotm2eul(Rtot,"XYZ");
    new_ax_list(t) = -euler_angles(1);
    new_ay_list(t) = -euler_angles(2);
    new_az_list(t) = euler_angles(3);
end

%% Plot Results

figure 
hold on
grid on
plot(os_time, new_ax_list*180/pi,'Color','blue','LineWidth',2)
plot(os_time, new_ay_list*180/pi,'Color','red','LineWidth',2)
plot(os_time, new_az_list*180/pi,'Color','green','LineWidth',2)
legend(["Rot1","Rot2","Rot3"])
hold off


%% Read Elbow Rotations OpenSim
elbow_flex = read_kinematic_data("elbow_flexion",OPENSIM_RESULTS_FILE,os_txtdata,delimiter);
elbow_rot = -read_kinematic_data("pro_sup",OPENSIM_RESULTS_FILE,os_txtdata,delimiter);

%% Save Extracted Euler Angles for BoB

file_name = "BoB_motion.txt";
data_shoulder = [os_time,new_ax_list*180/pi,new_ay_list*180/pi,new_az_list*180/pi];
data_elbow = [os_time,elbow_flex,elbow_rot];
segment_names = ["right_shoulder","right_elbow"];

save_txt_motion(file_name, data_shoulder, data_elbow, segment_names)

%% Save Extracted Euler Angles for Post-Processing

results_table = table(os_time, new_ax_list,new_ay_list,new_az_list,elbow_flex,elbow_rot);
writetable(results_table,output_file)

%% Plot Results over time

% Opensim Reference Axes
os_x_vector = [1;0;0];
os_y_vector = [0;1;0];
os_z_vector = [0;0;1];

% BoB Reference Axes
Rcorr = rotation_matrix(os_z_vector,pi);
bob_x_vector = Rcorr*[1;0;0];
bob_y_vector = Rcorr*[0;1;0];
bob_z_vector = Rcorr*[0;0;1];

fig = figure;

for t = 1 
    R = opensim_shoulder_rot(os_a1(t),os_a2(t),os_a3(t),os_a4(t),rot_axis1,rot_axis2,rot_axis3,rot_axis4);

    os_x_new = R*os_x_vector;
    os_y_new = R*os_y_vector;
    os_z_new = R*os_z_vector;

    Rcorr1 = rotation_matrix(os_x_new,pi/2);
    Rcorr2 = rotation_matrix(Rcorr1*os_z_new,pi);
    Rcorr = Rcorr2*Rcorr1;

    % bob_x_new = Rcorr*os_x_new;
    % bob_y_new = Rcorr*os_y_new;
    % bob_z_new = Rcorr*os_z_new;

    Rtot = Rcorr*R*(rotation_matrix(os_z_vector,pi)');
    bob_x_new = Rtot*bob_x_vector;
    bob_y_new = Rtot*bob_y_vector;
    bob_z_new = Rtot*bob_z_vector;

    clf(fig)
    view(45,45)
    title("Time is " + round(os_time(t),6,'significant'))
    xlim([-0.5,1.1]);
    ylim([-0.5,1.1]);
    zlim([-0.5,1.1]);
    xlabel("x")
    ylabel("y")
    zlabel("z")
    hold on

    % OpenSim
    plot3([0 0],[0 1],[0 0],'LineWidth',2,'Color','green')
    plot3([0 1],[0 0],[0 0],'LineWidth',2,'Color','red')
    plot3([0 0],[0 0],[0 1],'LineWidth',2,'Color','blue')
    plot3([0 os_x_new(1)],[0 os_x_new(2)],[0 os_x_new(3)],'Color',[0.5 0 0])
    plot3([0 os_y_new(1)],[0 os_y_new(2)],[0 os_y_new(3)],'Color',[0 0.5 0])
    plot3([0 os_z_new(1)],[0 os_z_new(2)],[0 os_z_new(3)],'Color',[0 0 0.5])

    % % BoB
    % plot3([0 bob_x_vector(1)],[0 bob_x_vector(2)],[0 bob_x_vector(3)],'LineWidth',2,'Color','blue','LineStyle','--')
    % plot3([0 bob_y_vector(1)],[0 bob_y_vector(2)],[0 bob_y_vector(3)],'LineWidth',2,'Color','green','LineStyle','--')
    % plot3([0 bob_z_vector(1)],[0 bob_z_vector(2)],[0 bob_z_vector(3)],'LineWidth',2,'Color','red','LineStyle','--')
    % plot3([0 bob_x_new(1)],[0 bob_x_new(2)],[0 bob_x_new(3)],'Color',[0.5 0.5 0.8],'LineStyle','--')
    % plot3([0 bob_y_new(1)],[0 bob_y_new(2)],[0 bob_y_new(3)],'Color',[0.5 0.8 0.5],'LineStyle','--')
    % plot3([0 bob_z_new(1)],[0 bob_z_new(2)],[0 bob_z_new(3)],'Color',[0.8 0.5 0.5],'LineStyle','--')
    hold off

    drawnow()
end


%% Functions

function R = opensim_shoulder_rot(os_a1,os_a2,os_a3,os_a4,rot_axis1,rot_axis2,rot_axis3,rot_axis4)

    % R1 rotation shoulder0 parent --> shoulder0 child
    R1 = rotation_matrix(rot_axis1,os_a1*pi/180);

    % R2 rotation shoulder1 parent --> shoulder1 child
    R2 = rotation_matrix(rot_axis2,os_a2*pi/180);
    
    % R3 rotation shoulder1 parent --> shoulder1 child
    R3 = rotation_matrix(rot_axis3,os_a3*pi/180);

    % R4 rotation shoulder2 parent --> shoulder2 child
    R4 = rotation_matrix(rot_axis4,os_a4*pi/180);
    
    % R rotation shoulder0 parent --> humerus frame
    R = R1*R2*R3*R4;
end

function R = rotation_matrix(r,angle)
% Compute rotation matrix around x/y/z axis over given angle
% Parameters
% ----------
% axis          : int       : 0/1/2 specifying axis of rotation (x/y/z)
% angle         : float     : rotation angle
%
% return
% R             : 3x3 array : rotation matrix
    
    S = antisymmetric(r);
    % angle = angle*-1;
    R = r*r' + (eye(3) - r*r')*cos(angle) + S*sin(angle);
end

function S = antisymmetric(r)
    S = [0 -r(3) r(2);r(3) 0 -r(1);-r(2) r(1) 0];
end


function [data_list,time_list,range] = read_kinematic_data(column_name,file_loc,txt_data,delimiter)
% sto, mot, csv (physio angles BoB)

    results = importdata(file_loc);

    if(txt_data == "colheaders")
        data_names = results.colheaders;
    elseif(txt_data == "txtdata") 
        data_names = results.textdata(end);
        data_names = string(data_names);
        data_names = strsplit(data_names,delimiter);
    else
        error("Unknown data format!")
    end

    data_index = find(data_names == column_name);
    if(isempty(data_index))
        disp(data_names)
        error("Column name should be inside this list!")
    end
    data_list = results.data(:,data_index);
    time_list = results.data(:,1);
    range = max(data_list)-min(data_list);
end


function save_txt_motion(file_name, data_shoulder, data_elbow, segment_names)
% Save data into .txt file for BoB
% Parameters
% ----------
% file_name         : string            : name of output file
% data              : Nx2 table         : time & torque to save 
% segment_name      : string            : BoB segment name to apply torque

    % Create new file
    output_file = fopen(file_name,'w');
    
    nmb_lines_tot = size(data_shoulder,2)+2;    % Total number of lines
    newLines = cell(nmb_lines_tot,1);    % Store new created lines
    
    % Right arm
    k=1;
    newLines{k} = char("%% OpenSim Extracted Euler Angles");
    k=k+1;
    newLines{k} = char("% Right Shoulder");
    k=k+1;
    newLines{k} = create_line(segment_names(1) + ".time",data_shoulder(:,1));
    k=k+1;
    newLines{k} = create_line(segment_names(1) + ".rotx",data_shoulder(:,2));
    k=k+1;
    newLines{k} = create_line(segment_names(1) + ".roty",data_shoulder(:,3));
    k=k+1;
    newLines{k} = create_line(segment_names(1) + ".rotz",data_shoulder(:,4));
    k=k+1;

    newLines{k} = char("% Right Elbow");
    k=k+1;
    newLines{k} = create_line(segment_names(2) + ".time",data_shoulder(:,1));
    k=k+1;
    newLines{k} = create_line(segment_names(2) + ".roty",data_elbow(:,2));
    k=k+1;
    newLines{k} = create_line(segment_names(2) + ".rotz",data_elbow(:,3));
   
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