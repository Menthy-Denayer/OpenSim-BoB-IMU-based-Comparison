clear all
clc
close all

%% Define Variables 

% BoB results location
BOB_RESULTS_FOLDER = "...";

% OpenSim results location
OPENSIM_RESULTS_FOLDER = "...";

%% Plot OpenSim Results

% Choose OpenSim File to Read
opensim_results_file = OPENSIM_RESULTS_FOLDER + "/" + "filteredIK_imu.sto";
opensim_column_name = "r_x";

plot_ylabel = "Torque [Nm]";
plot_name = "rx";

% Plot OpenSim Results
plot_data_opensim(opensim_results_file,opensim_column_name,plot_ylabel,plot_name)

%% Plot BoB Results

% Choose BoB File to Read
bob_results_file = BOB_RESULTS_FOLDER + "/" + "right_shoulder_torque.csv";
bob_column_name = "Torque magnitude[Nm]";
text_data = "textdata"; % specify where column headers are saved
data_delimiter = ","; % specify delimiter seperating data in file

plot_ylabel = "Torque [Nm]";
plot_name = "Shoulder Torque Magnitude";

% Plot BoB Results
plot_data_bob(bob_results_file, bob_column_name, plot_ylabel, plot_name, text_data, data_delimiter)

%% Compare Shoulder IK Results

% Plot generic & scaled segments in same graph?

% OpenSim IK result file (.csv)
opensim_ik_results_file = OPENSIM_RESULTS_FOLDER + "/" + "IK_os_imu_scaled_segments.csv";
opensim_ik_columns = ["new_ax_list","new_ay_list","new_az_list"];
os_factor = 180/pi;

% BoB IK result file (.csv)
bob_ik_results_file = BOB_RESULTS_FOLDER + "/" + "right_shoulder_angle_euler_scaled_segments.csv";
bob_ik_columns = ["Rot1","Rot2","Rot3"];
bob_txtdata = "textdata";
bob_delimiter = " ";

% Plot IK results for BoB - OpenSim
plot_name = "Shoulder Euler Angles";
compare_IK_data(opensim_ik_results_file,opensim_ik_columns,os_factor,bob_ik_results_file, bob_ik_columns, bob_txtdata, bob_delimiter, plot_name)

%% Compare Elbow IK Results

% Plot generic & scaled segments in same graph?

% OpenSim IK result file (.csv)
opensim_ik_results_file = OPENSIM_RESULTS_FOLDER + "/" + "IK_os_imu_scaled_segments.csv";
opensim_ik_columns = ["elbow_flex","elbow_rot"];
os_factor = 1;

% BoB IK result file (.csv)
bob_ik_results_file = BOB_RESULTS_FOLDER + "/" + "all_joints_physio_angles_scaled_segments.csv";
bob_ik_columns = ["Right elbow Flex/Ext","Right elbow Int/Ext"];
bob_txtdata = "colheaders";
bob_delimiter = "";

% Plot IK results for BoB - OpenSim
plot_name = "Elbow Angles";
compare_IK_data(opensim_ik_results_file, opensim_ik_columns, os_factor, bob_ik_results_file, bob_ik_columns, bob_txtdata, bob_delimiter, plot_name)

%% Compare Torque Magnitude Shoulder

% BoB ID result file (.csv)
bob_results_file = BOB_RESULTS_FOLDER + "/" + "right_shoulder_torque_scaled_segments.csv";

% OpenSim ID result files (.sto)
opensim_results_file1 = OPENSIM_RESULTS_FOLDER + "/" + "ID_results_fixed_thorax.sto";
opensim_results_file2 = OPENSIM_RESULTS_FOLDER + "/" + "ID_results_fixed_thorax_adjusted_mass.sto";
opensim_files = [opensim_results_file1,opensim_results_file2];

opensim_files = [OPENSIM_RESULTS_FOLDER+"/ID_imu_scaled_segments.sto"];

% OpenSim ID result column names
columns = ["elv_angle_moment","shoulder_elv_moment","shoulder_rot_moment"];

% Compare torque magnitude to BoB result
plot_name = "Shoulder Torque Magnitude";
% plot_labels = ["OpenSim Generic","OpenSim Adjusted","BoB"];
plot_labels = ["OpenSim","BoB"];
compare_torque_magn(opensim_files,columns,bob_results_file,plot_name,plot_labels,false)

%% Compare Torque Magnitude Elbow

% BoB ID result file (.csv)
bob_results_file = BOB_RESULTS_FOLDER + "/" + "right_elbow_torque_scaled_segments.csv";

% OpenSim ID result files (.sto)
opensim_results_file1 = OPENSIM_RESULTS_FOLDER + "/" + "ID_results_fixed_thorax.sto";
opensim_results_file2 = OPENSIM_RESULTS_FOLDER + "/" + "ID_results_fixed_thorax_adjusted_mass.sto";
opensim_files = [opensim_results_file1,opensim_results_file2];

opensim_files = [OPENSIM_RESULTS_FOLDER+"/ID_imu_scaled_segments.sto"];

% OpenSim ID result column names
columns = ["elbow_flexion_moment","pro_sup_moment"];

% Compare torque magnitude to BoB result
plot_name = "Elbow Torque Magnitude";
% plot_labels = ["OpenSim Generic","OpenSim Adjusted","BoB"];
plot_labels = ["OpenSim","BoB"];
compare_torque_magn(opensim_files,columns,bob_results_file,plot_name,plot_labels,false)

%% Compare Muscle Forces
% BoB SO result files (.csv)
bob_results_file1 = BOB_RESULTS_FOLDER + "/" + "all_muscle_forces_generic_segments.csv";
bob_results_file2 = BOB_RESULTS_FOLDER + "/" + "all_muscle_forces_scaled_segments.csv";

bob_result_files = [bob_results_file1,bob_results_file2];
bob_result_files = [bob_results_file2];
bob_columns = ["Triceps brachii long head right",...
    "Triceps brachii medial head right","Triceps brachii lateral head right"];

% OpenSim SO result file (.sto)
opensim_results_file = OPENSIM_RESULTS_FOLDER + "/" + "SO_imu_muscle_force_scaled_segments.sto";
opensim_columns = ["TRIlong","TRImed","TRIlat"];

% Compare muscle forces between BoB and OpenSim
os_error = 29.7444;
plot_name = "Triceps Force";
plot_labels = ["OpenSim (long)","OpenSim (medial)","OpenSim (lateral)","BoB (long)","BoB (medial)","BoB (lateral)"];
compare_muscle_forces(bob_result_files,bob_columns,opensim_results_file,opensim_columns,plot_name,plot_labels,os_error)

%% Compare Muscle Activations
% BoB SO result files (.csv)
bob_results_file1 = BOB_RESULTS_FOLDER + "/" + "all_muscle_forces_generic_segments.csv";
bob_results_file2 = BOB_RESULTS_FOLDER + "/" + "all_muscle_forces_scaled_segments.csv";

bob_result_files = [bob_results_file1,bob_results_file2];
bob_result_files = [bob_results_file2];
bob_columns = ["Triceps brachii long head right",...
    "Triceps brachii medial head right","Triceps brachii lateral head right"];

% Triceps Max Isometric Force
bob_max_force = [252.9,407.3,215.8]; % N

% OpenSim SO result file (.sto)
opensim_results_file = OPENSIM_RESULTS_FOLDER + "/" + "SO_imu_muscle_activations_scaled_segments.sto";
opensim_columns = ["TRIlong","TRImed","TRIlat"];

% Compare muscle forces between BoB and OpenSim
os_error = 29.7444;
plot_name = "Triceps Brachii Activations";
plot_labels = ["OpenSim (long)","OpenSim (medial)","OpenSim (lateral)","BoB (long)","BoB (medial)","BoB (lateral)"];
compare_muscle_activations(bob_result_files,bob_columns,bob_max_force,opensim_results_file,opensim_columns,plot_name,plot_labels,os_error)

%% Compute RMS SO Residuals
opensim_result_file = OPENSIM_RESULTS_FOLDER + "/" + "SO_imu_muscle_force_scaled_segments.sto";
opensim_column = "rz";
os_error = 29.7444;

rms_val = compute_rms_curve(opensim_result_file, opensim_column, os_error, "colheaders", "");

fprintf("RMS value for " + opensim_column + " is: " + rms_val +"\n")

%% Functions

% Plot OpenSim Data
function plot_data_opensim(file_loc, column_name, ylabel_txt, title_txt)
    [data_list, time_list] = read_data(column_name, file_loc, "colheaders", "");

    plot_comparison({time_list},{data_list},title_txt,ylabel_txt,'','','');
end

% Plot BoB Data
function plot_data_bob(file_loc, column_name, ylabel_txt, title_txt, textdata, delimiter)
    [data_list, time_list] = read_data(column_name, file_loc, textdata, delimiter);
    plot_comparison(time_list,data_list,title_txt,ylabel_txt,'','','');
end

% Compare BoB - OpenSim results and plot them
function compare_data(opensim_file, opensim_column, bob_file, bob_column, ylabel_txt, title_txt)
    [os_data_list, os_time_list,~] = read_data(opensim_column, opensim_file);
    [bob_data_list, bob_time_list,~] = read_data(bob_column, bob_file);
    
    plot_comparison({os_time_list,bob_time_list},{os_data_list,bob_data_list},title_txt,ylabel_txt,["red","blue"],'','')
end

% Compare BoB - OpenSim results for IK and plot them
function compare_IK_data(os_file, os_columns, os_factor, bob_file, bob_columns, bob_txtdata, bob_delimiter, title_txt)

    N_os_columns = size(os_columns,2);
    N_bob_columns = size(bob_columns,2);
    
    Ndata = N_bob_columns+N_os_columns;
    tot_data_list = cell(Ndata);
    tot_time_list = cell(Ndata);
    k=1;

    for i = 1:N_os_columns
        [os_data_list, os_time_list,~] = read_data(os_columns(i), os_file,"colheaders","");
        tot_time_list{k} = os_time_list;
        tot_data_list{k} = os_data_list*os_factor;
        k = k + 1;
    end

    for i = 1:N_bob_columns
        [bob_data_list, bob_time_list,~] = read_data(bob_columns(i), bob_file,bob_txtdata,bob_delimiter);
        
        if(bob_columns(i)=="Right elbow Int/Ext")
            bob_data_list(bob_data_list>0) = bob_data_list(bob_data_list>0)-360;
            bob_data_list = bob_data_list-bob_data_list(1);
        end
        
        tot_time_list{k} = bob_time_list;
        tot_data_list{k} = bob_data_list;
        k = k + 1;
    end

    % Compute correlation coefficient
    for i = 1:N_os_columns
        os_data = tot_data_list{i};
        bob_data = tot_data_list{i+N_os_columns};
        os_Ndata = length(os_data);
        bob_Ndata = length(bob_data);
        Nsample = os_Ndata/bob_Ndata;
        os_data_sampled = os_data(1:Nsample:end);
        fprintf("Computing Pearson Correlation Coefficient for " + os_columns(i) + "\n")
        r = pearson_corr_coeff(bob_data,os_data_sampled);
        fprintf("Computed Correlation Coefficient for " + os_columns(i) + " is " + r + "\n")

        fprintf("Computing RMSE for " + os_columns(i) + "\n")
        rmse_value = rmse(bob_data,os_data_sampled);
        fprintf("Computed RMSE for " + os_columns(i) + " is " + rmse_value + "\n")
    end
    
    if(N_os_columns == 3)
        colors = [[1,0,0];[0,1,0];[0,0,1];[1,0,0];[0,1,0];[0,0,1]];
        style = ["--","--","--","-","-","-"];
        labels = ["Eulerx (OpenSim)","Eulery (OpenSim)","Eulerz (OpenSim)",...
        "Eulerx (BoB)","Eulery (BoB)","Eulerz (BoB)"];
    else
        colors = [[1,0,0];[0,0,1];[1,0,0];[0,0,1]];
        style = ["--","--","-","-"];
        labels = ["Flex/Ext (OpenSim)","Ext/Int (OpenSim)","Flex/Ext (BoB)","Ext/Int (BoB)"];
    end
    
    plot_comparison(tot_time_list,tot_data_list,title_txt,"Angle [°]",colors,style,labels);
end

% Compare BoB - OpenSim Joint Torque Magnitudes and plot them
function compare_torque_magn(opensim_files,opensim_columns,bob_file,title_txt,label_txt,remove_offset)
    
    Ndata = size(opensim_files,2);
    tot_data_list = cell(Ndata+1,1);
    tot_time_list = cell(Ndata+1,1);

    [bob_data_list, bob_time_list,~] = read_data("Torque magnitude[Nm]", bob_file,"textdata",",");
    
    bob_Ndata = length(bob_data_list);
    
    tot_data_list{Ndata+1} = bob_data_list;
    tot_time_list{Ndata+1} = bob_time_list;

    for i = 1:Ndata
        [os_data_list, os_time_list] = compute_magnitude(opensim_columns,opensim_files(i),remove_offset);

        os_Ndata = length(os_data_list);
        Nsample = os_Ndata/bob_Ndata;
        os_data_list_sampled = os_data_list(1:Nsample:end);
        fprintf("Computing Pearson Correlation Coefficient for " + label_txt(i+1) + "\n")
        r = pearson_corr_coeff(bob_data_list,os_data_list_sampled);
        fprintf("Computed Correlation Coefficient for " + label_txt(i+1) + " is " + r + "\n")

        fprintf("Computing RMSE for " + label_txt(i+1) + "\n")
        rmse_value = rmse(bob_data_list,os_data_list_sampled);
        fprintf("Computed RMSE for " + label_txt(i+1) + " is " + rmse_value + "\n")
        
        tot_data_list{i} = os_data_list;
        tot_time_list{i} = os_time_list;
    end
    
    colors = [[0, 0, 0.73];[0.73, 0, 0]];

    plot_comparison(tot_time_list,tot_data_list,...
        title_txt,"Torque [Nm]",colors,'',label_txt)

end

% Compare muscle forces
function compare_muscle_forces(bob_files,bob_columns,opensim_file,opensim_columns,title_txt,label_txt,os_err)
    
    % Determine data size
    N_os_columns = size(opensim_columns,2);
    N_bob_files = size(bob_files,2);
    N_bob_columns = size(bob_columns,2);
    Ndata = N_os_columns + N_bob_files*N_bob_columns;

    tot_time_list = cell(Ndata,1);
    tot_data_list = cell(Ndata,1);
    k = 1;

    % Read OpenSim Results
    for i = 1:N_os_columns
        [os_data_list, os_time_list,~] = read_data(opensim_columns(i),opensim_file,"colheaders","");
        
        os_data_list = os_data_list(os_time_list<os_err);
        os_time_list = os_time_list(os_time_list<os_err);
        
        tot_time_list{k} = os_time_list;
        tot_data_list{k} = os_data_list;
        k = k + 1;
    end

    % Read BoB Data
    for i = 1:N_bob_files
        for j = 1:N_bob_columns
            [bob_data_list,bob_time_list,~] = read_data(bob_columns(j),bob_files(i),"textdata",",");
            
            bob_data_list = bob_data_list(bob_time_list < os_err);
            bob_time_list = bob_time_list(bob_time_list < os_err);
            tot_time_list{k} = bob_time_list;
            tot_data_list{k} = bob_data_list;
            k = k + 1;
    
            % Compute correlation coefficient
            os_Ndata = length(os_data_list);
            bob_Ndata = length(bob_data_list);
            Nsample = os_Ndata/bob_Ndata;
            os_data_list_sampled = os_data_list(1:Nsample:end);
            fprintf("Computing Pearson Correlation Coefficient for " + label_txt(i+N_os_columns) + "\n")
            r = pearson_corr_coeff(bob_data_list,os_data_list_sampled);
            fprintf("Computed Correlation Coefficient for " + label_txt(i+N_os_columns) + " is " + r + "\n")
        end
    end
    
    % colors = zeros(Ndata,3);
    % for i = 1:N_os_columns
    %     colors(i,1) = i/N_os_columns;
    % end
    % 
    % for i = N_os_columns+1:Ndata
    %     colors(i,3) = (i-N_os_columns)/(N_bob_files*N_bob_columns);
    % end

    % for i = 1:Ndata
    %     n_samples = length(tot_data_list{i});
    %     tot_data_list{Ndata+i} = mean(tot_data_list{i})*ones(n_samples,1);
    %     tot_time_list{Ndata+i} = tot_time_list{i};
    % end
    
    colors = [[0, 0, 0.73];[0.63,0.63,1];[0.45,0.45,1];[0.73, 0, 0];[1,0.63,0.63];[1, 0.45, 0.45]];
    % colors = [[0,0,0.73];[0.73,0,0];[0,0,0.73];[0.73,0,0]];
    style = ["-","-","-","--","--","--"];

    plot_comparison(tot_time_list,tot_data_list,...
        title_txt,"Force [N]",colors,style,label_txt)
end

% Compare muscle activations
function compare_muscle_activations(bob_files,bob_columns,bob_max_force,opensim_file,opensim_columns,title_txt,label_txt,os_err)
    
    % Determine data size
    N_os_columns = size(opensim_columns,2);
    N_bob_files = size(bob_files,2);
    N_bob_columns = size(bob_columns,2);
    Ndata = N_os_columns + N_bob_files*N_bob_columns;

    tot_time_list = cell(Ndata,1);
    tot_data_list = cell(Ndata,1);
    k = 1;

    % Read OpenSim Results
    for i = 1:N_os_columns
        [os_data_list, os_time_list,~] = read_data(opensim_columns(i),opensim_file,"colheaders","");
        
        os_data_list = os_data_list(os_time_list<os_err);
        os_time_list = os_time_list(os_time_list<os_err);
        
        tot_time_list{k} = os_time_list;
        tot_data_list{k} = os_data_list;
        k = k + 1;
    end

    % Read BoB Data
    for i = 1:N_bob_files
        for j = 1:N_bob_columns
            [bob_data_list,bob_time_list,~] = read_data(bob_columns(j),bob_files(i),"textdata",",");
            bob_data_list = bob_data_list/bob_max_force(j);
            bob_data_list = bob_data_list(bob_time_list < os_err);
            bob_time_list = bob_time_list(bob_time_list < os_err);
            tot_time_list{k} = bob_time_list;
            tot_data_list{k} = bob_data_list;
            k = k + 1;
    
            % Compute correlation coefficient
            os_Ndata = length(os_data_list);
            bob_Ndata = length(bob_data_list);
            Nsample = os_Ndata/bob_Ndata;
            os_data_list_sampled = os_data_list(1:Nsample:end);
            fprintf("Computing Pearson Correlation Coefficient for " + label_txt(i+N_os_columns) + "\n")
            r = pearson_corr_coeff(bob_data_list,os_data_list_sampled);
            fprintf("Computed Correlation Coefficient for " + label_txt(i+N_os_columns) + " is " + r + "\n")
        end
    end
    
    % colors = zeros(Ndata,3);
    % for i = 1:N_os_columns
    %     colors(i,1) = i/N_os_columns;
    % end
    % 
    % for i = N_os_columns+1:Ndata
    %     colors(i,3) = (i-N_os_columns)/(N_bob_files*N_bob_columns);
    % end

    % for i = 1:Ndata
    %     n_samples = length(tot_data_list{i});
    %     tot_data_list{Ndata+i} = mean(tot_data_list{i})*ones(n_samples,1);
    %     tot_time_list{Ndata+i} = tot_time_list{i};
    % end
    
    colors = [[0, 0, 0.73];[0.63,0.63,1];[0.45,0.45,1];[0.73, 0, 0];[1,0.63,0.63];[1, 0.45, 0.45]];
    % style = ["-","-","-","--","--","--"];

    plot_comparison(tot_time_list,tot_data_list,...
        title_txt,"Activation [-]",colors,'',label_txt)

end

% Read BoB - OpenSim data (.csv, .mot, .sto)
function [data_list,time_list,range] = read_data(column_name,file_loc,txt_data,delimiter)
    results = importdata(file_loc);
 
    if(txt_data == "colheaders")
        data_names = results.colheaders;
    elseif(txt_data == "textdata") 
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

% Compute signal magnitude based on data, possibly remove inital offset
function [data_list,time_list] = compute_magnitude(column_names,file_loc,remove_offset)
    
    time_list = read_data("time",file_loc,"colheaders","");
    Ndata = length(time_list);
    data_list = zeros(Ndata,1);
    Ncolumns = length(column_names);
    for i = 1:Ncolumns
        new_data = read_data(column_names(i),file_loc,"colheaders","");

        if(remove_offset)
            new_data = new_data - new_data(1);
        end

        data_list = data_list + new_data.^2;
    end

    data_list = sqrt(data_list);

end

% Plot multiple data series on same plot
function plot_comparison(time_list,data_list,title_txt,ylabel_txt,color,style,legend_txt)

    N = length(data_list);
    if(isempty(color))
        color = zeros(N,3);
        color(:,3) = 0.73*ones(N,1);
    end

    if(isempty(style))
        style = strings(N);
        style(:) = "-";
    end

    figure
    hold on
    grid on
    for i = 1:N
        plot(time_list{i},data_list{i},'LineWidth',1.5,'Color',color(i,:),'LineStyle',style(i))
    end
    
    if(~isempty(legend_txt))
        legend(legend_txt,'Location','best')
    end
    
    title(title_txt)
    ylabel(ylabel_txt)
    xlabel('Time [s]')
    xlim([0 time_list{1}(end)])
    hold off
end

% Compute Pearson Correlation Coefficient
function r = pearson_corr_coeff(data1,data2)
    n = length(data1);

    sumx = sum(data1);
    sumx2 = sum(data1.^2);
    sumy = sum(data2);
    sumy2 = sum(data2.^2);
    sumxy = sum(data1.*data2);

    num = n*sumxy - sumx*sumy;
    den2 = (n*sumx2 - (sumx)^2)*(n*sumy2 - (sumy)^2);

    r = num/sqrt(den2);

end

function rms_value = compute_rms_curve(file_loc, column_name, tlim, textdata, delimiter)

    [data,time,~] = read_data(column_name,file_loc,textdata,delimiter);
    rms_value = rms(data(time<tlim));

end


function rmse_value = compute_rmse(os_file_loc,os_column_name,bob_file_loc,bob_column_name)
    [os_data_list,~,~] = read_data(os_column_name,os_file_loc,"colheaders","");
    [bob_data_list,~,~] = read_data(bob_column_name,bob_file_loc,"colheaders","");

    os_Ndata = length(os_data_list);
    bob_Ndata = length(bob_data_list);
    Nsample = os_Ndata/bob_Ndata;
    os_data_list_sampled = os_data_list(1:Nsample:end);

    rmse_value = rmse(bob_data_list,os_data_list_sampled);
end