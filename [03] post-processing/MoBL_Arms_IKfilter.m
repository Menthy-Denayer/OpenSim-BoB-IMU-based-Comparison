%filter kinematic data
%edited by MEV 8/14/14
%edited by Menthy 19/01/2024
clear all;
close all;
clc;

%% Load Data
OPENSIM_RESULTS_FOLDER = "C:\Users\medenaye\Documents\2324\5. AI4Exo\3. OpenSim\2. results\BIOROB\";
filename = OPENSIM_RESULTS_FOLDER + "IK_imu_test.mot";

out_filename = OPENSIM_RESULTS_FOLDER + "filteredIK_test.sto";
reset_6dof = false;

if(reset_6dof)
    ncolumns=21;
else
    ncolumns=27;
end

delimiterIn = '\t';headerlinesIn = 11;
rawIK = importdata(filename,delimiterIn,headerlinesIn);
timevector = rawIK.data(:,1);

%% Filter design
windowSize = 30;
Hd=ones(1,windowSize)/windowSize;

%% Filter data 

if(reset_6dof)
    n_start = 8;
else
    n_start = 2;
end

for n=n_start:27
    filtered_IK(:,n) = filtfilt(Hd,1,rawIK.data(:,n));
end

%% Shift translation vector
% filtered_IK(:,6) = filtered_IK(:,6) - filtered_IK(1,6);

%% Filtered data with time into matrix
data_together = [timevector,filtered_IK(:,n_start:27)];

%% Create column headers without r_x,r_y,r_z,t_x,t_y,t_z, MARKERS ONLY
% new_textdata_headers=[rawIK.textdata(1,1);rawIK.textdata(2,1);rawIK.textdata(3,1);"nColumns="+ncolumns;rawIK.textdata(5,1);rawIK.textdata(8,1)];
% time_header=rawIK.textdata(9,1);
% dof_headers=rawIK.textdata(9,n_start:27);
% new_column_headers=[time_header,dof_headers];

%% Create column headers without r_x,r_y,r_z,t_x,t_y,t_z, IMU ONLY
new_textdata_headers=[rawIK.textdata(1,1);rawIK.textdata(2,1);rawIK.textdata(3,1);"nColumns="+ncolumns;rawIK.textdata(4,1);rawIK.textdata(6,1)];
headers=strsplit(string(rawIK.textdata(7,1)),"	");
time_header = headers(1);
dof_headers=headers(n_start:27);
new_column_headers=[time_header,dof_headers];

%% Export filtered data to .mot file

fileID=fopen(out_filename,'wt');
[M,N]=size(new_textdata_headers);
for i=1:M;
    for j=1:N;
fprintf(fileID,'%s\n',new_textdata_headers{i,j});
    end
end
[P,Q]=size(new_column_headers);
for k=1:P;
    for m=1:Q;
        fprintf(fileID,'%s\t',new_column_headers{k,m});
    end
end
fprintf(fileID,'\n');
fclose(fileID);
dlmwrite(out_filename,...
         data_together,'-append','delimiter','\t','newline','pc');
