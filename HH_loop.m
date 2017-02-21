%% Ashley Lyons, Heriot-Watt University Physics Department, created on 17/07/15, v1.0
% The following code calls the read_HH_t3 function in a parallel loop to analyse
% data collected from a scan of the piezo translation stage. File names are 
% expected to be the position of the translation stage in microns. See
% read_HH_t3_3_0.m for more details on the function.

clear
close all
clc

folder = 'C:\Users\Trial\Documents\MATLAB\LOOP_Graphene\30000'; % Directory of data files
folder_home = 'C:\Users\Trial\Documents\HOM_data';  % Directory of this file
cd(folder_home)
files = dir(folder);
names1 = {files.name}'; % Gather list of data file names

c1 = 1;
for bb = 1:length(names1)
    n1 = names1(bb);
    n2 = strrep(n1,'p','.');
    n4 = char(n2(1));
    n3 = str2double(n2);
    if length(n4) > 3 && strcmp(n4(end-2:end), 'out') == 1
        names(c1) = names1(bb);
        c1 = c1+1;
    else
        continue
    end
end
data_full = zeros(length(names),18); % Initialise array of correct size
% h = waitbar(0,'Calculating...');
window = 10e-9; % Coincidence window size (seconds)

cd(folder)

%% Loop through each data file and use analysis function. X_r gives singles 
% rates, XY_r gives coincidence rates, qc_XY gives Quantum Contrast (CAR),
% int_time gives integration time of measurement.

parfor aa = 1:length(names)
    % waitbar(aa/length(names),h)
    filename = char(names(aa));
    
        [A_r, B_r, C_r, D_r, AB_r, AC_r, AD_r, BC_r, BD_r, CD_r, qc_AB, ...
    qc_AC, qc_AD, qc_BC, qc_BD, qc_CD, int_time] = read_HH_t3_4_0(filename, window);
        pos1 = strrep(filename,'.out','');
        pos2 = strrep(pos1,'p','.');
        position = str2double(pos2); % Filename converted to position value
    
        data_full(aa,:) = [position, A_r, B_r, C_r, D_r, AB_r, AC_r, AD_r, BC_r, BD_r, CD_r, qc_AB, ...
    qc_AC, qc_AD, qc_BC, qc_BD, qc_CD, int_time];

end

%% Runs a python code to email user when completed, saves data and plays a sound file

cd(folder_home)
% gm = 'python C:\Users\Trial\Documents\MATLAB\gmail_codefinish.pyw HH_loop';
% gmstr = [gm,' ',char(39), mfilename,'_test_dip',char(39)];
% dos(gmstr);

save('C:\Users\Trial\Documents\MATLAB\LOOP\processed_19-08-16\30000_graphene\1','data_full')
[y,Fs] = audioread('Windows Notify.wav');
sound(y,Fs);