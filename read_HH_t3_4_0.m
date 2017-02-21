%% Ashley Lyons, Heriot-Watt University Physics Department, created on 17/07/15, v3.0

% The following code reads the data produced by the Hydraharp 400 using the
% TTTRmain.vi in T3 mode.
% Labview program in the demo libraries provided by PicoQuant.
% Data is produced in a 32 bit binary format where the following bit
% allocation is used in T3 mode:

% 1 - special (used to identify overflows & external markers)
% 2-7 - channel
% 8-22 - dtime (start-stop time between the sync and the event pulses)
% 32-32 - nsync (last logged sync number for the event)

% nsync has a counting limit of 1024 trigger signals. When this limit is
% reached an overflow marker (special = 1 & channel = 63) is sent out to
% signal a reset of the nsync counter.

% To get the global time of the event the number of overflows must be 
% multiplied by 1024 (number of syncs per counter) and added to the nsync 
% value and be combined with the sync repitition rate to get the trigger 
% arrival time and the dtime must be multiplied by the temporal resolution 
% of the Hydraharp to give the start-stop time in non-arbitrary units. 
% The dtime value given is always in units of the resolution (R) allowing 
% for times of up to R*2^15 from the sync event to be recorded.

% For 4 channel measurement the channel IDs are: ch1 = 0, ch2 = 32, ch3 = 16 & ch4 = 48

% Version 3.0 - Coincidences now calculated using temporal window as
% opposed to same sync number. Integration time approximated as time of
% arrival of last sync pulse. Rates now given instead of total Singles &
% coincidence counts. Coincidences no longer counted twice. Bit order now
% flipped for each bit allocation section. Speed improved by reducing the
% number of bi2de & de2bi calls.

% Version 4.0 - Data buffer included for counting coincidences. Rather than
% search through enitre data set for a second event this code now only 
% looks through a number of data points (equal to 2x buffer) surrounding the
% first event. Significant increased in computation time is gained this
% way. Start with large buffer to ensure no loss of coincidences then
% reduce down. Coincidence window has changed by a factor of 2, before a
% window of 10 could only detect events 5ns apart. Window now represents
% the maximum time difference between events allowed.

% Function inputs are data filename and coincidence window in seconds.
% 
function[A, B, C, D, AB, AC, AD, BC, BD, CD, qc_AB, ...
    qc_AC, qc_AD, qc_BC, qc_BD, qc_CD, int_time] = read_HH_t3_4_0(filename, window)

%% Load file
% filename = '44p760000';
fi = fopen(filename);
data = fread(fi,inf,'uint32');
fclose(fi);
savef = 0;       % True - saves processed data to ascii, False - no save

%% Constants
RR = 80e6;      % Laser Rep Rate (Hz)
res = 1e-12;    % Temporal resolution of Hydraharp (s)
% window = 12.5e-9; % Coincidence window (s)

%% Initialise loop values and data matrices
special = zeros(1,length(data));
channel = zeros(6,length(data));
channel_de = zeros(length(data),1);
data_de = zeros(length(data),3);
% counts = zeros(64,1);
ch_bi = fliplr(de2bi(0:63));
ovrflw = 0;
ovrflw_array = zeros(length(data),1);
markers = 0;

ch_63 = [1;1;1;1;1;1];

b = fliplr(de2bi(data,32)); % Convert data to binary for seperation

%% Loop through each 32 bit value, seperate by bit allocation and count single clicks
for ii = 1:length(data)
    
    
        special(1,ii)  = b(ii,1);      % Seperate by bit allocation
        channel(1:6,ii) = b(ii,2:7);
        
        if special(1,ii) == 1 && isequal(channel(1:6,ii),ch_63)   % Count overflows
            ovrflw = ovrflw + 1;
%             counts(64) = counts(64) + 1;
        elseif special(1,ii) == 1 && channel_de(ii) <= 15 && channel_de(ii) ~= 0    % Count external markers (0 unless specified in Labview VI)
            markers = markers + 1;
            
        elseif special(1,ii) == 0       % Count channel clicks & arrange decimal data
            
%             for mm = 1:64
%                 if isequal(channel(:,ii),ch_bi(mm,:)')
%                     counts(mm) = counts(mm) + 1;
%                 end
%             end
            

        end

ovrflw_array(ii) = ovrflw;
end

        data_de(:,1) = bi2de(fliplr(b(:,2:7)));     % Convert data back to decimal for processing
        data_de(:,3) = bi2de(fliplr(b(:,23:32))) + 1024*ovrflw_array;
        data_de(:,2) = res*bi2de(fliplr(b(:,8:22))) + (1/RR)*data_de(:,3);
        
        D_f = data_de(:,1) ~= 63;
        data_de2 = data_de(D_f,:);

int_time = max(data_de(:,3))*(RR)^(-1);     % Integration time calculated from the time of the last sync event (s)

%% Identify unique channel values and form a list of their counts
channel_de2 = unique(data_de(:,1));
channels = length(unique(channel_de2));
cnts = zeros(channels,2);

for x = 1:channels
    sings = length(find(data_de(:,1) == channel_de2(x)));
        cnts(x,1) = channel_de2(x);
        cnts(x,2) = sings;
end 

%% Coincidences
% Initialise single and multi-channel counts
ones = 0;
twos = 0;
threes = 0;
fours = 0;
times_ones = [];    % Timetags for which there is a single count
times_twos = [];      % Timetags for which there is a 2-way coincidence
times_threes = [];  % Timetags for which there is a 3-way coincidence
times_fours = [];   % Timetags for which there is a 4-way coincidence
chs_2 = [];     % Pairs of channel values that give a 2-way coincidence
i1 = 1;
i2 = 1;
i3 = 1;
i4 = 1;

timetags = data_de2(:,2);      % Removes data before first sync, external markers and overflows

% For each event find other events within the given window. This is done by
% searching through a buffer of events local to the first event.

buffer = 50;   % half the size of the number of local events to search through

len_counts = length(data_de2(:,2));     % Number of events

for jj = 1:length(timetags)
   time = timetags(jj);
   buff_start = jj-buffer+1;            % Index of start of buffer - Minimum of 1
   if buff_start <= 0
       buff_start = 1;
   end
   buff_end = jj+buffer-1;              % Index of end of buffer - Maximum of index of final event
   if buff_end-len_counts > 0
       buff_end = length(data_de2(:,2));
   end
   buff_array = (buff_start):(buff_end);        % Create array of indices to call the timetags for
   times2 = zeros(length(buff_array),2);
   times2(:,1) = buff_array;                % Stores in-buffer event locations 
   times2(:,2) = data_de2(buff_array,2) - time;     % Time difference between first event and in-buffer events
   times3 = abs(times2(:,2)) <= window;           % Checks against coincidence window
   z = find(times3);
   time_loc = times2(z,1);
    if length(z) == 1
        ones = ones + 1;
        times_ones(i1,:) = timetags(time_loc);
        i1 = i1+1;
    elseif length(z) == 2 && not(ismember(data_de2(time_loc(1,1),2),times_twos))     % ismember prevents double counting
        twos = twos + 1;
        times_twos(i2,:) = timetags(time_loc);
        chs_2(i2,:) = data_de2(time_loc,1);
        i2 = i2+1;
    elseif length(z) == 3 && not(ismember(data_de2(time_loc(1,1),2),times_threes))
        threes = threes + 1;
        times_threes(i3,:) = timetags(time_loc);
        i3 = i3+1;
    elseif length(z) == 4 && not(ismember(data_de2(time_loc(1,1),2),times_fours))
        times_fours(i4,:) = timetags(time_loc);
        fours = fours + 1;
        i4 = i4+1;
    end
end

%% Coincidence counting
% Initialisation

A_v = channel_de2(1);   % Values allocated to the channels (see preamble)
B_v = channel_de2(2);
C_v = channel_de2(3);
D_v = channel_de2(4);

AB = 0;     % Initialiased coincidence counts
AC = 0;
AD = 0;
BC = 0;
BD = 0;
CD = 0;

%% Coincidence conditionals

if ~isempty(chs_2)
for kk = 1:length(chs_2(:,1))
    p = chs_2(kk,1);
    q = chs_2(kk,2);
        if p == A_v
            if q == A_v
                continue
            elseif q == B_v
                AB = AB + 1;
            elseif q == C_v
                AC = AC + 1;
            elseif q == D_v
                AD = AD + 1;
            end
        elseif p == B_v
            if q == A_v
                AB = AB + 1;
            elseif q == B_v
                continue
            elseif q == C_v
                BC = BC + 1;
            elseif q == D_v
                BD = BD + 1;
            end
        elseif p == C_v
            if q == A_v
                AC = AC + 1;
            elseif q == B_v
                BC = BC + 1;
            elseif q == C_v
                continue
            elseif q == D_v
                CD = CD + 1;
            end
        elseif p == D_v
            if q == A_v
                AD = AD + 1;
            elseif q == B_v
                BD = BD + 1;
            elseif q == C_v;
                CD = CD + 1;
            elseif q == D_v
                continue
            end
        end
end
end

coin_tot = AB + AC + AD + BC + BD + CD;

% Singles counts values
A = cnts(cnts(:,1) == A_v,2);
B = cnts(cnts(:,1) == B_v,2);
C = cnts(cnts(:,1) == C_v,2);
D = cnts(cnts(:,1) == D_v,2);

% Count Rates
A_r = cnts(cnts(:,1) == A_v,2)/int_time;
B_r = cnts(cnts(:,1) == B_v,2)/int_time;
C_r = cnts(cnts(:,1) == C_v,2)/int_time;
D_r = cnts(cnts(:,1) == D_v,2)/int_time;

AB_r = AB/int_time;
AC_r = AC/int_time;
AD_r = AD/int_time;
BC_r = BC/int_time;
BD_r = BD/int_time;
CD_r = CD/int_time;

% Accidentals & Quantum Contrast
acc_AB = A_r*B_r*window;        % Accidentals
acc_AC = A_r*C_r*window;
acc_AD = A_r*D_r*window;
acc_BC = B_r*C_r*window;
acc_BD = B_r*D_r*window;
acc_CD = C_r*D_r*window;

qc_AB = AB_r/acc_AB;        % Quantum Contrast
qc_AC = AC_r/acc_AC;
qc_AD = AD_r/acc_AD;
qc_BC = BC_r/acc_BC;
qc_BD = BD_r/acc_BD;
qc_CD = CD_r/acc_CD;

sync_count_ovrflw = ovrflw;
external_markers = markers;

if savef == 1
    save_file = strrep(filename,'.out','_data');
    save(save_file,'data_de','-ascii');
end