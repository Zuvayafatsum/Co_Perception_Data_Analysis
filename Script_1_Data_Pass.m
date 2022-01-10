%Load ALLEEG (final), comments and behavioral (info) files to the workspace.

%% do it for the first one only !!
Subject_ID = 'CPS05V2';
data = {};
data.EEG = {};
data.coherence = [];
data.RT_Behavior = [];
data.Confidence = [];
data.answers = []; 
data.RT_EEG = [];
data.exp_condition = [];
data.visual_field = [];
data.subjects_ID = Subject_ID ;
data.time = ALLEEG.times;
data.block_num = [];
for i = 1:ALLEEG.nbchan;
data.Channel_Labels{i,1} = ALLEEG.chanlocs(1,i).labels;
end
%% Load data files (per block)
dir = 'D:\Program Files\MATLAB\eeglab2019_0';
Subject_ID = 'CPS05V2';
Block_Num =1;

%%Eeglab comment file
load([dir '\CP_ALLCom\' Subject_ID '_' num2str(Block_Num) '_comment.mat'])

%%Eeglab data file
load([dir '\CP_ALLEEG_Final\' Subject_ID '_' num2str(Block_Num) '_Final.mat'])

%%Behavioral data
load([dir '\CP_Behav\' Subject_ID '_' num2str(Block_Num) '.mat'])

%% Do these steps for each block seperately!!!

length_data = length(data.EEG); 

for i = length_data+1 : length_data+ length(ALLEEG.epoch)
data.EEG{i,1} = double(ALLEEG.data(:,:,i - length_data));
end

%% Save block number per trial
length_block_num = length(data.block_num);
for i = length_block_num+1 : length_block_num+ length(ALLEEG.epoch)
data.block_num(i,1) =Block_Num;
end
%% DÝKKAT
%!!!!!Find number of rejected trials from ALLCOM variable. jimmy script !!!!!!
cell_idx = [];
start_idx  = [];
end_idx  = [];
tmpidx  = [];
invalid_trials = [];

for i = 1:length(ALLCOM)
if strfind(ALLCOM{i},'pop_rejepoch') > 1;
cell_idx(i) = i;
else
cell_idx(i) = 0;
end
end
cell_idx = cell_idx(cell_idx~=0);

start_idx = strfind(ALLCOM{cell_idx}, '[');
end_idx = strfind(ALLCOM{cell_idx}, ']');

tmpidx = ALLCOM{cell_idx}(start_idx+1: end_idx-1);

invalid_trials = str2num(tmpidx);

%Binary coding of participants' reponses
answers = [];
answers = double(info.eval_answer);
for i = 1:length(answers);
    if answers(i) == 116;
        answers(i) = 2;
    else
        answers(i) = 1;
    end
end
answers(invalid_trials) = [];


%Pass binary response data to the data structure

length_answers = length(data.answers);

for i = length_answers + 1 : length_answers + length(answers);
    data.answers(i,1) = answers(i - length_answers);
end

%Take coherence, confidence and RT(behavioral) information and eleminate
%rejected trials
coh = [];
RT_Behavior = [];
conf = [];


coh = info.coh;
RT_Behavior = info.RT;
conf = info.Conf;

coh(invalid_trials) = [];
RT_Behavior(invalid_trials) = [];
conf(invalid_trials) = [];

%Pass this info to the data structure that we create


length_coherence = length(data.coherence);

for i = length_coherence+1 : length_coherence+ length(coh);
data.coherence(i,1) = double(coh(i-length_coherence));
data.RT_Behavior(i,1) =  double(RT_Behavior(i-length_coherence));
data.Confidence(i,1)=  double(abs(conf(i-length_coherence)));
end

%% Create a new vector containing info for publich/private trials
%!!! 1 = public trials and 2 = private trials!!!
length_exp_condition = length(data.exp_condition);
exp_condition = [];

if info.observerlocation == "R";
    exp_condition(info.RightSide,1) = 1; % 1 for public trails
    exp_condition(info.LeftSide,1) = 2; % 2 for private trials
else
    exp_condition(info.LeftSide,1) = 1; % 1 for public trials
    exp_condition(info.RightSide,1) = 2; %2 for private trials
end



exp_condition(invalid_trials) = [];

for i = length_exp_condition+1 : length_exp_condition+ length(ALLEEG.epoch);
    data.exp_condition(i,1) = exp_condition(i-length_exp_condition);
end

%% create vector containing which visual field the stimulus was presented 
% 1= Left, 2 = Right)
visual_field = [];
length_visual_field = length(data.visual_field);

visual_field(info.LeftSide) = 1;
visual_field(info.RightSide) = 2;

visual_field(invalid_trials) = [];

for i = length_visual_field+1 : length_visual_field+ length(ALLEEG.epoch);
    data.visual_field(i,1) = visual_field(i-length_visual_field);
end




%Take EEG data RT info from EEGLAB file and pass it to the new data
%structure

length_RT_EEG = length(data.RT_EEG)

for i = length_RT_EEG+1:length_RT_EEG + length(ALLEEG.epoch)
    if sum(size(ALLEEG.epoch(i-length_RT_EEG).eventlatency)) > 2
data.RT_EEG(i,1) = ALLEEG.epoch(i-length_RT_EEG).eventlatency{1, 2}
    else
        data.RT_EEG(i,1) = 5000;
end
end




clearvars -except data




