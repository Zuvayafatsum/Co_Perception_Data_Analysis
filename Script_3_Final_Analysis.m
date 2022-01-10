%% Load necessary files
load chanlocs.mat ;
load Corrdinates64.mat; 
%load subjects data (e.g., S01_combined.mat)
%% trial type vector
data.trial_type(find(data.coherence == 3.2 & data.exp_condition == 1)) = 1; % social and 3.2 coherence
data.trial_type(find(data.coherence == 6.4 & data.exp_condition == 1)) = 2; % social and 6.4 coherence
data.trial_type(find(data.coherence == 12.8 & data.exp_condition == 1)) = 3; % social and 12.8 coherence
data.trial_type(find(data.coherence == 51.2 & data.exp_condition == 1)) = 4; % social and 51.2 coherence

data.trial_type(find(data.coherence == 3.2 & data.exp_condition == 2)) = 5; % private and 3.2 coherence
data.trial_type(find(data.coherence == 6.4 & data.exp_condition == 2)) = 6; % private and 6.4 coherence
data.trial_type(find(data.coherence == 12.8 & data.exp_condition == 2)) = 7; % private and 12.8 coherence
data.trial_type(find(data.coherence == 51.2 & data.exp_condition == 2)) = 8; % private and 51.2 coherence

% channel indicies based on 64 electrode cap we use. (see image on desktop
% "acticap_image" 
left_chans = [11 12 14 43 44 45 47];
right_chans = [19 22 23 49 51 52 54];
all_cp_chans = [11 12 14 43 44 45 47 19 22 23 49 51 52 54 13 48 53];

data.left_chans = left_chans;
data.right_chans = right_chans;
data.all_cp_chans = all_cp_chans;

% fixate color scale per coherence level
color_map = [0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980; 0.6350 0.0780 0.1840; 0.4940 0.1840 0.5560];
data.color_map = color_map;

% combine trial information for visual field of stimulus presentation
all_trial_inf(:,1) = data.trial_type; % 1 = coh1 social, 5 =  coh1 priv
all_trial_inf(:,2) = data.visual_field; % 1 == left 2 == right
all_trial_inf(:,3) = data.coherence;
data.all_trial_inf = all_trial_inf;

data.Fsample = 500;
clearvars -except data X Y Z chanlocs group_data

%% Reorder EEG matrix and apply Laplacian filtering

data.EEG = cell2mat(arrayfun(@(x)permute(x{:},[3 2 1]),data.EEG,'UniformOutput',false));
data.EEG = permute(data.EEG,[3 2 1]);

Laplacian_filtered = zeros(size(data.EEG));
for i=1:size(data.EEG,3)
    tmpdata = data.EEG(:,:,i);
    Laplacian_filtered(:,:,i) = laplacian_perrinX(tmpdata,X,Y,Z);
end

data.Laplacian_filtered = Laplacian_filtered;
 
clearvars -except data X Y Z chanlocs group_data

%% Baseline correction ( all trials regardless of conditions etc.)
baseline_start = -200; %in miliseconds
baseline_end = 0; %in miliseconds
bc_start_indx = find(data.time == baseline_start);
bc_end_indx = find(data.time == baseline_end);
%calculate mean baseline activity per chan and per trial

bc_per_chan_trial = zeros(size(data.EEG,1),size(data.EEG,3));

for i = 1:size(data.EEG,1);
    for j = 1:size(data.EEG,3);
        bc_per_chan_trial(i,j) = mean(data.Laplacian_filtered(i,bc_start_indx:bc_end_indx,j),2);
    end
end
bc_per_chan_trial = double(bc_per_chan_trial);

%substract baseline activity for each chan and trial from epoch
data.bc_corrected = zeros(size(data.EEG));
for i = 1:size(data.EEG,1); 
    for j = 1:size(data.EEG,3);
        data.bc_corrected(i,:,j) = data.Laplacian_filtered(i,:,j) - bc_per_chan_trial(i,j);
    end
end

clearvars -except data X Y Z chanlocs group_data

% sanity check for baseline correction
start_plot = -200;
end_plot = 1300;
start_indx = find(data.time == start_plot);
end_indx = find(data.time == end_plot);

set(0,'DefaultLineLineWidth',1.5)
plot(data.time(start_indx:end_indx),(mean(mean(data.Laplacian_filtered(data.all_cp_chans,start_indx:end_indx,:),3),1)),'color','b')
hold on
plot(data.time(start_indx:end_indx),(mean(mean(data.bc_corrected(data.all_cp_chans,start_indx:end_indx,:),3),1)))
hold on
yline([0],'--k')
xline([0],'--k')
legend('Laplacian','bc corrected')

%% Calculate mean CPP per trial (bilaterally over all centro parietal region)
data.overall_CPP = zeros(size(data.bc_corrected,3),size(data.bc_corrected,2));

for i = 1:size(data.Laplacian_filtered,3);
    data.overall_CPP(i,:) = mean(data.bc_corrected(data.all_cp_chans,:,i),1);
end

% sanity check for CPP
start_plot = -200;
end_plot = 1500;
start_indx = find(data.time == start_plot);
end_indx = find(data.time == end_plot);

set(0,'DefaultLineLineWidth',1.5)
plot(data.time(start_indx:end_indx),mean(data.overall_CPP(:,start_indx:end_indx),1))
hold on
yline([0],'--k')
xline([0],'--k')
legend('CPP All conditions')
ylim([-max(mean(data.overall_CPP(:,start_indx:end_indx),1)) max(mean(data.overall_CPP(:,start_indx:end_indx),1))])
%% Calculate CPP moving average (based on CPP signal we just created)
movmean_start = -200;
start_indx = find(data.time == movmean_start);

window_size = 25;

data.CPP_movmean = zeros(size(data.overall_CPP));
data.CPP_movmean_plot  = zeros(size(data.overall_CPP));

for i = 1:size(data.bc_corrected,3)
data.CPP_movmean(i,start_indx:end) = movmean(data.overall_CPP(i,start_indx:end),[0 window_size]);
data.CPP_movmean_plot(i,start_indx:end) = data.CPP_movmean(i,start_indx:end) - mean(data.CPP_movmean(i,start_indx:751),2);
end

%check moving average
start_plot = -200;
end_plot =1000;
start_indx = find(data.time == start_plot);
end_indx = find(data.time == end_plot);

figure(1)
plot(data.time(start_indx:end_indx),mean(data.CPP_movmean(:,start_indx:end_indx),1))
hold on
yline([0],'--k')
xline([0],'--k')
legend('CPP All conditions')
ylim([-max(mean(data.CPP_movmean(:,start_indx:end_indx),1)) max(mean(data.CPP_movmean(:,start_indx:end_indx),1))])

figure(2)
set(0,'DefaultLineLineWidth',1.5)
plot(data.time(start_indx:end_indx),mean(data.CPP_movmean(data.trial_type == 1,start_indx:end_indx),1),'color',data.color_map(1,:))
hold on
plot(data.time(start_indx:end_indx),mean(data.CPP_movmean(data.trial_type == 2,start_indx:end_indx),1),'color',data.color_map(2,:))
hold on
plot(data.time(start_indx:end_indx),mean(data.CPP_movmean(data.trial_type == 3,start_indx:end_indx),1),'color',data.color_map(3,:))
hold on
plot(data.time(start_indx:end_indx),mean(data.CPP_movmean(data.trial_type == 4,start_indx:end_indx),1),'color',data.color_map(4,:))
hold on
plot(data.time(start_indx:end_indx),mean(data.CPP_movmean(data.trial_type == 5,start_indx:end_indx),1),'--','color',data.color_map(1,:))
hold on
plot(data.time(start_indx:end_indx),mean(data.CPP_movmean(data.trial_type == 6,start_indx:end_indx),1),'--','color',data.color_map(2,:))
hold on
plot(data.time(start_indx:end_indx),mean(data.CPP_movmean(data.trial_type == 7,start_indx:end_indx),1),'--','color',data.color_map(3,:))
hold on
plot(data.time(start_indx:end_indx),mean(data.CPP_movmean(data.trial_type == 8,start_indx:end_indx),1),'--','color',data.color_map(4,:))
hold on
legend('S1','S2','S3','S4','P1','P2','P3','P4')
legend('AutoUpdate','off')
sgtitle('Social vs Private Evidence Accumulation Per Coherence')
yline([0],'--k')
xline([0],'--k')
ylabel('Magnitude (mv/s2)')
xlabel('Time (ms)')



data.movemean_window_size = window_size;

clearvars -except data X Y Z chanlocs group_data
%% Calculate diff waves per coherence level & Visual field of stimulus presentation
% find trial indexes satisfying specified condition
leftvf_coh1_social = find(data.all_trial_inf(:,1) == 1 & data.all_trial_inf(:,2) == 1); % Social & 3.2 & Left
leftvf_coh2_social = find(data.all_trial_inf(:,1) == 2 & data.all_trial_inf(:,2) == 1); % Social & 6.4 & Left
leftvf_coh3_social = find(data.all_trial_inf(:,1) == 3 & data.all_trial_inf(:,2) == 1); % Social & 12.8 & Left
leftvf_coh4_social = find(data.all_trial_inf(:,1) == 4 & data.all_trial_inf(:,2) == 1); % Social & 51.2 & Left

rightvf_coh1_social = find(data.all_trial_inf(:,1) == 1 & data.all_trial_inf(:,2) == 2); % Social & 3.2 & Right
rightvf_coh2_social = find(data.all_trial_inf(:,1) == 2 & data.all_trial_inf(:,2) == 2); % Social & 6.4 & Right
rightvf_coh3_social = find(data.all_trial_inf(:,1) == 3 & data.all_trial_inf(:,2) == 2); % Social & 12.8 & Right
rightvf_coh4_social = find(data.all_trial_inf(:,1) == 4 & data.all_trial_inf(:,2) == 2); % Social & 51.2 & Right

leftvf_coh1_priv = find(data.all_trial_inf(:,1) == 5 & data.all_trial_inf(:,2) == 1); % Priv & 3.2 & Left
leftvf_coh2_priv = find(data.all_trial_inf(:,1) == 6 & data.all_trial_inf(:,2) == 1); % Priv & 6.4 & Left
leftvf_coh3_priv = find(data.all_trial_inf(:,1) == 7 & data.all_trial_inf(:,2) == 1); % Priv & 12.8 & Left
leftvf_coh4_priv = find(data.all_trial_inf(:,1) == 8 & data.all_trial_inf(:,2) == 1); % Priv & 51.2 & Left

rightvf_coh1_priv = find(data.all_trial_inf(:,1) == 5 & data.all_trial_inf(:,2) == 2); % Priv & 3.2 & Right
rightvf_coh2_priv = find(data.all_trial_inf(:,1) == 6 & data.all_trial_inf(:,2) == 2); % Priv & 3.2 & Right
rightvf_coh3_priv = find(data.all_trial_inf(:,1) == 7 & data.all_trial_inf(:,2) == 2); % Priv & 3.2 & Right
rightvf_coh4_priv = find(data.all_trial_inf(:,1) == 8 & data.all_trial_inf(:,2) == 2); % Priv & 3.2 & Right

% mean for each coherence level in SOCIAL condition for LEFT visual field
% presentations
coh1_leftv_soc_mean = mean(mean(data.bc_corrected(data.right_chans,:,leftvf_coh1_social),3),1);
coh2_leftv_soc_mean = mean(mean(data.bc_corrected(data.right_chans,:,leftvf_coh2_social),3),1);
coh3_leftv_soc_mean = mean(mean(data.bc_corrected(data.right_chans,:,leftvf_coh3_social),3),1);
coh4_leftv_soc_mean = mean(mean(data.bc_corrected(data.right_chans,:,leftvf_coh4_social),3),1);

% mean for each coherence level in SOCIAL condition for RIGHT visual field
% presentations
coh1_rightv_soc_mean = mean(mean(data.bc_corrected(data.left_chans,:,rightvf_coh1_social),3),1);
coh2_rightv_soc_mean = mean(mean(data.bc_corrected(data.left_chans,:,rightvf_coh2_social),3),1);
coh3_rightv_soc_mean = mean(mean(data.bc_corrected(data.left_chans,:,rightvf_coh3_social),3),1);
coh4_rightv_soc_mean = mean(mean(data.bc_corrected(data.left_chans,:,rightvf_coh4_social),3),1);

% mean for each coherence level in PRIVATE condition for LEFT visual field
% presentations
coh1_leftv_priv_mean = mean(mean(data.bc_corrected(data.right_chans,:,leftvf_coh1_priv),3),1);
coh2_leftv_priv_mean = mean(mean(data.bc_corrected(data.right_chans,:,leftvf_coh2_priv),3),1);
coh3_leftv_priv_mean = mean(mean(data.bc_corrected(data.right_chans,:,leftvf_coh3_priv),3),1);
coh4_leftv_priv_mean = mean(mean(data.bc_corrected(data.right_chans,:,leftvf_coh4_priv),3),1);

% mean for each coherence level in PRIVATE condition for RIGHT visual field
% presentations
coh1_rightv_priv_mean = mean(mean(data.bc_corrected(data.left_chans,:,rightvf_coh1_priv),3),1);
coh2_rightv_priv_mean = mean(mean(data.bc_corrected(data.left_chans,:,rightvf_coh2_priv),3),1);
coh3_rightv_priv_mean = mean(mean(data.bc_corrected(data.left_chans,:,rightvf_coh3_priv),3),1);
coh4_rightv_priv_mean = mean(mean(data.bc_corrected(data.left_chans,:,rightvf_coh4_priv),3),1);

%Calculate difference waves per coherence level and per visual field, e.g.
%calculate social - private for 3.2 coherence and left visual field
%presentations and do the same for right visual field presentations. Then
%sum these two difference waves per each coherence level to have total diff
%wave for 3.2 (SOCIAL - PRIV)

coh1_diff_left_visual = coh1_leftv_soc_mean - coh1_leftv_priv_mean;
coh2_diff_left_visual = coh2_leftv_soc_mean - coh2_leftv_priv_mean;
coh3_diff_left_visual = coh3_leftv_soc_mean - coh3_leftv_priv_mean;
coh4_diff_left_visual = coh4_leftv_soc_mean - coh4_leftv_priv_mean;

coh1_diff_right_visual = coh1_rightv_soc_mean - coh1_rightv_priv_mean;
coh2_diff_right_visual = coh2_rightv_soc_mean - coh2_rightv_priv_mean;
coh3_diff_right_visual = coh3_rightv_soc_mean - coh3_rightv_priv_mean;
coh4_diff_right_visual = coh4_rightv_soc_mean - coh4_rightv_priv_mean;

%Calculate total diff per coherenc level
coh1_total_diff = coh1_diff_left_visual + coh1_diff_right_visual;
coh2_total_diff = coh2_diff_left_visual + coh2_diff_right_visual;
coh3_total_diff = coh3_diff_left_visual + coh3_diff_right_visual;
coh4_total_diff = coh4_diff_left_visual + coh4_diff_right_visual;

data.coh1_total_diff = coh1_total_diff;
data.coh2_total_diff = coh2_total_diff;
data.coh3_total_diff = coh3_total_diff;
data.coh4_total_diff = coh4_total_diff;


%% Plot difference waves per coherence levels
start_plot = -200;
end_plot = 1300;
start_indx = find(data.time == start_plot);
end_indx = find(data.time == end_plot);
a = data.coh1_total_diff(start_indx:end_indx);
b = data.coh2_total_diff(start_indx:end_indx);
c = data.coh3_total_diff(start_indx:end_indx);
d = data.coh4_total_diff(start_indx:end_indx);

plot(data.time(start_indx:end_indx),a,'color',data.color_map(1,:))
hold on
plot(data.time(start_indx:end_indx),b,'color',data.color_map(2,:))
hold on
plot(data.time(start_indx:end_indx),c,'color',data.color_map(3,:))
hold on
plot(data.time(start_indx:end_indx),d,'color',data.color_map(4,:))
hold on
legend('Coh1 diff','Coh2 diff','Coh3 diff','Coh4 diff')
legend('AutoUpdate','off')
hold on
xline([0],'--k')
yline([0],'--k')
sgtitle("'(Social - Private)' over Centro-parietal channels")
ylabel('Magnitude (mv/s2)')
xlabel('Time (ms)')
clearvars -except data X Y Z chanlocs group_data
%% Calculate difference waves Social - Private regardless of coherence level
% calculate averages based on overall CPP data
overall_social_mean = mean(data.overall_CPP(find(data.exp_condition == 1),:),1);
overall_priv_mean = mean(data.overall_CPP(find(data.exp_condition == 2),:),1);
% calculate difference
overall_diff = overall_social_mean - overall_priv_mean;

data.overall_diff = overall_diff;
data.overall_social_mean = overall_social_mean;
data.overall_priv_mean = overall_priv_mean;


%plot data
start_plot = -200;
end_plot = 1300;
start_indx = find(data.time == start_plot);
end_indx = find(data.time == end_plot);

figure(1)
plot(data.time(start_indx:end_indx),data.overall_social_mean(start_indx:end_indx))
hold on
plot(data.time(start_indx:end_indx),data.overall_priv_mean (start_indx:end_indx))
hold on
plot(data.time(start_indx:end_indx),data.overall_diff(start_indx:end_indx),'color',data.color_map(1,:))
legend('Social','Private','Social - Private')
legend('AutoUpdate','off')
hold on
xline([0],'--k')
yline([0],'--k')
ylabel('Magnitude (mV)')
xlabel('Time (ms)')
sgtitle("Overall Social and Private over Centro-Parietal Channels (All Coherences)")

figure(2)
plot(data.time(start_indx:end_indx),mean(data.CPP_movmean_plot(data.exp_condition == 1,start_indx:end_indx),1))
hold on
plot(data.time(start_indx:end_indx),mean(data.CPP_movmean_plot(data.exp_condition == 2,start_indx:end_indx),1))
legend('Social','Private')
legend('AutoUpdate','off')
xline([0],'--k')
yline([0],'--k')
sgtitle('Overall Social vs Private Evidence Accumulation')
ylabel('Magnitude (mv/s2)')
xlabel('Time (ms)')

clearvars -except data X Y Z chanlocs group_data

%% Calculate overall (all channels) CPP ERPs based on coherence level and exp condition
coh1_social_CPP = mean(data.overall_CPP(data.coh1_social,:),1);
coh2_social_CPP = mean(data.overall_CPP(data.coh2_social,:),1);
coh3_social_CPP = mean(data.overall_CPP(data.coh3_social,:),1);
coh4_social_CPP = mean(data.overall_CPP(data.coh4_social,:),1);

coh1_priv_CPP = mean(data.overall_CPP(data.coh1_private,:),1);
coh2_priv_CPP = mean(data.overall_CPP(data.coh2_private,:),1);
coh3_priv_CPP = mean(data.overall_CPP(data.coh3_private,:),1);
coh4_priv_CPP = mean(data.overall_CPP(data.coh4_private,:),1);

data.coh1_social_CPP = coh1_social_CPP;
data.coh2_social_CPP = coh2_social_CPP;
data.coh3_social_CPP = coh3_social_CPP;
data.coh4_social_CPP = coh4_social_CPP;

data.coh1_priv_CPP = coh1_priv_CPP;
data.coh2_priv_CPP = coh2_priv_CPP;
data.coh3_priv_CPP = coh3_priv_CPP;
data.coh4_priv_CPP = coh4_priv_CPP;

% Plot signals
start_plot = -200;
end_plot = 1200;
start_indx = find(data.time == start_plot);
end_indx = find(data.time == end_plot);

plot(data.time(start_indx:end_indx),data.coh1_social_CPP(start_indx:end_indx),'color',data.color_map(1,:))
hold on
plot(data.time(start_indx:end_indx),data.coh2_social_CPP(start_indx:end_indx),'color',data.color_map(2,:))
hold on
plot(data.time(start_indx:end_indx),data.coh3_social_CPP(start_indx:end_indx),'color',data.color_map(3,:))
hold on
plot(data.time(start_indx:end_indx),data.coh4_social_CPP(start_indx:end_indx),'color',data.color_map(4,:))
hold on
plot(data.time(start_indx:end_indx),data.coh1_priv_CPP(start_indx:end_indx),':','color',data.color_map(1,:))
hold on
plot(data.time(start_indx:end_indx),data.coh2_priv_CPP(start_indx:end_indx),':','color',data.color_map(2,:))
hold on
plot(data.time(start_indx:end_indx),data.coh3_priv_CPP(start_indx:end_indx),':','color',data.color_map(3,:))
hold on
plot(data.time(start_indx:end_indx),data.coh4_priv_CPP(start_indx:end_indx),':','color',data.color_map(4,:))

sgtitle("'(Social - Private)' over Centro-parietal channels")
legend('Coh1 Social','Coh2 Social','Coh3 Social','Coh4 Social','Coh1 Priv','Coh2 Priv','Coh3 Priv','Coh4 Priv')
legend('AutoUpdate','off')
hold on
xline([0],'--k')
yline([0],'--k')
ylabel('Magnitude (mv/s2)')
xlabel('Time (ms)')



clearvars -except data X Y Z chanlocs group_data
%% Overall data based on coherence regardless of Exp condition

overall_coh1_CPP = mean(data.overall_CPP(find(data.coherence == 3.2),:),1);
overall_coh2_CPP = mean(data.overall_CPP(find(data.coherence == 6.4),:),1);
overall_coh3_CPP = mean(data.overall_CPP(find(data.coherence == 12.8),:),1);
overall_coh4_CPP = mean(data.overall_CPP(find(data.coherence == 51.2),:),1);

data.overall_coh1_CPP = overall_coh1_CPP;
data.overall_coh2_CPP = overall_coh2_CPP;
data.overall_coh3_CPP = overall_coh3_CPP;
data.overall_coh4_CPP = overall_coh4_CPP;


% plot signals
start_plot = -200;
end_plot = 1300;
start_indx = find(data.time == start_plot);
end_indx = find(data.time == end_plot);

plot(data.time(start_indx:end_indx),data.overall_coh1_CPP(start_indx:end_indx),'color',data.color_map(1,:))
hold on
plot(data.time(start_indx:end_indx),data.overall_coh2_CPP(start_indx:end_indx),'color',data.color_map(2,:))
hold on
plot(data.time(start_indx:end_indx),data.overall_coh3_CPP(start_indx:end_indx),'color',data.color_map(3,:))
hold on
plot(data.time(start_indx:end_indx),data.overall_coh4_CPP(start_indx:end_indx),'color',data.color_map(4,:))
legend('Coh1','Coh2','Coh3','Coh4')
legend('AutoUpdate','off')
hold on
xline([0],'--k')
yline([0],'--k')
sgtitle('CPPs per coherence level, regardless of Experimental Condition')
ylabel('Magnitude (mv/s2)')
xlabel('Time (ms)')
clearvars -except data X Y Z chanlocs group_data

%% Find best window for regression fitting
time_step = 1; % in samples !!
start_time = 0; % in miliseconds !
start_indx = find(data.time == start_time);
end_time = 2000; %in milisecond !
end_indx = find(data.time == end_time);
window_size = 25; %in samples !!

b = start_indx: time_step : end_indx - window_size ;

for i = 1:length(b)-1;
time_intervals(i,1) = b(i);
time_intervals(i,2) = time_intervals(i,1) + window_size;
end

window_mean = zeros(size(data.overall_CPP,1),length(time_intervals));
for j = 1:size(data.overall_CPP,1)
for i = 1:length(time_intervals);
window_mean(j,i) = mean(data.CPP_movmean(j,time_intervals(i,1) : time_intervals(i,2)));
end
end

plot(window_mean(1,:))

model = [];
CI_coeffs = [];
coeffs_per_window = zeros(size(window_mean,2)-1,1);

for i = 1:size(window_mean,2);
model = fitlm(data.coherence,window_mean(:,i));
coeffs_per_window(i,1) = model.Coefficients.Estimate(2,1);
CI_coeffs =  coefCI(model);
coeffs_per_window(i,2) = CI_coeffs(2,1);
coeffs_per_window(i,3) = CI_coeffs(2,2);
CI_coeffs = [];
model = [];
end

max_index = find(coeffs_per_window(:,1) == max(coeffs_per_window(:,1)))


sensitive_period_time = [data.time(time_intervals(max_index,1)) data.time(time_intervals(max_index,2))];
sensitive_index = [time_intervals(max_index,1) time_intervals(max_index,2)];

ciplot(coeffs_per_window(:,2),coeffs_per_window(:,3),1:length(coeffs_per_window),'k',.2)
hold on
plot(coeffs_per_window(:,1),'color',data.color_map(1,:))
hold on
max_value = max(coeffs_per_window(:,1));
xline([find(coeffs_per_window(:,1) == max_value)],'--r')
legend('Coeff CI','Coefficient per time interval','Max Coeff')

data.coeffs_per_window = coeffs_per_window ;
data.sensitive_period_time = sensitive_period_time;
data.sensitive_index = sensitive_index;

clearvars -except data X Y Z chanlocs group_data
%% Calculate slopes per each trial on MOVMEAN CPP data
% use sensitive time interval's median value for end of regression
start_reg_time = 0;
start_indx = find(data.time == start_reg_time);
end_regres_indx = round(median(data.sensitive_index));

x_vector = 0:1:end_regres_indx - start_indx;


slope_coeffs = zeros(size(data.CPP_movmean,1),1);
model = [];
for i = 1:size(data.CPP_movmean,1);
model = fitlm(x_vector,data.CPP_movmean(i,start_indx:end_regres_indx));
slope_coeffs(i) = model.Coefficients.Estimate(2,1);
model = [];
end

slope_coeffs(:,2) = data.trial_type(1,:);
slope_coeffs(:,3) = data.coherence(:,1);
slope_coeffs(:,4) = data.exp_condition(:,1);
slope_coeffs(:,5) = data.block_num;
slope_coeffs(:,6) = data.subject;
slope_coeffs(:,7) = data.RT_EEG;
slope_coeffs(:,8) = data.Confidence;
data.slope_coeffs = slope_coeffs;



mean_g1 = mean(data.slope_coeffs(find(data.slope_coeffs(:,2) == 1),1));
std_g1 = std(data.slope_coeffs(find(data.slope_coeffs(:,2) == 1),1))/ sqrt(length(find(data.slope_coeffs(:,2) == 1)));

mean_g2 = mean(data.slope_coeffs(find(data.slope_coeffs(:,2) == 2),1));
std_g2 = std(data.slope_coeffs(find(data.slope_coeffs(:,2) == 2),1))/ sqrt(length(find(data.slope_coeffs(:,2) == 3)));

mean_g3 = mean(data.slope_coeffs(find(data.slope_coeffs(:,2) == 4),1));
std_g3 = std(data.slope_coeffs(find(data.slope_coeffs(:,2) == 4),1))/ sqrt(length(find(data.slope_coeffs(:,2) == 3)));

mean_g4 = mean(data.slope_coeffs(find(data.slope_coeffs(:,2) == 3),1));
std_g4 = std(data.slope_coeffs(find(data.slope_coeffs(:,2) == 3),1))/ sqrt(length(find(data.slope_coeffs(:,2) == 5)));

mean_g5 = mean(data.slope_coeffs(find(data.slope_coeffs(:,2) == 5),1));
std_g5 = std(data.slope_coeffs(find(data.slope_coeffs(:,2) == 5),1))/ sqrt(length(find(data.slope_coeffs(:,2) == 5)));

mean_g6 = mean(data.slope_coeffs(find(data.slope_coeffs(:,2) == 6),1));
std_g6 = std(data.slope_coeffs(find(data.slope_coeffs(:,2) == 6),1))/ sqrt(length(find(data.slope_coeffs(:,2) == 6)));

mean_g7 = mean(data.slope_coeffs(find(data.slope_coeffs(:,2) == 7),1));
std_g7 = std(data.slope_coeffs(find(data.slope_coeffs(:,2) == 7),1))/ sqrt(length(find(data.slope_coeffs(:,2) == 7)));

mean_g8 = mean(data.slope_coeffs(find(data.slope_coeffs(:,2) == 8),1));
std_g8 = std(data.slope_coeffs(find(data.slope_coeffs(:,2) == 8),1))/ sqrt(length(find(data.slope_coeffs(:,2) == 8)));


slope_by_condition_mean = [mean_g1 mean_g5; mean_g2 mean_g6; mean_g3 mean_g7; mean_g4 mean_g8];
slope_by_condition_SE = [std_g1 std_g5; std_g2 std_g6; std_g3 std_g7; std_g4 std_g8];

ngroups = size(slope_by_condition_mean, 1);
nbars = size(slope_by_condition_mean, 2);
bar(slope_by_condition_mean)
hold on
set(gca,'xticklabel',{'3.2','6.4','12.8','51.2'});
xlabel('Coherence Level')
ylabel('beta Values')
legend('Social','Private')
legend('AutoUpdate','off')
sgtitle('Regression Coefficients per Coherence Level and Condition')
hold on
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, slope_by_condition_mean(:,i), slope_by_condition_SE(:,i), '.','Color', 'k');
end

time = data.time; 

clearvars -except data X Y Z chanlocs group_data








