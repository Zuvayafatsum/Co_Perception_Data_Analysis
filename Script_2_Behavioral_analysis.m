%%Load subjects .full data and exclude trials with otulier RTs.

%find trials that RTs are neither 5000 nor 0
plot(data.RT_EEG)

hit_trials = find(~(data.RT_EEG == 5000 | data.RT_EEG == 0));

data.EEG = data.EEG(hit_trials);
data.coherence = data.coherence(hit_trials);
data.RT_Behavior = data.RT_Behavior(hit_trials);
data.Confidence = data.Confidence(hit_trials);
data.answers = data.answers(hit_trials);
data.RT_EEG = data.RT_EEG(hit_trials);
data.exp_condition = data.exp_condition(hit_trials);
data.visual_field = data.visual_field(hit_trials);
data.block_num = data.block_num(hit_trials);

plot(data.RT_EEG)

hist(data.block_num)
%% Calculate the mean and std of the remaining RTs. Exclude +- 3SD RTs from data
std_rt = std(data.RT_EEG);
mean_rt = mean(data.RT_EEG);

valid_RT = find( (mean_rt + 3*std_rt) > data.RT_EEG &  data.RT_EEG  > (mean_rt - 3*std_rt));


data.EEG = data.EEG(valid_RT);
data.coherence = data.coherence(valid_RT);
data.RT_Behavior = data.RT_Behavior(valid_RT);
data.Confidence = data.Confidence(valid_RT);
data.answers = data.answers(valid_RT);
data.RT_EEG = data.RT_EEG(valid_RT);
data.exp_condition = data.exp_condition(valid_RT);
data.visual_field = data.visual_field(valid_RT);
data.block_num = data.block_num(valid_RT);

plot(data.RT_EEG)
hist(data.block_num)

clearvars -except data group_data
%% Index trials based on coherence levels and experimental conditions
coh1_social = find(data.coherence == 3.2 & data.exp_condition == 1);
coh1_private = find(data.coherence == 3.2 & data.exp_condition == 2);
coh2_social = find(data.coherence == 6.4 & data.exp_condition == 1);
coh2_private = find(data.coherence == 6.4 & data.exp_condition == 2);
coh3_social = find(data.coherence == 12.8 & data.exp_condition == 1);
coh3_private = find(data.coherence == 12.8 & data.exp_condition == 2);
coh4_social = find(data.coherence == 51.2 & data.exp_condition == 1);
coh4_private = find(data.coherence == 51.2 & data.exp_condition == 2);

%Save trial indexes in data file to be used in EEG analysis
data.coh1_social = coh1_social;
data.coh1_private = coh1_private;
data.coh2_social = coh2_social;
data.coh2_private  = coh2_private ;
data.coh3_social = coh3_social;
data.coh3_private = coh3_private;
data.coh4_social = coh4_social;
data.coh4_private = coh4_private;
%Extract mean scores for each behavioral measure based on trial indexes

rt_social = [mean(data.RT_Behavior(coh1_social)) mean(data.RT_Behavior(coh2_social)) mean(data.RT_Behavior(coh3_social)) mean(data.RT_Behavior(coh4_social))];
rt_private = [mean(data.RT_Behavior(coh1_private)) mean(data.RT_Behavior(coh2_private)) mean(data.RT_Behavior(coh3_private)) mean(data.RT_Behavior(coh4_private))];

conf_social = [mean(data.Confidence(coh1_social)) mean(data.Confidence(coh2_social)) mean(data.Confidence(coh3_social)) mean(data.Confidence(coh4_social))];
conf_private = [mean(data.Confidence(coh1_private)) mean(data.Confidence(coh2_private)) mean(data.Confidence(coh3_private)) mean(data.Confidence(coh4_private))];

acc_social = [sum((data.answers(coh1_social)) == 2)/length(data.answers(coh1_social)) sum((data.answers(coh2_social)) == 2)/length(data.answers(coh2_social)).......
    sum((data.answers(coh3_social)) == 2)/length(data.answers(coh3_social)) sum((data.answers(coh4_social)) == 2)/length(data.answers(coh4_social)) ];
acc_private = [sum((data.answers(coh1_private)) == 2)/length(data.answers(coh1_private)) sum((data.answers(coh2_private)) == 2)/length(data.answers(coh2_private)).......
    sum((data.answers(coh3_private)) == 2)/length(data.answers(coh3_private)) sum((data.answers(coh4_private)) == 2)/length(data.answers(coh4_private)) ];


% Visualize the Results
figure(1)
set(0, 'DefaultLineLineWidth', 1);
subplot(221)
plot(1:4,rt_social,'-o')
hold on
plot(1:4,rt_private,'-o')
legend('Social','Private','Location','North')
xlabel('Coherence Level')
ylabel('Mean Reaction Times')
subplot(222)
plot(1:4,conf_social,'-o')
hold on
plot(1:4,conf_private,'-o')
legend('Social','Private','Location','North')
xlabel('Coherence Level')
ylabel('Mean Confidence Judgements')
subplot(223)
plot(1:4,acc_social,'-o')
hold on
plot(1:4,acc_private,'-o')
legend('Social','Private','Location','North')
xlabel('Coherence Level')
ylabel('Mean Accuracy (%)')
sgtitle(char(data.subjects_ID))

%% Save the data set
clearvars -except data group_data
%% save data for group behavioral analysis
group_data = {};
group_data.coherence = [];



% change per subject



group_data.color_map = data.color_map;
group_data.time = data.time;
group_data.left_chans = data.left_chans;
group_data.right_chans = data.right_chans;
group_data.all_cp_chans = data.all_cp_chans;
group_data.Fsample = data.Fsample;


subject_num = 5;
trial_length = length(group_data.coherence);

for i = (trial_length +1) : (trial_length + length(data.coherence))
    
group_data.RT_EEG(i,1) = data.RT_EEG(i-trial_length,1);    
group_data.coherence(i,1) = data.coherence(i-trial_length,1);
group_data.Confidence(i,1) = data.Confidence(i-trial_length,1);
group_data.answers(i,1) = data.answers(i-trial_length,1);
group_data.exp_condition(i,1)  = data.exp_condition(i-trial_length,1);
group_data.block_num(i,1) = data.block_num(i-trial_length,1);
group_data.visual_field(i,1) = data.visual_field(i-trial_length,1);
group_data.subject(i,1) = subject_num;
group_data.EEG(:,:,i)= data.EEG(:,:,i- trial_length);
group_data.Laplacian_filtered(:,:,i) = data.Laplacian_filtered(:,:,i-trial_length);
group_data.bc_corrected(:,:,i) = data.bc_corrected(:,:,i - trial_length);
end


clearvars -except  group_data




group_data.coh1_social = find(group_data.coherence == 3.2 & group_data.exp_condition == 1);
group_data.coh1_private = find(group_data.coherence == 3.2 & group_data.exp_condition == 2);
group_data.coh2_social = find(group_data.coherence == 6.4 & group_data.exp_condition == 1);
group_data.coh2_private = find(group_data.coherence == 6.4 & group_data.exp_condition == 2);
group_data.coh3_social = find(group_data.coherence == 12.8 & group_data.exp_condition == 1);
group_data.coh3_private = find(group_data.coherence == 12.8 & group_data.exp_condition == 2);
group_data.coh4_social = find(group_data.coherence == 51.2 & group_data.exp_condition == 1);
group_data.coh4_private = find(group_data.coherence == 51.2 & group_data.exp_condition == 2);


%% Calculate behavioral metrics and Visualize group data
%for RTs
rt_social = [mean(group_data.RT(group_data.coh1_social)) mean(group_data.RT(group_data.coh2_social)) mean(group_data.RT(group_data.coh3_social)) mean(group_data.RT(group_data.coh4_social))];
rt_private = [mean(group_data.RT(group_data.coh1_private)) mean(group_data.RT(group_data.coh2_private)) mean(group_data.RT(group_data.coh3_private)) mean(group_data.RT(group_data.coh4_private))];
%for Confidences
conf_social = [mean(group_data.confidence(group_data.coh1_social)) mean(group_data.confidence(group_data.coh2_social)) mean(group_data.confidence(group_data.coh3_social)) mean(group_data.confidence(group_data.coh4_social))];
conf_private = [mean(group_data.confidence(group_data.coh1_private)) mean(group_data.confidence(group_data.coh2_private)) mean(group_data.confidence(group_data.coh3_private)) mean(group_data.confidence(group_data.coh4_private))];
%for Accuracy
acc_social = [sum((group_data.answers(group_data.coh1_social)) == 2)/length(group_data.answers(group_data.coh1_social)) sum((group_data.answers(group_data.coh2_social)) == 2)/length(group_data.answers(group_data.coh2_social)).......
    sum((group_data.answers(group_data.coh3_social)) == 2)/length(group_data.answers(group_data.coh3_social)) sum((group_data.answers(group_data.coh4_social)) == 2)/length(group_data.answers(group_data.coh4_social)) ];
acc_private = [sum((group_data.answers(group_data.coh1_private)) == 2)/length(group_data.answers(group_data.coh1_private)) sum((group_data.answers(group_data.coh2_private)) == 2)/length(group_data.answers(group_data.coh2_private)).......
    sum((group_data.answers(group_data.coh3_private)) == 2)/length(group_data.answers(group_data.coh3_private)) sum((group_data.answers(group_data.coh4_private)) == 2)/length(group_data.answers(group_data.coh4_private)) ];


%Visualize 
figure(1)
set(0, 'DefaultLineLineWidth', 1);
subplot(221)
plot(1:4,rt_social,'-o')
hold on
plot(1:4,rt_private,'-o')
legend('Social','Private','Location','North')
xlabel('Coherence Level')
ylabel('Mean Reaction Times')
subplot(222)
plot(1:4,conf_social,'-o')
hold on
plot(1:4,conf_private,'-o')
legend('Social','Private','Location','North')
xlabel('Coherence Level')
ylabel('Mean Confidence Judgements')
subplot(223)
plot(1:4,acc_social,'-o')
hold on
plot(1:4,acc_private,'-o')
legend('Social','Private','Location','North')
xlabel('Coherence Level')
ylabel('Mean Accuracy (%)')
sgtitle('Group Means')


%% Testing for learning effects
% rt by block number
rt_by_block = fitlm(group_data.block_num,group_data.RT);
rt_by_block.Coefficients
plot(rt_by_block)

% rt by order AND exp_condition

tbl = table(group_data.RT,group_data.block_num,group_data.exp_condition,'VariableNames',{'RT','Block','Condition'});
rt_by_block_and_cond = fitlm(tbl,'RT~Block+Condition');
rt_by_block_and_cond.Coefficients

% rt by order AND exp_condition AND Interaction
rt_by_block_and_cond_and_interaction = fitlm(tbl,'RT ~ Block + Condition + Block*Condition ');
rt_by_block_and_cond_and_interaction.Coefficients



