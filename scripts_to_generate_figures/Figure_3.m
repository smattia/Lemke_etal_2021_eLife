%% SCRIPT TO GENERATE FIGURE 3 & FIGURE 3 FIGURE SUPPLEMENTS
% 
% to generate figures, download data from https://doi.org/10.7272/Q6KK9927

clc; clear all; close all;
data_path = 'D:\Lemke_etal_2021\Lemke_eLife_2021_data';

%% FIG 3 A | EXAMPLE HIGH LFP COHERENCE PAIRS ACROSS THREE SESSION

    day_blocks = {[1 2],[3 4 5 6 8],[9 10 11 12 13 14 15],[16 17 18 19 20 21 22 24 25],[26 27 30 31 36 39],[41 42 43 47 49 51 52]};
    pre_sleep_blocks = {[1],[1],[1],[3],[1],[1]};
    post_sleep_blocks = {[2],[4],[6],[8],[5],[4]};

    load([data_path '\LC4\sleep_activity\LFP_coherence.mat'])

    fig = figure;

    for m1_chan_ind = 1:16

        subplot(1,3,1)
        hold on;
        pair_count = 1;
        for m1_chan = 1:16
            for dls_chan = 1:16
                if mean(pre_LFP_coh{pair_count}(3,48:113))>0.65 & m1_chan==m1_chan_ind
                    plot([m1_chan dls_chan],[10 1],'LineWidth', 1,'color',[0 0 0])
                end
                pair_count = pair_count + 1;
            end
        end
        xlim([0 17])

        subplot(1,3,2)
        hold on;
        pair_count = 1;
        for m1_chan = 1:16
            for dls_chan = 1:16
                if mean(post_LFP_coh{pair_count}(3,48:113))>0.65 & m1_chan==m1_chan_ind
                    plot([m1_chan dls_chan],[10 1],'LineWidth', 1,'color',[0 0 0])
                end
                pair_count = pair_count + 1;
            end
        end
        xlim([0 17])

        subplot(1,3,3)
        hold on;
        pair_count = 1;
        for m1_chan = 1:16
            for dls_chan = 1:16
                if mean(pre_LFP_coh{pair_count}(4,48:113))>0.65 & m1_chan==m1_chan_ind
                    plot([m1_chan dls_chan],[10 1],'LineWidth', 1,'color',[0 0 0])
                end
                pair_count = pair_count + 1;
            end
        end
        xlim([0 17])

    end
  
%% FIG 3 B | EXAMPLE COMPARISON OF LFP COHERENCE SPECTRUMS ACROSS TWO SLEEP SESSIONS

    load([data_path '\LC3\sleep_activity\LFP_coherence.mat'])

    coh_1 = [];
    coh_2 = [];
    coh_3 = [];
    coh_4 = [];
    for pairs = [18 32 41 49 50 54 57 60 63 73]
        coh_1 = [coh_1; pre_LFP_coh{pairs}(2,:)];
        coh_2 = [coh_2; post_LFP_coh{pairs}(2,:)];
        coh_3 = [coh_3; pre_LFP_coh{pairs}(5,:)];
        coh_4 = [coh_4; post_LFP_coh{pairs}(5,:)];
    end
    fig = figure
    hold on
    plot(mean(coh_1),'color',[0 0 0]);
    plot(mean(coh_2),'color',[0 0 0]);
    plot(mean(coh_3),'color',[1 0 0]);
    plot(mean(coh_4),'color',[1 0 0]);   
    
%% FIG 3 E | sFIG 1 | TIME COURSE OF COHERENCE IN ALL INDIVIDUAL ANIMALS

all_animal_coh = [];
all_animal_corr = [];

%%% LC3

    %%% BEHAVIOR

        xvel = cell(1,8);
        yvel = cell(1,8);
        day_count = 1;

        for day = 1:8

            files = dir([data_path '\LC3\reach_trajectories\day_' num2str(day) '\*.csv']);

            %%% find mean pellet position across trials
            pellet_position = [];
            for file = 1:length(files)
                trial_kin = readmatrix([data_path '\LC3\reach_trajectories\day_' num2str(day) '\' files(file).name]);
                tmp_pellet_x = trial_kin(1:25,5);
                tmp_pellet_y = trial_kin(1:25,6);
                tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
                pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                pellet_position = [pellet_position; pellet_location];
            end
            pellet_position = nanmean(pellet_position);

            %%% get velocity profiles time locked to frame with paw closest to pellet
            all_diffx = [];
            all_diffy = [];

            for file = 1:length(files)

                trial_kin = readmatrix([data_path '\LC3\reach_trajectories\day_' num2str(day) '\' files(file).name]);

                %%% get paw trajectory
                tmp_x = trial_kin(:,2);
                tmp_y = trial_kin(:,3);
                tmp_x(trial_kin(:,4)<.99) = NaN;
                tmp_y(trial_kin(:,4)<.99) = NaN;

                %%% get paw velocity
                diff_x = diff(tmp_x);
                diff_y = diff(tmp_y);

                %%% find frame with paw closest to pellet
                [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));

                %%% skip bad trials
                if ~(i<31 || i>size(trial_kin,1)-31)
                    all_diffx = [all_diffx; diff_x(i-30:i+30)'];
                    all_diffy = [all_diffy; diff_y(i-30:i+30)'];
                end

            end

            xvel{day_count} = all_diffx;
            yvel{day_count} = all_diffy;
            day_count = day_count + 1;

        end

        %%% calculate correlation between each day's trial-averaged velocity profile and trial-averaged velocity profile on last day
        x_corr = [];
        y_corr = [];
        for n = 1:length(xvel)
            x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
            y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
        end
        animal_corr = mean([x_corr; y_corr]);

    %%% SLEEP COHERENCE

        load([data_path '\LC3\sleep_activity\LFP_coherence.mat'])

        increase_channels = [];
        for pairs=1:256
            if mean(mean([pre_LFP_coh{pairs}(8,48:113); post_LFP_coh{pairs}(8,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)])))+0.025
                increase_channels = [increase_channels 1];
            else
                increase_channels = [increase_channels 0];
            end
        end

        day_pre = [];
        day_post = [];
        for pairs = find(increase_channels)
            tmp_pre = [];
            tmp_post = [];
            for day=[1:8]
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,48:113))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,48:113))];
            end
            day_pre = [day_pre; tmp_pre];
            day_post = [day_post; tmp_post];
        end

    %%% PLOT

        figure
        for day=1:8
            
            subplot(2,1,1); hold on;
                errorbar(day-0.1,mean(day_pre(:,day)),std(day_pre(:,day))/sqrt(length(day_pre(:,day))),'color',[0.5 0.5 0.5])
                errorbar(day+0.1,mean(day_post(:,day)),std(day_post(:,day))/sqrt(length(day_post(:,day))),'color',[0 0 0])
                plot([day-0.1 day+0.1],[mean(day_pre(:,day)) mean(day_post(:,day))],'color',[0.5 0.5 0.5])
                if day<8
                    plot([day+0.1 day+0.9],[mean(day_post(:,day)) mean(day_pre(:,day+1))],'color',[0 0 0])
                end
                xlim([0 9])
                
            subplot(2,1,2); hold on;
                scatter(day,animal_corr(day),50,[1 0 0],'filled')
                xlim([0 9])
        end

        all_animal_coh = [all_animal_coh zscore(mean([day_pre; day_post]))];
        all_animal_corr = [all_animal_corr zscore(animal_corr)];
        
%%% LC4

    %%% BEHAVIOR

        xvel = cell(1,5);
        yvel = cell(1,5);
        day_count = 1;

        for day = 2:6

            files = dir([data_path '\LC4\reach_trajectories\day_' num2str(day) '\*.csv']);

            %%% find mean pellet position across trials
            pellet_position = [];
            for file = 1:length(files)
                trial_kin = readmatrix([data_path '\LC4\reach_trajectories\day_' num2str(day) '\' files(file).name]);
                tmp_pellet_x = trial_kin(1:25,5);
                tmp_pellet_y = trial_kin(1:25,6);
                tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
                pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                pellet_position = [pellet_position; pellet_location];
            end
            pellet_position = nanmean(pellet_position);

            %%% get velocity profiles time locked to frame with paw closest to pellet
            all_diffx = [];
            all_diffy = [];

            for file = 1:length(files)

                trial_kin = readmatrix([data_path '\LC4\reach_trajectories\day_' num2str(day) '\' files(file).name]);

                %%% get paw trajectory
                tmp_x = trial_kin(:,2);
                tmp_y = trial_kin(:,3);
                tmp_x(trial_kin(:,4)<.99) = NaN;
                tmp_y(trial_kin(:,4)<.99) = NaN;

                %%% get paw velocity
                diff_x = diff(tmp_x);
                diff_y = diff(tmp_y);

                %%% find frame with paw closest to pellet
                [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));

                %%% skip bad trials
                if ~(i<31 || i>size(trial_kin,1)-31)
                    all_diffx = [all_diffx; diff_x(i-30:i+30)'];
                    all_diffy = [all_diffy; diff_y(i-30:i+30)'];
                end

            end

            xvel{day_count} = all_diffx;
            yvel{day_count} = all_diffy;
            day_count = day_count + 1;

        end

        %%% calculate correlation between each day's trial-averaged velocity profile and trial-averaged velocity profile on last day
        x_corr = [];
        y_corr = [];
        for n = 1:length(xvel)
            x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
            y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
        end
        animal_corr = mean([x_corr; y_corr]);

    %%% SLEEP COHERENCE

        load([data_path '\LC4\sleep_activity\LFP_coherence.mat'])

        increase_channels = [];
        for pairs=1:256
            if mean(mean([pre_LFP_coh{pairs}(6,48:113); post_LFP_coh{pairs}(6,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(2,48:113); post_LFP_coh{pairs}(2,48:113)])))+0.025
                increase_channels = [increase_channels 1];
            else
                increase_channels = [increase_channels 0];
            end
        end

        day_pre = [];
        day_post = [];
        for pairs = find(increase_channels)
            tmp_pre = [];
            tmp_post = [];
            for day=[2:6]
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,48:113))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,48:113))];
            end
            day_pre = [day_pre; tmp_pre];
            day_post = [day_post; tmp_post];
        end

    %%% PLOT

        figure
        for day=1:5
            
            subplot(2,1,1); hold on;
                errorbar(day-0.1,mean(day_pre(:,day)),std(day_pre(:,day))/sqrt(length(day_pre(:,day))),'color',[0.5 0.5 0.5])
                errorbar(day+0.1,mean(day_post(:,day)),std(day_post(:,day))/sqrt(length(day_post(:,day))),'color',[0 0 0])
                plot([day-0.1 day+0.1],[mean(day_pre(:,day)) mean(day_post(:,day))],'color',[0.5 0.5 0.5])
                if day<5
                    plot([day+0.1 day+0.9],[mean(day_post(:,day)) mean(day_pre(:,day+1))],'color',[0 0 0])
                end
                xlim([0 6])
                
            subplot(2,1,2); hold on;
                scatter(day,animal_corr(day),50,[1 0 0],'filled')
                xlim([0 6])
        end
        
        all_animal_coh = [all_animal_coh zscore(mean([day_pre; day_post]))];
        all_animal_corr = [all_animal_corr zscore(animal_corr)];
        
%%% LC6

    %%% BEHAVIOR

        xvel = cell(1,9);
        yvel = cell(1,9);
        reach_blocks = {[],[7],[12 13],[17 18],[21 22],[26],[29],[33 34],[38 39],[43 44]};
        day_count = 1;

        for day = 2:10
    
            files = dir([data_path '\LC6\reach_trajectories\day_' num2str(day) '\*.csv']);
            files = sort_nat({files.name});
            vid_timestamp = dir([data_path '\LC6\reach_trajectories\day_' num2str(day) '\*.mat']); 
            vid_timestamp = sort_nat({vid_timestamp.name});

            all_diffx = [];
            all_diffy = [];

            for block = 1:length(files)

                %%% load all trial kinematics
                kin_all = readmatrix([data_path '\LC6\reach_trajectories\day_' num2str(day) '\' files{block}]);
                load([data_path '\LC6\reach_trajectories\day_' num2str(day) '\' vid_timestamp{block}]);

                %%% find individual trials 
                [~,i] = find(Vid.scalars.Vid0.data>1e8);
                for n_i = 1:length(i)
                    Vid.scalars.Vid0.data(i(n_i)) = round(mean(Vid.scalars.Vid0.data(i(n_i)-1):Vid.scalars.Vid0.data(i(n_i)+1)));
                end
                Vid.scalars.Vid0.data = Vid.scalars.Vid0.data+1;
                first_frames = find(diff(Vid.scalars.Vid0.ts)>8);
                trial_starts_idx = [Vid.scalars.Vid0.data(1) Vid.scalars.Vid0.data(first_frames+1)];
                trial_ends_idx = [Vid.scalars.Vid0.data(first_frames) Vid.scalars.Vid0.data(end)];

                %%% find mean pellet position across trials
                pellet_position = [];
                for n_trial = 1:length(trial_starts_idx)-1

                    %%% get timesstamps (ts) for frames in trial
                    trial_frames_idx = Vid.scalars.Vid0.data(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));
                    trial_frames_ts = Vid.scalars.Vid0.ts(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));

                    %%% interpolate to get missing ts
                    interp_trial_ts = interp1(trial_frames_idx,trial_frames_ts,trial_starts_idx(n_trial):trial_ends_idx(n_trial));

                    %%% get trial kinematics
                    trial_kin = kin_all(trial_starts_idx(n_trial):trial_ends_idx(n_trial),:);

                    %%% get pellet location from first 25 frames
                    tmp_pellet_x = trial_kin(1:25,5);
                    tmp_pellet_y = trial_kin(1:25,6);
                    tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
                    tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
                    pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                    pellet_position = [pellet_position; pellet_location];

                end
                pellet_position = nanmean(pellet_position);

                for n_trial = 1:length(trial_starts_idx)-1

                    %%% get timesstamps (ts) for frames in trial
                    trial_frames_idx = Vid.scalars.Vid0.data(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));
                    trial_frames_ts = Vid.scalars.Vid0.ts(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));

                    %%% interpolate to get missing ts
                    interp_trial_ts = interp1(trial_frames_idx,trial_frames_ts,trial_starts_idx(n_trial):trial_ends_idx(n_trial));

                    %%% get trial kinematics
                    trial_kin = kin_all(trial_starts_idx(n_trial):trial_ends_idx(n_trial),:);

                    %%% get paw trajectory
                    tmp_x = trial_kin(:,2);
                    tmp_y = trial_kin(:,3);
                    tmp_x(trial_kin(:,4)<.99) = NaN;
                    tmp_y(trial_kin(:,4)<.99) = NaN;

                    %%% get paw velocity
                    diff_x = diff(tmp_x);
                    diff_y = diff(tmp_y);

                    %%% find frame with paw closest to pellet
                    [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));

                    %%% skip bad trials
                    if ~(i<31 || i>size(trial_kin,1)-31)
                        all_diffx = [all_diffx; diff_x(i-30:i+30)'];
                        all_diffy = [all_diffy; diff_y(i-30:i+30)'];
                    end

                end
            end

            xvel{day_count} = all_diffx;
            yvel{day_count} = all_diffy;
            day_count = day_count + 1;

        end
        
        %%% calculate correlation between each day's trial-averaged velocity profile and trial-averaged velocity profile on last day
        x_corr = [];
        y_corr = [];
        for n = 1:length(xvel)
            x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
            y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
        end
        animal_corr = mean([x_corr; y_corr]);

    %%% SLEEP COHERENCE

        load([data_path '\LC6\sleep_activity\LFP_coherence.mat'])

        increase_channels = [];
        for pairs=1:256
            if mean(mean([pre_LFP_coh{pairs}(9,48:113); post_LFP_coh{pairs}(9,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)])))+0.025
                increase_channels = [increase_channels 1];
            else
                increase_channels = [increase_channels 0];
            end
        end

        day_pre = [];
        day_post = [];
        for pairs = find(increase_channels)
            tmp_pre = [];
            tmp_post = [];
            for day=[1:9]
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,48:113))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,48:113))];
            end
            day_pre = [day_pre; tmp_pre];
            day_post = [day_post; tmp_post];
        end

    %%% PLOT

        figure
        for day=1:9
            
            subplot(2,1,1); hold on;
                errorbar(day-0.1,mean(day_pre(:,day)),std(day_pre(:,day))/sqrt(length(day_pre(:,day))),'color',[0.5 0.5 0.5])
                errorbar(day+0.1,mean(day_post(:,day)),std(day_post(:,day))/sqrt(length(day_post(:,day))),'color',[0 0 0])
                plot([day-0.1 day+0.1],[mean(day_pre(:,day)) mean(day_post(:,day))],'color',[0.5 0.5 0.5])
                if day<9
                    plot([day+0.1 day+0.9],[mean(day_post(:,day)) mean(day_pre(:,day+1))],'color',[0 0 0])
                end
                xlim([0 10])
                
            subplot(2,1,2); hold on;
                scatter(day,animal_corr(day),50,[1 0 0],'filled')
                xlim([0 10])
        end

        all_animal_coh = [all_animal_coh zscore(mean([day_pre; day_post]))];
        all_animal_corr = [all_animal_corr zscore(animal_corr)];
        
%%% LC5

    %%% BEHAVIOR

        xvel = cell(1,14);
        yvel = cell(1,14);
        reach_blocks = {[],[],[],[9],[12],[16],[19],[22 23 24],[28],[31],[35],[39],[42],[46],[50],[54],[59]};
        day_count = 1;
 
        for day = 4:17
    
            files = dir([data_path '\LC5\reach_trajectories\day_' num2str(day) '\*.csv']);
            files = sort_nat({files.name});
            vid_timestamp = dir([data_path '\LC5\reach_trajectories\day_' num2str(day) '\*.mat']); 
            vid_timestamp = sort_nat({vid_timestamp.name});

            all_diffx = [];
            all_diffy = [];

            for block = 1:length(files)

                %%% load all trial kinematics
                kin_all = readmatrix([data_path '\LC5\reach_trajectories\day_' num2str(day) '\' files{block}]);
                load([data_path '\LC5\reach_trajectories\day_' num2str(day) '\' vid_timestamp{block}]);

                %%% find individual trials 
                [~,i] = find(Vid.scalars.Vid0.data>1e8);
                for n_i = 1:length(i)
                    Vid.scalars.Vid0.data(i(n_i)) = round(mean(Vid.scalars.Vid0.data(i(n_i)-1):Vid.scalars.Vid0.data(i(n_i)+1)));
                end
                Vid.scalars.Vid0.data = Vid.scalars.Vid0.data+1;
                first_frames = find(diff(Vid.scalars.Vid0.ts)>8);
                trial_starts_idx = [Vid.scalars.Vid0.data(1) Vid.scalars.Vid0.data(first_frames+1)];
                trial_ends_idx = [Vid.scalars.Vid0.data(first_frames) Vid.scalars.Vid0.data(end)];

                %%% find mean pellet position across trials
                pellet_position = [];
                for n_trial = 1:length(trial_starts_idx)-1

                    %%% get timesstamps (ts) for frames in trial
                    trial_frames_idx = Vid.scalars.Vid0.data(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));
                    trial_frames_ts = Vid.scalars.Vid0.ts(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));

                    %%% interpolate to get missing ts
                    interp_trial_ts = interp1(trial_frames_idx,trial_frames_ts,trial_starts_idx(n_trial):trial_ends_idx(n_trial));

                    %%% get trial kinematics
                    trial_kin = kin_all(trial_starts_idx(n_trial):trial_ends_idx(n_trial),:);

                    %%% get pellet location from first 25 frames
                    tmp_pellet_x = trial_kin(1:25,5);
                    tmp_pellet_y = trial_kin(1:25,6);
                    tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
                    tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
                    pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                    pellet_position = [pellet_position; pellet_location];

                end
                pellet_position = nanmean(pellet_position);

                for n_trial = 1:length(trial_starts_idx)-1

                    %%% get timesstamps (ts) for frames in trial
                    trial_frames_idx = Vid.scalars.Vid0.data(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));
                    trial_frames_ts = Vid.scalars.Vid0.ts(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));

                    %%% interpolate to get missing ts
                    interp_trial_ts = interp1(trial_frames_idx,trial_frames_ts,trial_starts_idx(n_trial):trial_ends_idx(n_trial));

                    %%% get trial kinematics
                    trial_kin = kin_all(trial_starts_idx(n_trial):trial_ends_idx(n_trial),:);

                    %%% get paw trajectory
                    tmp_x = trial_kin(:,2);
                    tmp_y = trial_kin(:,3);
                    tmp_x(trial_kin(:,4)<.99) = NaN;
                    tmp_y(trial_kin(:,4)<.99) = NaN;

                    %%% get paw velocity
                    diff_x = diff(tmp_x);
                    diff_y = diff(tmp_y);

                    %%% find frame with paw closest to pellet
                    [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));

                    %%% skip bad trials
                    if ~(i<31 || i>size(trial_kin,1)-31)
                        all_diffx = [all_diffx; diff_x(i-30:i+30)'];
                        all_diffy = [all_diffy; diff_y(i-30:i+30)'];
                    end

                end
            end

            xvel{day_count} = all_diffx;
            yvel{day_count} = all_diffy;
            day_count = day_count + 1;

        end
      
        %%% calculate correlation between each day's trial-averaged velocity profile and trial-averaged velocity profile on last day
        x_corr = [];
        y_corr = [];
        for n = 1:length(xvel)
            x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
            y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
        end
        animal_corr = mean([x_corr; y_corr]);

    %%% SLEEP COHERENCE

        load([data_path '\LC5\sleep_activity\LFP_coherence.mat'])

        increase_channels = [];
        for pairs=1:512
            if mean(mean([pre_LFP_coh{pairs}(15,48:113); post_LFP_coh{pairs}(15,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(2,48:113); post_LFP_coh{pairs}(2,48:113)])))+0.025
                increase_channels = [increase_channels 1];
            else
                increase_channels = [increase_channels 0];
            end
        end

        day_pre = [];
        day_post = [];
        for pairs = find(increase_channels)
            tmp_pre = [];
            tmp_post = [];
            for day=[2:15]
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,48:113))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,48:113))];
            end
            day_pre = [day_pre; tmp_pre];
            day_post = [day_post; tmp_post];
        end

    %%% PLOT

        figure
        for day=1:14
            
            subplot(2,1,1); hold on;
                errorbar(day-0.1,mean(day_pre(:,day)),std(day_pre(:,day))/sqrt(length(day_pre(:,day))),'color',[0.5 0.5 0.5])
                errorbar(day+0.1,mean(day_post(:,day)),std(day_post(:,day))/sqrt(length(day_post(:,day))),'color',[0 0 0])
                plot([day-0.1 day+0.1],[mean(day_pre(:,day)) mean(day_post(:,day))],'color',[0.5 0.5 0.5])
                if day<14
                    plot([day+0.1 day+0.9],[mean(day_post(:,day)) mean(day_pre(:,day+1))],'color',[0 0 0])
                end
                xlim([0 15])
                
            subplot(2,1,2); hold on;
                scatter(day,animal_corr(day),50,[1 0 0],'filled')
                xlim([0 15])
        end

        all_animal_coh = [all_animal_coh zscore(mean([day_pre; day_post]))];
        all_animal_corr = [all_animal_corr zscore(animal_corr)];
        
%%% LC2

    %%% BEHAVIOR

        xvel = cell(1,10);
        yvel = cell(1,10);
        day_count = 1;
        for day = 4:13

            files = dir([data_path '\LC2\reach_trajectories\day_' num2str(day) '\*.csv']);

            %%% find mean pellet position across trials
            pellet_position = [];
            for file = 1:length(files)
                trial_kin = readmatrix([data_path '\LC2\reach_trajectories\day_' num2str(day) '\' files(file).name]);
                tmp_pellet_x = trial_kin(1:25,8);
                tmp_pellet_y = trial_kin(1:25,9);
                tmp_pellet_x(trial_kin(1:25,10)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,10)<.99) = NaN;
                pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                pellet_position = [pellet_position; pellet_location];
            end
            pellet_position = nanmean(pellet_position);

            %%% get velocity profiles time locked to frame with paw closest to pellet
            all_diffx = [];
            all_diffy = [];

            for file = 1:length(files)

                trial_kin = readmatrix([data_path '\LC2\reach_trajectories\day_' num2str(day) '\' files(file).name]);

                %%% get paw trajectory
                tmp_x = trial_kin(:,2);
                tmp_y = trial_kin(:,3);
                tmp_x(trial_kin(:,4)<.99) = NaN;
                tmp_y(trial_kin(:,4)<.99) = NaN;

                %%% get paw velocity
                diff_x = diff(tmp_x);
                diff_y = diff(tmp_y);

                %%% find frame with paw closest to pellet
                [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));

                %%% skip bad trials
                if ~(i<31 || i>size(trial_kin,1)-31)
                    all_diffx = [all_diffx; diff_x(i-30:i+30)'];
                    all_diffy = [all_diffy; diff_y(i-30:i+30)'];
                end

            end

            xvel{day_count} = all_diffx;
            yvel{day_count} = all_diffy;
            day_count = day_count + 1;

        end
        
        %%% calculate correlation between each day's trial-averaged velocity profile and trial-averaged velocity profile on last day
        x_corr = [];
        y_corr = [];
        for n = 1:length(xvel)
            x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
            y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
        end
        animal_corr = mean([x_corr; y_corr]);

    %%% SLEEP COHERENCE

        load([data_path '\LC2\sleep_activity\LFP_coherence.mat'])

        increase_channels = [];
        for pairs=1:1024
            if mean(mean([pre_LFP_coh{pairs}(10,40:95); post_LFP_coh{pairs}(10,40:95)])) > (mean(mean([pre_LFP_coh{pairs}(1,40:95); post_LFP_coh{pairs}(1,40:95)])))+0.025
                increase_channels = [increase_channels 1];
            else
                increase_channels = [increase_channels 0];
            end
        end

        day_pre = [];
        day_post = [];
        for pairs = find(increase_channels)
            tmp_pre = [];
            tmp_post = [];
            for day=1:10
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,40:95))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,40:95))];
            end
            day_pre = [day_pre; tmp_pre];
            day_post = [day_post; tmp_post];
        end

    %%% PLOT

        figure
        for day=1:10
            
            subplot(2,1,1); hold on;
                errorbar(day-0.1,mean(day_pre(:,day)),std(day_pre(:,day))/sqrt(length(day_pre(:,day))),'color',[0.5 0.5 0.5])
                errorbar(day+0.1,mean(day_post(:,day)),std(day_post(:,day))/sqrt(length(day_post(:,day))),'color',[0 0 0])
                plot([day-0.1 day+0.1],[mean(day_pre(:,day)) mean(day_post(:,day))],'color',[0.5 0.5 0.5])
                if day<10
                    plot([day+0.1 day+0.9],[mean(day_post(:,day)) mean(day_pre(:,day+1))],'color',[0 0 0])
                end
                xlim([0 11])
                
            subplot(2,1,2); hold on;
                scatter(day,animal_corr(day),50,[1 0 0],'filled')
                xlim([0 11])
        end

        all_animal_coh = [all_animal_coh zscore(mean([day_pre; day_post]))];
        all_animal_corr = [all_animal_corr zscore(animal_corr)];
        
%%% LC1

    %%% BEHAVIOR

        xvel = cell(1,8);
        yvel = cell(1,8);
        day_count = 1;
        for day = [1:8]

            files = dir([data_path '\LC1\reach_trajectories\day_' num2str(day) '\*.csv']);

            %%% find mean pellet position across trials
            pellet_position = [];
            for file = 1:length(files)
                trial_kin = readmatrix([data_path '\LC1\reach_trajectories\day_' num2str(day) '\' files(file).name]);
                tmp_pellet_x = trial_kin(1:25,8);
                tmp_pellet_y = trial_kin(1:25,9);
                tmp_pellet_x(trial_kin(1:25,10)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,10)<.99) = NaN;
                pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                pellet_position = [pellet_position; pellet_location];
            end
            pellet_position = nanmean(pellet_position);

            %%% get velocity profiles time locked to frame with paw closest to pellet
            all_diffx = [];
            all_diffy = [];

            for file = 1:length(files)

                trial_kin = readmatrix([data_path '\LC1\reach_trajectories\day_' num2str(day) '\' files(file).name]);

                %%% get paw trajectory
                tmp_x = trial_kin(:,2);
                tmp_y = trial_kin(:,3);
                tmp_x(trial_kin(:,4)<.99) = NaN;
                tmp_y(trial_kin(:,4)<.99) = NaN;

                %%% get paw velocity
                diff_x = diff(tmp_x);
                diff_y = diff(tmp_y);

                %%% find frame with paw closest to pellet
                [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));

                %%% skip bad trials
                if ~(i<31 || i>size(trial_kin,1)-31)
                    all_diffx = [all_diffx; diff_x(i-30:i+30)'];
                    all_diffy = [all_diffy; diff_y(i-30:i+30)'];
                end

            end

            xvel{day_count} = all_diffx;
            yvel{day_count} = all_diffy;
            day_count = day_count + 1;

        end
        
        %%% calculate correlation between each day's trial-averaged velocity profile and trial-averaged velocity profile on last day
        x_corr = [];
        y_corr = [];
        for n = 1:length(xvel)
            x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
            y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
        end
        animal_corr = mean([x_corr; y_corr]);

    %%% SLEEP COHERENCE

        load([data_path '\LC1\sleep_activity\LFP_coherence.mat'])

        increase_channels = [];
        for pairs=1:1406
            if mean(mean([pre_LFP_coh{pairs}(7,40:95); post_LFP_coh{pairs}(7,40:95)])) > (mean(mean([pre_LFP_coh{pairs}(1,40:95); post_LFP_coh{pairs}(1,40:95)])))+0.025
                increase_channels = [increase_channels 1];
            else
                increase_channels = [increase_channels 0];
            end
        end

        day_pre = [];
        day_post = [];
        for pairs = find(increase_channels)
            tmp_pre = [];
            tmp_post = [];
            for day=[1:7]
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,40:95))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,40:95))];
            end
            day_pre = [day_pre; tmp_pre];
            day_post = [day_post; tmp_post];
        end

    %%% PLOT

        figure
        for day=1:8
            
            if day<8
            subplot(2,1,1); hold on;
                errorbar(day-0.1,mean(day_pre(:,day)),std(day_pre(:,day))/sqrt(length(day_pre(:,day))),'color',[0.5 0.5 0.5])
                errorbar(day+0.1,mean(day_post(:,day)),std(day_post(:,day))/sqrt(length(day_post(:,day))),'color',[0 0 0])
                plot([day-0.1 day+0.1],[mean(day_pre(:,day)) mean(day_post(:,day))],'color',[0.5 0.5 0.5])
                if day<7
                    plot([day+0.1 day+0.9],[mean(day_post(:,day)) mean(day_pre(:,day+1))],'color',[0 0 0])
                end
                xlim([0 9])
            end
                
            subplot(2,1,2); hold on;
                scatter(day,animal_corr(day),50,[1 0 0],'filled')
                xlim([0 9])
        end

        all_animal_coh = [all_animal_coh zscore(mean([day_pre; day_post]))];
        all_animal_corr = [all_animal_corr zscore(animal_corr(1:7))];
        
fig = figure;
    hold on;
    scatter(all_animal_coh,all_animal_corr,30,[0 0 0],'filled')
    [R,P] = corrcoef(all_animal_coh,all_animal_corr);
    [p] = polyfit(all_animal_coh,all_animal_corr,1);
    x=[-2.5:.1:2.5];
    xlim([-2.5 2.5]);
    y = x*p(1)+p(2);
    plot(x,y,'color','r','LineStyle','--','LineWidth',3)
    title(['BEHAVIOR VS LFP COH- R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);

%% FIG 3 D | HISTOGRAM OF MEAN ONLINE AND OFFLINE 4-8HZ LFP COHERENCE CHANGE ACROSS ALL M1-DLS ELECTRODE PAIRS

    all_animal_day_change = [];
    all_animal_night_change = [];
    all_increase_chan = [];
    all_animal_mean_day_change = [];
    all_animal_mean_night_change = [];

    %%% LC3

        load([data_path '\LC3\sleep_activity\LFP_coherence.mat'])

        increase_channels = [];
        for pairs=1:256
            if mean(mean([pre_LFP_coh{pairs}(8,48:113); post_LFP_coh{pairs}(8,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)])))+0.025
                increase_channels = [increase_channels 1];
            else
                increase_channels = [increase_channels 0];
            end
        end
        all_increase_chan = [all_increase_chan sum(increase_channels==1)];

        day_change = [];
        for pairs = find(increase_channels)
            tmp_day_change = [];
            for day=1:8
                tmp_day_change = [tmp_day_change mean(post_LFP_coh{pairs}(day,48:113))-mean(pre_LFP_coh{pairs}(day,48:113))];
            end
            day_change = [day_change mean(tmp_day_change)];
        end

        night_change = [];
        for pairs = find(increase_channels)
            tmp_night_change = [];
            for day=1:7
                tmp_night_change = [tmp_night_change mean(pre_LFP_coh{pairs}(day+1,48:113))-mean(post_LFP_coh{pairs}(day,48:113))];
            end
            night_change = [night_change mean(tmp_night_change)];
        end

        figure;
            subplot(1,2,1);
                histogram(day_change,[-.05:0.005:.05],'DisplayStyle','Stairs')
                [h p] = ttest(day_change);
                title(num2str(p))                
            subplot(1,2,2);        
                histogram(night_change,[-.05:0.005:.05],'DisplayStyle','Stairs')
                [h p] = ttest(night_change);
                title(num2str(p))                 
                
        all_animal_day_change = [all_animal_day_change day_change];
        all_animal_night_change = [all_animal_night_change night_change];
        all_animal_mean_day_change = [all_animal_mean_day_change mean(day_change)];
        all_animal_mean_night_change = [all_animal_mean_night_change mean(night_change)];
        
    %%% LC4

        load([data_path '\LC4\sleep_activity\LFP_coherence.mat'])

        increase_channels = [];
        for pairs=1:256
            if mean(mean([pre_LFP_coh{pairs}(6,48:113); post_LFP_coh{pairs}(6,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(2,48:113); post_LFP_coh{pairs}(2,48:113)])))+0.025
                increase_channels = [increase_channels 1];
            else
                increase_channels = [increase_channels 0];
            end
        end
        all_increase_chan = [all_increase_chan sum(increase_channels==1)];

        day_change = [];
        for pairs = find(increase_channels)
            tmp_day_change = [];
            for day=2:6
                tmp_day_change = [tmp_day_change mean(post_LFP_coh{pairs}(day,48:113))-mean(pre_LFP_coh{pairs}(day,48:113))];
            end
            day_change = [day_change mean(tmp_day_change)];
        end

        night_change = [];
        for pairs = find(increase_channels)
            tmp_night_change = [];
            for day=2:5
                tmp_night_change = [tmp_night_change mean(pre_LFP_coh{pairs}(day+1,48:113))-mean(post_LFP_coh{pairs}(day,48:113))];
            end
            night_change = [night_change mean(tmp_night_change)];
        end

        figure;
            subplot(1,2,1);
                histogram(day_change,[-.1:0.005:.1],'DisplayStyle','Stairs')
                [h p] = ttest(day_change);
                title(num2str(p)) 
            subplot(1,2,2);        
                histogram(night_change,[-.1:0.005:.1],'DisplayStyle','Stairs')
                [h p] = ttest(night_change);
                title(num2str(p)) 

        all_animal_day_change = [all_animal_day_change day_change];
        all_animal_night_change = [all_animal_night_change night_change];
        all_animal_mean_day_change = [all_animal_mean_day_change mean(day_change)];
        all_animal_mean_night_change = [all_animal_mean_night_change mean(night_change)];

    %%% LC6

        load([data_path '\LC6\sleep_activity\LFP_coherence.mat'])

        increase_channels = [];
        for pairs=1:256
            if mean(mean([pre_LFP_coh{pairs}(9,48:113); post_LFP_coh{pairs}(9,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)])))+0.025
                increase_channels = [increase_channels 1];
            else
                increase_channels = [increase_channels 0];
            end
        end
        all_increase_chan = [all_increase_chan sum(increase_channels==1)];

        day_change = [];
        for pairs = find(increase_channels)
            tmp_day_change = [];
            for day=1:9
                tmp_day_change = [tmp_day_change mean(post_LFP_coh{pairs}(day,48:113))-mean(pre_LFP_coh{pairs}(day,48:113))];
            end
            day_change = [day_change mean(tmp_day_change)];
        end

        night_change = [];
        for pairs = find(increase_channels)
            tmp_night_change = [];
            for day=1:8
                tmp_night_change = [tmp_night_change mean(pre_LFP_coh{pairs}(day+1,48:113))-mean(post_LFP_coh{pairs}(day,48:113))];
            end
            night_change = [night_change mean(tmp_night_change)];
        end

        figure;
            subplot(1,2,1);
                histogram(day_change,[-.05:0.005:.05],'DisplayStyle','Stairs')
                [h p] = ttest(day_change);
                title(num2str(p)) 
            subplot(1,2,2);        
                histogram(night_change,[-.05:0.005:.05],'DisplayStyle','Stairs')
                [h p] = ttest(night_change);
                title(num2str(p)) 
 
        all_animal_day_change = [all_animal_day_change day_change];
        all_animal_night_change = [all_animal_night_change night_change];
        all_animal_mean_day_change = [all_animal_mean_day_change mean(day_change)];
        all_animal_mean_night_change = [all_animal_mean_night_change mean(night_change)];

    %%% LC5

        load([data_path '\LC5\sleep_activity\LFP_coherence.mat'])

        increase_channels = [];

        for pairs=1:512
            if mean(mean([pre_LFP_coh{pairs}(15,48:113); post_LFP_coh{pairs}(15,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(2,48:113); post_LFP_coh{pairs}(2,48:113)])))+0.025
                increase_channels = [increase_channels 1];
            else
                increase_channels = [increase_channels 0];
            end
        end
        all_increase_chan = [all_increase_chan sum(increase_channels==1)];

        day_change = [];
        for pairs = find(increase_channels)
            tmp_day_change = [];
            for day=2:15
                tmp_day_change = [tmp_day_change mean(post_LFP_coh{pairs}(day,48:113))-mean(pre_LFP_coh{pairs}(day,48:113))];
            end
            day_change = [day_change mean(tmp_day_change)];
        end

        night_change = [];
        for pairs = find(increase_channels)
            tmp_night_change = [];
            for day=2:14
                tmp_night_change = [tmp_night_change mean(pre_LFP_coh{pairs}(day+1,48:113))-mean(post_LFP_coh{pairs}(day,48:113))];
            end
            night_change = [night_change mean(tmp_night_change)];
        end

        figure;
            subplot(1,2,1);
                histogram(day_change,[-.03:0.005:.03],'DisplayStyle','Stairs')
                [h p] = ttest(day_change);
                title(num2str(p)) 
            subplot(1,2,2);        
                histogram(night_change,[-.03:0.005:.03],'DisplayStyle','Stairs')
                [h p] = ttest(night_change);
                title(num2str(p)) 

        all_animal_day_change = [all_animal_day_change day_change];
        all_animal_night_change = [all_animal_night_change night_change];
        all_animal_mean_day_change = [all_animal_mean_day_change mean(day_change)];
        all_animal_mean_night_change = [all_animal_mean_night_change mean(night_change)];

    %%% LC2

        load([data_path '\LC2\sleep_activity\LFP_coherence.mat'])

        increase_channels = [];
        for pairs=1:1024
            if mean(mean([pre_LFP_coh{pairs}(10,40:95); post_LFP_coh{pairs}(10,40:95)])) > (mean(mean([pre_LFP_coh{pairs}(1,40:95); post_LFP_coh{pairs}(1,40:95)])))+0.025
                increase_channels = [increase_channels 1];
            else
                increase_channels = [increase_channels 0];
            end
        end
        all_increase_chan = [all_increase_chan sum(increase_channels==1)];

        day_change = [];
        for pairs = find(increase_channels)
            tmp_day_change = [];
            for day=1:10
                tmp_day_change = [tmp_day_change mean(post_LFP_coh{pairs}(day,40:95))-mean(pre_LFP_coh{pairs}(day,40:95))];
            end
            day_change = [day_change mean(tmp_day_change)];
        end

        night_change = [];
        for pairs = find(increase_channels)
            tmp_night_change = [];
            for day=1:9
                tmp_night_change = [tmp_night_change mean(pre_LFP_coh{pairs}(day+1,40:95))-mean(post_LFP_coh{pairs}(day,40:95))];
            end
            night_change = [night_change mean(tmp_night_change)];
        end

        figure;
            subplot(1,2,1);
                histogram(day_change,[-.1:0.005:.1],'DisplayStyle','Stairs')
                [h p] = ttest(day_change);
                title(num2str(p)) 
            subplot(1,2,2);        
                histogram(night_change,[-.1:0.005:.1],'DisplayStyle','Stairs')
                [h p] = ttest(night_change);
                title(num2str(p)) 

        all_animal_day_change = [all_animal_day_change day_change];
        all_animal_night_change = [all_animal_night_change night_change];
        all_animal_mean_day_change = [all_animal_mean_day_change mean(day_change)];
        all_animal_mean_night_change = [all_animal_mean_night_change mean(night_change)];
        
    %%% LC1
    
        load([data_path '\LC1\sleep_activity\LFP_coherence.mat'])

        increase_channels = [];
        for pairs=1:1406
            if mean(mean([pre_LFP_coh{pairs}(7,40:95); post_LFP_coh{pairs}(7,40:95)])) > (mean(mean([pre_LFP_coh{pairs}(1,40:95); post_LFP_coh{pairs}(1,40:95)])))+0.025
                increase_channels = [increase_channels 1];
            else
                increase_channels = [increase_channels 0];
            end
        end
        all_increase_chan = [all_increase_chan sum(increase_channels==1)];

        day_change = [];
        for pairs = find(increase_channels)
            tmp_day_change = [];
            for day=1:7
                tmp_day_change = [tmp_day_change mean(post_LFP_coh{pairs}(day,40:95))-mean(pre_LFP_coh{pairs}(day,40:95))];
            end
            day_change = [day_change mean(tmp_day_change)];
        end

        night_change = [];
        for pairs = find(increase_channels)
            tmp_night_change = [];
            for day=1:6
                tmp_night_change = [tmp_night_change mean(pre_LFP_coh{pairs}(day+1,40:95))-mean(post_LFP_coh{pairs}(day,40:95))];
            end
            night_change = [night_change mean(tmp_night_change)];
        end

        figure;
            subplot(1,2,1);
                histogram(day_change,[-.1:0.005:.1],'DisplayStyle','Stairs')
                [h p] = ttest(day_change);
                title(num2str(p)) 
            subplot(1,2,2);        
                histogram(night_change,[-.1:0.005:.1],'DisplayStyle','Stairs')
                [h p] = ttest(night_change);
                title(num2str(p)) 

        all_animal_day_change = [all_animal_day_change day_change];
        all_animal_night_change = [all_animal_night_change night_change];
        all_animal_mean_day_change = [all_animal_mean_day_change mean(day_change)];
        all_animal_mean_night_change = [all_animal_mean_night_change mean(night_change)];
        
    %%% PLOT ALL

        fig = figure;
        subplot(1,2,1); hold on;
        for n = 1:length(all_animal_mean_day_change)
            plot([0 all_animal_mean_day_change(n)], [n n],'k')
        end
        xlim([-0.03 0.03])
        subplot(1,2,2); hold on;
        for n = 1:length(all_animal_mean_day_change)
            plot([0 all_animal_mean_night_change(n)], [n n],'k')
        end
        xlim([-0.03 0.03])

        fig = figure;
            subplot(1,2,1)
                hold on;
                histogram(all_animal_day_change,[-0.1:0.005:0.1],'DisplayStyle','Stairs','normalization','probability')
                ylim([0 0.2])
                plot([0 0],[0 0.2])
                plot([mean(all_animal_day_change) mean(all_animal_day_change)],[0 0.5],'color','r')
                xlim([-.1 .1])
                [h p ci stats] = ttest(all_animal_day_change)
                title(num2str(p))
            subplot(1,2,2)
                hold on;
                histogram(all_animal_night_change,[-0.1:0.005:0.1],'DisplayStyle','Stairs','normalization','probability')
                ylim([0 0.2])
                plot([0 0],[0 0.2])
                plot([mean(all_animal_night_change) mean(all_animal_night_change)],[0 0.5],'color','r')
                xlim([-.1 .1])
                [h p ci stats] = ttest(all_animal_night_change)
                title(num2str(p))

%% CALCULATE PARTIAL CORRELATION BETWEEEN LFP COHERENCE, VEL PROF CORR, AND MAX VEL

all_animal_coh = [];
all_animal_corr = [];
all_animal_vel = [];
        
%%% LC3

    %%% BEHAVIOR

        xvel = cell(1,8);
        yvel = cell(1,8);
        day_velocity = cell(1,8);
        
        day_count = 1;

        for day = 1:8

            files = dir([data_path '\LC3\reach_trajectories\day_' num2str(day) '\*.csv']);

            %%% find mean pellet position across trials
            pellet_position = [];
            for file = 1:length(files)
                trial_kin = readmatrix([data_path '\LC3\reach_trajectories\day_' num2str(day) '\' files(file).name]);
                tmp_pellet_x = trial_kin(1:25,5);
                tmp_pellet_y = trial_kin(1:25,6);
                tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
                pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                pellet_position = [pellet_position; pellet_location];
            end
            pellet_position = nanmean(pellet_position);

            %%% get velocity profiles time locked to frame with paw closest to pellet
            all_diffx = [];
            all_diffy = [];

            for file = 1:length(files)

                trial_kin = readmatrix([data_path '\LC3\reach_trajectories\day_' num2str(day) '\' files(file).name]);

                %%% get paw trajectory
                tmp_x = trial_kin(:,2);
                tmp_y = trial_kin(:,3);
                tmp_x(trial_kin(:,4)<.99) = NaN;
                tmp_y(trial_kin(:,4)<.99) = NaN;

                %%% get paw velocity
                diff_x = diff(tmp_x);
                diff_y = diff(tmp_y);

                %%% find frame with paw closest to pellet
                [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));

                %%% skip bad trials
                if ~(i<31 || i>size(trial_kin,1)-31)
                    all_diffx = [all_diffx; diff_x(i-30:i+30)'];
                    all_diffy = [all_diffy; diff_y(i-30:i+30)'];
                end

            end
            
            day_velocity{day_count} = (max(((abs(all_diffx)+abs(all_diffy))/2)')*0.0711)*30;
            
            xvel{day_count} = all_diffx;
            yvel{day_count} = all_diffy;
            day_count = day_count + 1;

        end

        %%% calculate correlation between each day's trial-averaged velocity profile and trial-averaged velocity profile on last day
        x_corr = [];
        y_corr = [];
        animal_vel = [];
        for n = 1:length(xvel)
            x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
            y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
            animal_vel = [animal_vel nanmean(day_velocity{n})];
        end
        animal_corr = mean([x_corr; y_corr]);

    %%% SLEEP COHERENCE

        load([data_path '\LC3\sleep_activity\LFP_coherence.mat'])

        increase_channels = [];
        for pairs=1:256
            if mean(mean([pre_LFP_coh{pairs}(8,48:113); post_LFP_coh{pairs}(8,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)])))+0.025
                increase_channels = [increase_channels 1];
            else
                increase_channels = [increase_channels 0];
            end
        end

        day_pre = [];
        day_post = [];
        for pairs = find(increase_channels)
            tmp_pre = [];
            tmp_post = [];
            for day=[1:8]
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,48:113))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,48:113))];
            end
            day_pre = [day_pre; tmp_pre];
            day_post = [day_post; tmp_post];
        end

    %%% PLOT

        figure
        for day=1:8
            
            subplot(3,1,1); hold on;
                errorbar(day-0.1,mean(day_pre(:,day)),std(day_pre(:,day))/sqrt(length(day_pre(:,day))),'color',[0.5 0.5 0.5])
                errorbar(day+0.1,mean(day_post(:,day)),std(day_post(:,day))/sqrt(length(day_post(:,day))),'color',[0 0 0])
                plot([day-0.1 day+0.1],[mean(day_pre(:,day)) mean(day_post(:,day))],'color',[0.5 0.5 0.5])
                if day<8
                    plot([day+0.1 day+0.9],[mean(day_post(:,day)) mean(day_pre(:,day+1))],'color',[0 0 0])
                end
                xlim([0 9])
                
            subplot(3,1,2); hold on;
                scatter(day,animal_corr(day),50,[1 0 0],'filled')
                xlim([0 9])
                
            subplot(3,1,3); hold on;
                scatter(day,animal_vel(day),50,[1 0 0],'filled')
                xlim([0 9])
        end

        all_animal_coh = [all_animal_coh zscore(mean([day_pre]))];
        all_animal_corr = [all_animal_corr zscore(animal_corr)];
        all_animal_vel = [all_animal_vel zscore(animal_vel)];
        
%%% LC4

    %%% BEHAVIOR

        xvel = cell(1,5);
        yvel = cell(1,5);
        day_velocity = cell(1,5);
        day_count = 1;

        for day = 2:6

            files = dir([data_path '\LC4\reach_trajectories\day_' num2str(day) '\*.csv']);

            %%% find mean pellet position across trials
            pellet_position = [];
            for file = 1:length(files)
                trial_kin = readmatrix([data_path '\LC4\reach_trajectories\day_' num2str(day) '\' files(file).name]);
                tmp_pellet_x = trial_kin(1:25,5);
                tmp_pellet_y = trial_kin(1:25,6);
                tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
                pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                pellet_position = [pellet_position; pellet_location];
            end
            pellet_position = nanmean(pellet_position);

            %%% get velocity profiles time locked to frame with paw closest to pellet
            all_diffx = [];
            all_diffy = [];

            for file = 1:length(files)

                trial_kin = readmatrix([data_path '\LC4\reach_trajectories\day_' num2str(day) '\' files(file).name]);

                %%% get paw trajectory
                tmp_x = trial_kin(:,2);
                tmp_y = trial_kin(:,3);
                tmp_x(trial_kin(:,4)<.99) = NaN;
                tmp_y(trial_kin(:,4)<.99) = NaN;

                %%% get paw velocity
                diff_x = diff(tmp_x);
                diff_y = diff(tmp_y);

                %%% find frame with paw closest to pellet
                [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));

                %%% skip bad trials
                if ~(i<31 || i>size(trial_kin,1)-31)
                    all_diffx = [all_diffx; diff_x(i-30:i+30)'];
                    all_diffy = [all_diffy; diff_y(i-30:i+30)'];
                end

            end
            
            day_velocity{day_count} = (max(((abs(all_diffx)+abs(all_diffy))/2)')*0.0711)*30;
            
            xvel{day_count} = all_diffx;
            yvel{day_count} = all_diffy;
            day_count = day_count + 1;

        end

        %%% calculate correlation between each day's trial-averaged velocity profile and trial-averaged velocity profile on last day
        x_corr = [];
        y_corr = [];
        animal_vel = [];
        for n = 1:length(xvel)
            x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
            y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
            animal_vel = [animal_vel nanmean(day_velocity{n})];
        end
        animal_corr = mean([x_corr; y_corr]);

    %%% SLEEP COHERENCE

        load([data_path '\LC4\sleep_activity\LFP_coherence.mat'])

        increase_channels = [];
        for pairs=1:256
            if mean(mean([pre_LFP_coh{pairs}(6,48:113); post_LFP_coh{pairs}(6,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(2,48:113); post_LFP_coh{pairs}(2,48:113)])))+0.025
                increase_channels = [increase_channels 1];
            else
                increase_channels = [increase_channels 0];
            end
        end

        day_pre = [];
        day_post = [];
        for pairs = find(increase_channels)
            tmp_pre = [];
            tmp_post = [];
            for day=[2:6]
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,48:113))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,48:113))];
            end
            day_pre = [day_pre; tmp_pre];
            day_post = [day_post; tmp_post];
        end

    %%% PLOT

        figure
        for day=1:5
            
            subplot(3,1,1); hold on;
                errorbar(day-0.1,mean(day_pre(:,day)),std(day_pre(:,day))/sqrt(length(day_pre(:,day))),'color',[0.5 0.5 0.5])
                errorbar(day+0.1,mean(day_post(:,day)),std(day_post(:,day))/sqrt(length(day_post(:,day))),'color',[0 0 0])
                plot([day-0.1 day+0.1],[mean(day_pre(:,day)) mean(day_post(:,day))],'color',[0.5 0.5 0.5])
                if day<5
                    plot([day+0.1 day+0.9],[mean(day_post(:,day)) mean(day_pre(:,day+1))],'color',[0 0 0])
                end
                xlim([0 6])
                
            subplot(3,1,2); hold on;
                scatter(day,animal_corr(day),50,[1 0 0],'filled')
                xlim([0 6])
                
            subplot(3,1,3); hold on;
                scatter(day,animal_vel(day),50,[1 0 0],'filled')
                xlim([0 6])
        end

        all_animal_coh = [all_animal_coh zscore(mean([day_pre]))];
        all_animal_corr = [all_animal_corr zscore(animal_corr)];
        all_animal_vel = [all_animal_vel zscore(animal_vel)];
        
%%% LC6

    %%% BEHAVIOR

        xvel = cell(1,9);
        yvel = cell(1,9);
        reach_blocks = {[],[7],[12 13],[17 18],[21 22],[26],[29],[33 34],[38 39],[43 44]};
        day_velocity = cell(1,9);
        day_count = 1;

        for day = 2:10
    
            files = dir([data_path '\LC6\reach_trajectories\day_' num2str(day) '\*.csv']);
            files = sort_nat({files.name});
            vid_timestamp = dir([data_path '\LC6\reach_trajectories\day_' num2str(day) '\*.mat']); 
            vid_timestamp = sort_nat({vid_timestamp.name});

            all_diffx = [];
            all_diffy = [];

            for block = 1:length(files)

                %%% load all trial kinematics
                kin_all = readmatrix([data_path '\LC6\reach_trajectories\day_' num2str(day) '\' files{block}]);
                load([data_path '\LC6\reach_trajectories\day_' num2str(day) '\' vid_timestamp{block}]);

                %%% find individual trials 
                [~,i] = find(Vid.scalars.Vid0.data>1e8);
                for n_i = 1:length(i)
                    Vid.scalars.Vid0.data(i(n_i)) = round(mean(Vid.scalars.Vid0.data(i(n_i)-1):Vid.scalars.Vid0.data(i(n_i)+1)));
                end
                Vid.scalars.Vid0.data = Vid.scalars.Vid0.data+1;
                first_frames = find(diff(Vid.scalars.Vid0.ts)>8);
                trial_starts_idx = [Vid.scalars.Vid0.data(1) Vid.scalars.Vid0.data(first_frames+1)];
                trial_ends_idx = [Vid.scalars.Vid0.data(first_frames) Vid.scalars.Vid0.data(end)];

                %%% find mean pellet position across trials
                pellet_position = [];
                for n_trial = 1:length(trial_starts_idx)-1

                    %%% get timesstamps (ts) for frames in trial
                    trial_frames_idx = Vid.scalars.Vid0.data(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));
                    trial_frames_ts = Vid.scalars.Vid0.ts(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));

                    %%% interpolate to get missing ts
                    interp_trial_ts = interp1(trial_frames_idx,trial_frames_ts,trial_starts_idx(n_trial):trial_ends_idx(n_trial));

                    %%% get trial kinematics
                    trial_kin = kin_all(trial_starts_idx(n_trial):trial_ends_idx(n_trial),:);

                    %%% get pellet location from first 25 frames
                    tmp_pellet_x = trial_kin(1:25,5);
                    tmp_pellet_y = trial_kin(1:25,6);
                    tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
                    tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
                    pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                    pellet_position = [pellet_position; pellet_location];

                end
                pellet_position = nanmean(pellet_position);

                for n_trial = 1:length(trial_starts_idx)-1

                    %%% get timesstamps (ts) for frames in trial
                    trial_frames_idx = Vid.scalars.Vid0.data(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));
                    trial_frames_ts = Vid.scalars.Vid0.ts(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));

                    %%% interpolate to get missing ts
                    interp_trial_ts = interp1(trial_frames_idx,trial_frames_ts,trial_starts_idx(n_trial):trial_ends_idx(n_trial));

                    %%% get trial kinematics
                    trial_kin = kin_all(trial_starts_idx(n_trial):trial_ends_idx(n_trial),:);

                    %%% get paw trajectory
                    tmp_x = trial_kin(:,2);
                    tmp_y = trial_kin(:,3);
                    tmp_x(trial_kin(:,4)<.99) = NaN;
                    tmp_y(trial_kin(:,4)<.99) = NaN;

                    %%% get paw velocity
                    diff_x = diff(tmp_x);
                    diff_y = diff(tmp_y);

                    %%% find frame with paw closest to pellet
                    [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));

                    %%% skip bad trials
                    if ~(i<31 || i>size(trial_kin,1)-31)
                        all_diffx = [all_diffx; diff_x(i-30:i+30)'];
                        all_diffy = [all_diffy; diff_y(i-30:i+30)'];
                    end

                end
            end
            
            day_velocity{day_count} = (max(((abs(all_diffx)+abs(all_diffy))/2)')*0.0222)*60;

            xvel{day_count} = all_diffx;
            yvel{day_count} = all_diffy;
            day_count = day_count + 1;

        end

        %%% calculate correlation between each day's trial-averaged velocity profile and trial-averaged velocity profile on last day
        x_corr = [];
        y_corr = [];
        animal_vel = [];
        for n = 1:length(xvel)
            x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
            y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
            animal_vel = [animal_vel nanmean(day_velocity{n})];
        end
        animal_corr = mean([x_corr; y_corr]);

    %%% SLEEP COHERENCE

        load([data_path '\LC6\sleep_activity\LFP_coherence.mat'])

        increase_channels = [];
        for pairs=1:256
            if mean(mean([pre_LFP_coh{pairs}(9,48:113); post_LFP_coh{pairs}(9,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)])))+0.025
                increase_channels = [increase_channels 1];
            else
                increase_channels = [increase_channels 0];
            end
        end

        day_pre = [];
        day_post = [];
        for pairs = find(increase_channels)
            tmp_pre = [];
            tmp_post = [];
            for day=[1:9]
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,48:113))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,48:113))];
            end
            day_pre = [day_pre; tmp_pre];
            day_post = [day_post; tmp_post];
        end

    %%% PLOT

        figure
        for day=1:9
            
            subplot(3,1,1); hold on;
                errorbar(day-0.1,mean(day_pre(:,day)),std(day_pre(:,day))/sqrt(length(day_pre(:,day))),'color',[0.5 0.5 0.5])
                errorbar(day+0.1,mean(day_post(:,day)),std(day_post(:,day))/sqrt(length(day_post(:,day))),'color',[0 0 0])
                plot([day-0.1 day+0.1],[mean(day_pre(:,day)) mean(day_post(:,day))],'color',[0.5 0.5 0.5])
                if day<9
                    plot([day+0.1 day+0.9],[mean(day_post(:,day)) mean(day_pre(:,day+1))],'color',[0 0 0])
                end
                xlim([0 10])
                
            subplot(3,1,2); hold on;
                scatter(day,animal_corr(day),50,[1 0 0],'filled')
                xlim([0 10])
            subplot(3,1,3); hold on;
                scatter(day,animal_vel(day),50,[1 0 0],'filled')
                xlim([0 10])
        end

        all_animal_coh = [all_animal_coh zscore(mean([day_pre]))];
        all_animal_corr = [all_animal_corr zscore(animal_corr)];
        all_animal_vel = [all_animal_vel zscore(animal_vel)];
        
%%% LC5

    %%% BEHAVIOR

        xvel = cell(1,14);
        yvel = cell(1,14);
        reach_blocks = {[],[],[],[9],[12],[16],[19],[22 23 24],[28],[31],[35],[39],[42],[46],[50],[54],[59]};
        day_velocity = cell(1,14);
        day_count = 1;
 
        for day = 4:17
    
            files = dir([data_path '\LC5\reach_trajectories\day_' num2str(day) '\*.csv']);
            files = sort_nat({files.name});
            vid_timestamp = dir([data_path '\LC5\reach_trajectories\day_' num2str(day) '\*.mat']); 
            vid_timestamp = sort_nat({vid_timestamp.name});

            all_diffx = [];
            all_diffy = [];

            for block = 1:length(files)

                %%% load all trial kinematics
                kin_all = readmatrix([data_path '\LC5\reach_trajectories\day_' num2str(day) '\' files{block}]);
                load([data_path '\LC5\reach_trajectories\day_' num2str(day) '\' vid_timestamp{block}]);

                %%% find individual trials 
                [~,i] = find(Vid.scalars.Vid0.data>1e8);
                for n_i = 1:length(i)
                    Vid.scalars.Vid0.data(i(n_i)) = round(mean(Vid.scalars.Vid0.data(i(n_i)-1):Vid.scalars.Vid0.data(i(n_i)+1)));
                end
                Vid.scalars.Vid0.data = Vid.scalars.Vid0.data+1;
                first_frames = find(diff(Vid.scalars.Vid0.ts)>8);
                trial_starts_idx = [Vid.scalars.Vid0.data(1) Vid.scalars.Vid0.data(first_frames+1)];
                trial_ends_idx = [Vid.scalars.Vid0.data(first_frames) Vid.scalars.Vid0.data(end)];

                %%% find mean pellet position across trials
                pellet_position = [];
                for n_trial = 1:length(trial_starts_idx)-1

                    %%% get timesstamps (ts) for frames in trial
                    trial_frames_idx = Vid.scalars.Vid0.data(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));
                    trial_frames_ts = Vid.scalars.Vid0.ts(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));

                    %%% interpolate to get missing ts
                    interp_trial_ts = interp1(trial_frames_idx,trial_frames_ts,trial_starts_idx(n_trial):trial_ends_idx(n_trial));

                    %%% get trial kinematics
                    trial_kin = kin_all(trial_starts_idx(n_trial):trial_ends_idx(n_trial),:);

                    %%% get pellet location from first 25 frames
                    tmp_pellet_x = trial_kin(1:25,5);
                    tmp_pellet_y = trial_kin(1:25,6);
                    tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
                    tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
                    pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                    pellet_position = [pellet_position; pellet_location];

                end
                pellet_position = nanmean(pellet_position);

                for n_trial = 1:length(trial_starts_idx)-1

                    %%% get timesstamps (ts) for frames in trial
                    trial_frames_idx = Vid.scalars.Vid0.data(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));
                    trial_frames_ts = Vid.scalars.Vid0.ts(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));

                    %%% interpolate to get missing ts
                    interp_trial_ts = interp1(trial_frames_idx,trial_frames_ts,trial_starts_idx(n_trial):trial_ends_idx(n_trial));

                    %%% get trial kinematics
                    trial_kin = kin_all(trial_starts_idx(n_trial):trial_ends_idx(n_trial),:);

                    %%% get paw trajectory
                    tmp_x = trial_kin(:,2);
                    tmp_y = trial_kin(:,3);
                    tmp_x(trial_kin(:,4)<.99) = NaN;
                    tmp_y(trial_kin(:,4)<.99) = NaN;

                    %%% get paw velocity
                    diff_x = diff(tmp_x);
                    diff_y = diff(tmp_y);

                    %%% find frame with paw closest to pellet
                    [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));

                    %%% skip bad trials
                    if ~(i<31 || i>size(trial_kin,1)-31)
                        all_diffx = [all_diffx; diff_x(i-30:i+30)'];
                        all_diffy = [all_diffy; diff_y(i-30:i+30)'];
                    end

                end
            end
            
            day_velocity{day_count} = (max(((abs(all_diffx)+abs(all_diffy))/2)')*0.0254)*60;
            
            xvel{day_count} = all_diffx;
            yvel{day_count} = all_diffy;
            day_count = day_count + 1;

        end

        %%% calculate correlation between each day's trial-averaged velocity profile and trial-averaged velocity profile on last day
        x_corr = [];
        y_corr = [];
        animal_vel = [];
        for n = 1:length(xvel)
            x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
            y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
            animal_vel = [animal_vel nanmean(day_velocity{n})];
        end
        animal_corr = mean([x_corr; y_corr]);

    %%% SLEEP COHERENCE

        load([data_path '\LC5\sleep_activity\LFP_coherence.mat'])

        increase_channels = [];
        for pairs=1:512
            if mean(mean([pre_LFP_coh{pairs}(15,48:113); post_LFP_coh{pairs}(15,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(2,48:113); post_LFP_coh{pairs}(2,48:113)])))+0.025
                increase_channels = [increase_channels 1];
            else
                increase_channels = [increase_channels 0];
            end
        end

        day_pre = [];
        day_post = [];
        for pairs = find(increase_channels)
            tmp_pre = [];
            tmp_post = [];
            for day=[2:15]
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,48:113))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,48:113))];
            end
            day_pre = [day_pre; tmp_pre];
            day_post = [day_post; tmp_post];
        end

    %%% PLOT

        figure
        for day=1:14
            
            subplot(3,1,1); hold on;
                errorbar(day-0.1,mean(day_pre(:,day)),std(day_pre(:,day))/sqrt(length(day_pre(:,day))),'color',[0.5 0.5 0.5])
                errorbar(day+0.1,mean(day_post(:,day)),std(day_post(:,day))/sqrt(length(day_post(:,day))),'color',[0 0 0])
                plot([day-0.1 day+0.1],[mean(day_pre(:,day)) mean(day_post(:,day))],'color',[0.5 0.5 0.5])
                if day<14
                    plot([day+0.1 day+0.9],[mean(day_post(:,day)) mean(day_pre(:,day+1))],'color',[0 0 0])
                end
                xlim([0 15])
                
            subplot(3,1,2); hold on;
                scatter(day,animal_corr(day),50,[1 0 0],'filled')
                xlim([0 15])
            subplot(3,1,3); hold on;
                scatter(day,animal_vel(day),50,[1 0 0],'filled')
                xlim([0 15])
        end

        all_animal_coh = [all_animal_coh zscore(mean([day_pre]))];
        all_animal_corr = [all_animal_corr zscore(animal_corr)];
        all_animal_vel = [all_animal_vel zscore(animal_vel)];
        
%%% LC2

    %%% BEHAVIOR

        xvel = cell(1,10);
        yvel = cell(1,10);
        day_velocity = cell(1,10);
        day_count = 1;
        
        for day = 4:13

            files = dir([data_path '\LC2\reach_trajectories\day_' num2str(day) '\*.csv']);

            %%% find mean pellet position across trials
            pellet_position = [];
            for file = 1:length(files)
                trial_kin = readmatrix([data_path '\LC2\reach_trajectories\day_' num2str(day) '\' files(file).name]);
                tmp_pellet_x = trial_kin(1:25,8);
                tmp_pellet_y = trial_kin(1:25,9);
                tmp_pellet_x(trial_kin(1:25,10)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,10)<.99) = NaN;
                pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                pellet_position = [pellet_position; pellet_location];
            end
            pellet_position = nanmean(pellet_position);

            %%% get velocity profiles time locked to frame with paw closest to pellet
            all_diffx = [];
            all_diffy = [];

            for file = 1:length(files)

                trial_kin = readmatrix([data_path '\LC2\reach_trajectories\day_' num2str(day) '\' files(file).name]);

                %%% get paw trajectory
                tmp_x = trial_kin(:,2);
                tmp_y = trial_kin(:,3);
                tmp_x(trial_kin(:,4)<.99) = NaN;
                tmp_y(trial_kin(:,4)<.99) = NaN;

                %%% get paw velocity
                diff_x = diff(tmp_x);
                diff_y = diff(tmp_y);

                %%% find frame with paw closest to pellet
                [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));

                %%% skip bad trials
                if ~(i<31 || i>size(trial_kin,1)-31)
                    all_diffx = [all_diffx; diff_x(i-30:i+30)'];
                    all_diffy = [all_diffy; diff_y(i-30:i+30)'];
                end

            end

            day_velocity{day_count} = (max(((abs(all_diffx)+abs(all_diffy))/2)')*0.0237)*75;
            
            xvel{day_count} = all_diffx;
            yvel{day_count} = all_diffy;
            day_count = day_count + 1;

        end

        %%% calculate correlation between each day's trial-averaged velocity profile and trial-averaged velocity profile on last day
        x_corr = [];
        y_corr = [];
        animal_vel = [];
        for n = 1:length(xvel)
            x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
            y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
            animal_vel = [animal_vel nanmean(day_velocity{n})];
        end
        animal_corr = mean([x_corr; y_corr]);

    %%% SLEEP COHERENCE

        load([data_path '\LC2\sleep_activity\LFP_coherence.mat'])

        increase_channels = [];
        for pairs=1:1024
            if mean(mean([pre_LFP_coh{pairs}(10,40:95); post_LFP_coh{pairs}(10,40:95)])) > (mean(mean([pre_LFP_coh{pairs}(1,40:95); post_LFP_coh{pairs}(1,40:95)])))+0.025
                increase_channels = [increase_channels 1];
            else
                increase_channels = [increase_channels 0];
            end
        end

        day_pre = [];
        day_post = [];
        for pairs = find(increase_channels)
            tmp_pre = [];
            tmp_post = [];
            for day=1:10
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,40:95))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,40:95))];
            end
            day_pre = [day_pre; tmp_pre];
            day_post = [day_post; tmp_post];
        end

    %%% PLOT

        figure
        for day=1:10
            
            subplot(3,1,1); hold on;
                errorbar(day-0.1,mean(day_pre(:,day)),std(day_pre(:,day))/sqrt(length(day_pre(:,day))),'color',[0.5 0.5 0.5])
                errorbar(day+0.1,mean(day_post(:,day)),std(day_post(:,day))/sqrt(length(day_post(:,day))),'color',[0 0 0])
                plot([day-0.1 day+0.1],[mean(day_pre(:,day)) mean(day_post(:,day))],'color',[0.5 0.5 0.5])
                if day<10
                    plot([day+0.1 day+0.9],[mean(day_post(:,day)) mean(day_pre(:,day+1))],'color',[0 0 0])
                end
                xlim([0 11])
                
            subplot(3,1,2); hold on;
                scatter(day,animal_corr(day),50,[1 0 0],'filled')
                xlim([0 11])
            subplot(3,1,3); hold on;
                scatter(day,animal_vel(day),50,[1 0 0],'filled')
                xlim([0 11])
        end

        all_animal_coh = [all_animal_coh zscore(mean([day_pre]))];
        all_animal_corr = [all_animal_corr zscore(animal_corr)];
        all_animal_vel = [all_animal_vel zscore(animal_vel)];
        
%%% LC1

    %%% BEHAVIOR

        xvel = cell(1,8);
        yvel = cell(1,8);
        day_velocity = cell(1,8);
        day_count = 1;
        
        for day = [1:8]

            files = dir([data_path '\LC1\reach_trajectories\day_' num2str(day) '\*.csv']);

            %%% find mean pellet position across trials
            pellet_position = [];
            for file = 1:length(files)
                trial_kin = readmatrix([data_path '\LC1\reach_trajectories\day_' num2str(day) '\' files(file).name]);
                tmp_pellet_x = trial_kin(1:25,8);
                tmp_pellet_y = trial_kin(1:25,9);
                tmp_pellet_x(trial_kin(1:25,10)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,10)<.99) = NaN;
                pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                pellet_position = [pellet_position; pellet_location];
            end
            pellet_position = nanmean(pellet_position);

            %%% get velocity profiles time locked to frame with paw closest to pellet
            all_diffx = [];
            all_diffy = [];

            for file = 1:length(files)

                trial_kin = readmatrix([data_path '\LC1\reach_trajectories\day_' num2str(day) '\' files(file).name]);

                %%% get paw trajectory
                tmp_x = trial_kin(:,2);
                tmp_y = trial_kin(:,3);
                tmp_x(trial_kin(:,4)<.99) = NaN;
                tmp_y(trial_kin(:,4)<.99) = NaN;

                %%% get paw velocity
                diff_x = diff(tmp_x);
                diff_y = diff(tmp_y);

                %%% find frame with paw closest to pellet
                [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));

                %%% skip bad trials
                if ~(i<31 || i>size(trial_kin,1)-31)
                    all_diffx = [all_diffx; diff_x(i-30:i+30)'];
                    all_diffy = [all_diffy; diff_y(i-30:i+30)'];
                end

            end
            
            day_velocity{day_count} = (max(((abs(all_diffx)+abs(all_diffy))/2)')*0.0127)*100;

            xvel{day_count} = all_diffx;
            yvel{day_count} = all_diffy;
            day_count = day_count + 1;

        end

        %%% calculate correlation between each day's trial-averaged velocity profile and trial-averaged velocity profile on last day
        x_corr = [];
        y_corr = [];
        animal_vel = [];
        for n = 1:length(xvel)
            x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
            y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
            animal_vel = [animal_vel nanmean(day_velocity{n})];
        end
        animal_corr = mean([x_corr; y_corr]);

    %%% SLEEP COHERENCE

        load([data_path '\LC1\sleep_activity\LFP_coherence.mat'])

        increase_channels = [];
        for pairs=1:1406
            if mean(mean([pre_LFP_coh{pairs}(7,40:95); post_LFP_coh{pairs}(7,40:95)])) > (mean(mean([pre_LFP_coh{pairs}(1,40:95); post_LFP_coh{pairs}(1,40:95)])))+0.025
                increase_channels = [increase_channels 1];
            else
                increase_channels = [increase_channels 0];
            end
        end

        day_pre = [];
        day_post = [];
        for pairs = find(increase_channels)
            tmp_pre = [];
            tmp_post = [];
            for day=[1:7]
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,40:95))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,40:95))];
            end
            day_pre = [day_pre; tmp_pre];
            day_post = [day_post; tmp_post];
        end

    %%% PLOT

        figure
        for day=1:8
            
            if day<8
            subplot(3,1,1); hold on;
                errorbar(day-0.1,mean(day_pre(:,day)),std(day_pre(:,day))/sqrt(length(day_pre(:,day))),'color',[0.5 0.5 0.5])
                errorbar(day+0.1,mean(day_post(:,day)),std(day_post(:,day))/sqrt(length(day_post(:,day))),'color',[0 0 0])
                plot([day-0.1 day+0.1],[mean(day_pre(:,day)) mean(day_post(:,day))],'color',[0.5 0.5 0.5])
                if day<7
                    plot([day+0.1 day+0.9],[mean(day_post(:,day)) mean(day_pre(:,day+1))],'color',[0 0 0])
                end
                xlim([0 9])
            end
                
            subplot(3,1,2); hold on;
                scatter(day,animal_corr(day),50,[1 0 0],'filled')
                xlim([0 9])
            subplot(3,1,3); hold on;
                scatter(day,animal_vel(day),50,[1 0 0],'filled')
                xlim([0 9])
        end

        all_animal_coh = [all_animal_coh zscore(mean([day_pre]))];
        all_animal_corr = [all_animal_corr zscore(animal_corr(1:end-1))];
        all_animal_vel = [all_animal_vel zscore(animal_vel(1:end-1))];
        
        %%% CORRELATION AND PARTIAL CORRELATION
        [R,P] = corrcoef(all_animal_coh,all_animal_corr)
        [r p] = partialcorr([all_animal_coh' all_animal_corr'],all_animal_vel')

        [R,P] = corrcoef(all_animal_coh,all_animal_vel)
        [r p] = partialcorr([all_animal_coh' all_animal_vel'],all_animal_corr')

%% CALCULTE AMOUNT OF 4-8HZ LFP COHRENCE INCREASE CHANNEL PAIRS OVER LEARNING

increase_channels = [];

%%% LC3

    load([data_path '\LC3\sleep_activity\LFP_coherence.mat'])
    for pairs=1:256
        increase_channels = [increase_channels mean(mean([pre_LFP_coh{pairs}(8,48:113); post_LFP_coh{pairs}(8,48:113)]))-mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)]))];
    end
    
%%% LC4

    load([data_path '\LC4\sleep_activity\LFP_coherence.mat'])
    for pairs=1:256
        increase_channels = [increase_channels mean(mean([pre_LFP_coh{pairs}(6,48:113); post_LFP_coh{pairs}(6,48:113)]))-mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)]))];
    end

%%% LC6

    load([data_path '\LC6\sleep_activity\LFP_coherence.mat'])
    for pairs=1:256
        increase_channels = [increase_channels mean(mean([pre_LFP_coh{pairs}(9,48:113); post_LFP_coh{pairs}(9,48:113)]))-mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)]))];
    end

%%% LC5

    load([data_path '\LC5\sleep_activity\LFP_coherence.mat'])
    for pairs=1:512
        increase_channels = [increase_channels mean(mean([pre_LFP_coh{pairs}(15,48:113); post_LFP_coh{pairs}(15,48:113)]))-mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)]))];
    end

%%% LC2

    load([data_path '\LC2\sleep_activity\LFP_coherence.mat'])
    for pairs=1:1024
        increase_channels = [increase_channels mean(mean([pre_LFP_coh{pairs}(10,40:95); post_LFP_coh{pairs}(10,40:95)]))-mean(mean([pre_LFP_coh{pairs}(1,40:95); post_LFP_coh{pairs}(1,40:95)]))];
    end

%%% LC1

    load([data_path '\LC1\sleep_activity\LFP_coherence.mat'])
    for pairs=1:1406
        increase_channels = [increase_channels mean(mean([pre_LFP_coh{pairs}(8,40:95); post_LFP_coh{pairs}(8,40:95)]))-mean(mean([pre_LFP_coh{pairs}(1,40:95); post_LFP_coh{pairs}(1,40:95)]))];
    end

%%% PLOT

    sum(increase_channels>-0.025 & increase_channels<0.025)/length(increase_channels)
    sum(increase_channels>0.025)/length(increase_channels)
    sum(increase_channels<-0.025)/length(increase_channels)
    