%% SCRIPT TO GENERATE FIGURE 1 & FIGURE 1 FIGURE SUPPLEMENTS
% 
% to generate figures, download data from https://doi.org/10.7272/Q6KK9927

clc; clear all; close all;
data_path = 'D:\Lemke_etal_2021\Lemke_eLife_2021_data';

%% FIG 1 C & D | EXAMPLE REACH TRAJECTORY AND VELOCITY PROFILE ACROSS LEARNING 

    figure(1);
    figure(2);

    xvel = cell(1,8);
    yvel = cell(1,8);
    
    for day = 1:8

        figure(1);
        subplot(2,4,day); hold on;
        title(['day ' num2str(day)])
        
        %%% get DLC trajectory files
        files = dir([data_path '\LC1\reach_trajectories\day_' num2str(day) '\*LC1*']);

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

        %%% plot reach trajectories time locked to frame with paw closest to pellet 
        all_x = [];
        all_y = [];
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
            if ~(i<51 || i>size(trial_kin,1)-31)

                %%% plot individual trial reach trajectory                
                plot(tmp_x(i-30:i+30)'-pellet_position(1),tmp_y(i-30:i+30)'-pellet_position(2),'color',[0.5 0.5 0.5]);

                all_diffx = [all_diffx; diff_x(i-50:i+30)'];
                all_diffy = [all_diffy; diff_y(i-50:i+30)'];
                all_x = [all_x; tmp_x(i-30:i+30)'-pellet_position(1)];
                all_y = [all_y; tmp_y(i-30:i+30)'-pellet_position(2)];

            end

        end

        %%% plot mean reach trajectory
        plot(nanmean(all_x,1),nanmean(all_y,1),'color',[0+(day/8) 0 1-(day/8)],'LineWidth',1);
        xlim([-400 100])
        ylim([-110 175])
        set(gca, 'YDir','reverse')
        
        %%% plot mean velocity trajectory
        figure(2);
        subplot(1,2,1); hold on;
            plot(nanmean(all_diffx,1),'color',[0+(day/8) 0 1-(day/8)],'LineWidth',1);
            title('x velocity')
            xlim([1 80])
        subplot(1,2,2); hold on;
            plot(nanmean(all_diffy,1),'color',[0+(day/8) 0 1-(day/8)],'LineWidth',1);
            title('y velocity')
            xlim([1 80])
        
        xvel{day} = all_diffx;
        yvel{day} = all_diffy;    

    end

    %%% calculate correlation between each day's trial-averaged velocity profile and trial-averaged velocity profile on last day
    x_corr = [];
    y_corr = [];
    for n = 1:length(xvel)
        x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
        y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
    end
    figure(3)
    hold on
    plot(mean([x_corr; y_corr]),'k')
    for day = 1:8
        scatter(day,mean([x_corr(day) y_corr(day)]),50,[0 0 0],'filled')
    end
    ylim([0 1])
    xlabel('days')
    ylabel('correlation')
    title('velocity profile correlation')

%% FIG 1 D | LEARNING COHORT | VELOCITY PROFILE CORRELATION (ALL ANIMALS)

    learning_corr_day = cell(1,6);
    animalCount = 1;
    
    animals = {'LC3','LC4'};
    days = {1:8,2:6};
    for animal = 1:length(animals)

        animal_days = days{animal};
        xvel = cell(1,length(animal_days));
        yvel = cell(1,length(animal_days));

        day_count = 1;
        for day = animal_days

            files = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.csv']);

            %%% find mean pellet position across trials
            pellet_position = [];
            for file = 1:length(files)
                trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);
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

                trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);

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

        learning_corr_day{animalCount} = animal_corr;
        animalCount = animalCount + 1;
    end

    animals = {'LC2','LC1'};
    days = {4:13,1:8};
    for animal = 1:length(animals)

        animal_days = days{animal};
        xvel = cell(1,length(animal_days));
        yvel = cell(1,length(animal_days));

        day_count = 1;
        for day = animal_days

            files = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.csv']);

            %%% find mean pellet position across trials
            pellet_position = [];
            for file = 1:length(files)
                trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);
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

                trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);

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

        learning_corr_day{animalCount} = animal_corr;
        animalCount = animalCount + 1;
    end 

    animals = {'LC6','LC5'};
    days = {2:10,4:17};
    for animal = 1:length(animals)

        animal_days = days{animal};
        xvel = cell(1,length(animal_days));
        yvel = cell(1,length(animal_days));

        day_count = 1;
        for day = animal_days

            files = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.csv']);
            files = sort_nat({files.name});
            vid_timestamp = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.mat']); 
            vid_timestamp = sort_nat({vid_timestamp.name});

            all_diffx = [];
            all_diffy = [];

            for block = 1:length(files)

                %%% load all trial kinematics
                kin_all = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files{block}]);
                load([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' vid_timestamp{block}]);

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

        learning_corr_day{animalCount} = animal_corr;
        animalCount = animalCount + 1;
    end

    figure; hold on;
        for animal = 1:6
            day_value = [];
            for days = 1:8
                if days<numel(learning_corr_day{animal})
                    day_value = [day_value learning_corr_day{animal}(days)];
                end        
            end
            plot(day_value,'color','k')
             scatter([1:length(learning_corr_day{animal})-1],learning_corr_day{animal}(1:end-1),30,[0 0 0],'filled')
        end
        ylim([0 1])
        xlim([0.5 8.5])

%% FIG 1 E | sFIG 4 & 5 | AP5/SALINE COHORTS | VELOCITY PROFILE CORRELATION VS. SINGLE TRIAL PEAK VELOCITY 

%%% compute velocity profile correlation and single trial peak velocity 

    AP5_corr_day = [];
    saline_corr_day = [];
    AP5_corr_sum = [];
    saline_corr_sum = [];

    saline_vel_day = [];
    AP5_vel_day = [];
    saline_vel_sum = [];
    AP5_vel_sum = [];

    for animal = 1:6

        for day = 1:10

            files = dir([data_path '\AP' num2str(animal) '\behavioral_data\day_' num2str(day) '\*csv']);

            %%% FIND PELLET POSITION
            pellet_position = [];
            for file = 1:length(files)
                trial_kin = readmatrix([data_path '\AP' num2str(animal) '\behavioral_data\day_' num2str(day) '\' files(file).name]);
                tmp_pellet_x = trial_kin(1:25,5);
                tmp_pellet_y = trial_kin(1:25,6);
                tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
                pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                pellet_position = [pellet_position; pellet_location];
            end
            pellet_position = nanmean(pellet_position);

            %%% GET TRAJECTORIES AND VELOCITY PROFILES
            all_diffx = [];
            all_diffy = [];
            all_x = [];
            all_y = [];
            for file = 1:length(files)
                trial_kin = readmatrix([data_path '\AP' num2str(animal) '\behavioral_data\day_' num2str(day) '\' files(file).name]);
                %%% get paw trajectory and velocity
                tmp_x = trial_kin(:,2);
                tmp_y = trial_kin(:,3);
                tmp_x(trial_kin(:,4)<.99) = NaN;
                tmp_y(trial_kin(:,4)<.99) = NaN;
                diff_x = diff(tmp_x);
                diff_y = diff(tmp_y);
                %%% find frame with paw closest to pellet
                [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));
                %%% skip bad trials
                if ~(i<31 || i>size(trial_kin,1)-31)
                    all_diffx = [all_diffx; diff_x(i-30:i+30)'];
                    all_diffy = [all_diffy; diff_y(i-30:i+30)'];
                    all_x = [all_x; tmp_x(i-30:i+30)'-pellet_position(1)];
                    all_y = [all_y; tmp_y(i-30:i+30)'-pellet_position(2)];
                end
            end
            AP5_x{1,day} = all_x;
            AP5_y{1,day} = all_y;
            AP5_xvel{1,day} = all_diffx;
            AP5_yvel{1,day} = all_diffy;
            AP5_maxVel{1,day} = (max(((abs(all_diffx)+abs(all_diffy))/2)')*0.0237)*75;

        end

        AP5_x_corr = [];
        AP5_y_corr = [];
        for n = 1:length(AP5_xvel)
            AP5_x_corr = [AP5_x_corr corr(smooth(nanmean(AP5_xvel{1,n}(:,:)),1),smooth(nanmean(AP5_xvel{1,end}(:,:)),1))];
            AP5_y_corr = [AP5_y_corr corr(smooth(nanmean(AP5_yvel{1,n}(:,:)),1),smooth(nanmean(AP5_yvel{1,end}(:,:)),1))];
        end
        AP5_corr = mean([AP5_x_corr; AP5_y_corr]);

        tmp_velocity = [];
        for day = 1:length(AP5_maxVel)
            tmp_velocity = [tmp_velocity nanmean(AP5_maxVel{day})];
        end

        if (animal == 1 | animal==3 | animal==6)
            AP5_corr_day = [AP5_corr_day; AP5_corr(6:9) nan];
            saline_corr_day = [saline_corr_day; AP5_corr(1:5)];
            AP5_corr_sum = [AP5_corr_sum; sum(diff(AP5_corr(6:9)))];
            saline_corr_sum = [saline_corr_sum; sum(diff(AP5_corr(1:5)))];      

            AP5_vel_day = [AP5_vel_day; tmp_velocity(6:10)];
            saline_vel_day = [saline_vel_day; tmp_velocity(1:5)];      
            AP5_vel_sum = [AP5_vel_sum; sum(diff(tmp_velocity(6:10)))];
            saline_vel_sum = [saline_vel_sum; sum(diff(tmp_velocity(1:5)))];                 
        else
            saline_corr_day = [saline_corr_day; AP5_corr(6:9) nan];
            AP5_corr_day = [AP5_corr_day; AP5_corr(1:5)];
            saline_corr_sum = [saline_corr_sum; sum(diff(AP5_corr(6:9)))];
            AP5_corr_sum = [AP5_corr_sum; sum(diff(AP5_corr(1:5)))];      

            saline_vel_day = [saline_vel_day; tmp_velocity(6:10)];
            AP5_vel_day = [AP5_vel_day; tmp_velocity(1:5)];
            saline_vel_sum = [saline_vel_sum; sum(diff(tmp_velocity(6:10)))];
            AP5_vel_sum = [AP5_vel_sum; sum(diff(tmp_velocity(1:5)))];               
        end

    end

%%% PLOT ALL ANIMAL AP5/SALINE COHORT VELOCITY PROFILE CORRELATION AND SINGLE TRIAL PEAK VELOCITY

    figure; hold on;
    for day = 1:5
        scatter(day*ones(1,6),AP5_corr_day(:,day),30,[0 0 0],'filled')
        plot(AP5_corr_day','color','k')
    end
    ylim([0 1])

    figure; hold on;
    for day = 1:5
        scatter(day*ones(1,6),saline_corr_day(:,day),30,[0 0 0],'filled')
        plot(saline_corr_day','color','k')
    end
    ylim([0 1])

    figure; hold on;
    for day = 1:5
        scatter(day*ones(1,6),AP5_vel_day(:,day),30,[0 0 0],'filled')
        plot(AP5_vel_day','color','k')
    end
    ylim([20 55])

    figure; hold on;
    for day = 1:5
        scatter(day*ones(1,6),saline_vel_day(:,day),30,[0 0 0],'filled')
        plot(saline_vel_day','color','k')
    end
    ylim([20 55])

%%% PLOT AP5/SALINE INFUSION VELOCITY PROFILE CORRELATION COMPARISON

    figure;
        hold on;
        scatter(ones(1,6),saline_vel_sum,30,[.5 .5 .5],'LineWidth',2)
        errorbar(1.25,mean(saline_vel_sum),std(saline_vel_sum)/sqrt(6),'color','k','LineWidth',2)
        scatter(2*ones(1,6),AP5_vel_sum,30,[.5 .5 .5],'LineWidth',2)
        errorbar(1.75,mean(AP5_vel_sum),std(AP5_vel_sum)/sqrt(6),'color','k','LineWidth',2)
        plot([1.25 1.75],[mean(saline_vel_sum) mean(AP5_vel_sum)],'color','k','LineWidth',2)
        xlim([.5 2.5])
        [p h] = ranksum(saline_vel_sum,AP5_vel_sum)
        title(num2str(p))
        ylim([-20 20])

%%% PLOT INDIVIDUAL ANIMAL AP5/SALINE VELOCITY PROFILE CORRELATION AND SINGLE TRIAL PEAK VELOCITY
    
    for n = [1 3 6]

        figure;
            subplot(2,1,1); hold on;
                scatter([1:10], [saline_vel_day(n,:) AP5_vel_day(n,:)],30,[0 0 0],'filled');
                plot([saline_vel_day(n,:) AP5_vel_day(n,:)],'k')
                if n == 1
                    ylim([25 55])
                elseif n == 3
                    ylim([25 50])
                elseif n == 6
                    ylim([30 55])
                end
            subplot(2,1,2); hold on;
                scatter([1:10], [saline_corr_day(n,:) AP5_corr_day(n,:)],30,[0 0 0],'filled');        
                plot([saline_corr_day(n,:) AP5_corr_day(n,:)],'k')

    end

    for n = [2 4 5]

        figure;
            subplot(2,1,1); hold on;
                scatter([1:10], [AP5_vel_day(n,:) saline_vel_day(n,:)],30,[0 0 0],'filled');        
                plot([AP5_vel_day(n,:) saline_vel_day(n,:)],'k')
                if n == 2
                    ylim([25 45])
                elseif n == 4
                    ylim([30 55])
                elseif n == 5
                    ylim([20 40])
                end
            subplot(2,1,2); hold on;
                scatter([1:10], [AP5_corr_day(n,:) saline_corr_day(n,:)],30,[0 0 0],'filled');                
                plot([AP5_corr_day(n,:) saline_corr_day(n,:)],'k')

    end

%%% VELOCITY PROFILE CORRELATION HISTOGRAMS

    AP5_day2day = [];
    AP5_perAnimal = [];
    saline_day2day = [];
    saline_perAnimal = [];

    for animal = 1:6

        tmp = diff(AP5_corr_day(animal,:));
        AP5_day2day = [AP5_day2day tmp(~isnan(tmp))];
        AP5_perAnimal = [AP5_perAnimal mean(tmp(~isnan(tmp)))];

        tmp = diff(saline_corr_day(animal,:));
        saline_day2day = [saline_day2day tmp(~isnan(tmp))];
        saline_perAnimal = [saline_perAnimal mean(tmp(~isnan(tmp)))];

    end

    figure;
        subplot(2,2,1); hold on;
            for animal = 1:6
                plot([0 AP5_perAnimal(animal)],[animal animal],'k');
            end
            xlim([-.15 .15])
        subplot(2,2,2); hold on;
            for animal = 1:6
                plot([0 saline_perAnimal(animal)],[animal animal],'k');
            end
            xlim([-.15 .15])
        subplot(2,2,3); hold on;
            histogram(AP5_day2day,[-.4:.1:.4],'DisplayStyle','Stairs')
            ylim([0 12])            
        subplot(2,2,4); hold on;
            histogram(saline_day2day,[-.4:.1:.4],'DisplayStyle','Stairs')
            ylim([0 12])
            
%%% SINGLE TRIAL PEAK VELOCITY HISTOGRAMS

    AP5_day2day = [];
    AP5_perAnimal = [];
    saline_day2day = [];
    saline_perAnimal = [];

    for animal = 1:6

        tmp = diff(AP5_vel_day(animal,:));
        AP5_day2day = [AP5_day2day tmp(~isnan(tmp))];
        AP5_perAnimal = [AP5_perAnimal mean(tmp(~isnan(tmp)))];

        tmp = diff(saline_vel_day(animal,:));
        saline_day2day = [saline_day2day tmp(~isnan(tmp))];
        saline_perAnimal = [saline_perAnimal mean(tmp(~isnan(tmp)))];

    end

    figure;
        subplot(2,2,1); hold on;
            for animal = 1:6
                plot([0 AP5_perAnimal(animal)],[animal animal],'k');
            end
            xlim([-5 5])
        subplot(2,2,2); hold on;
            for animal = 1:6
                plot([0 saline_perAnimal(animal)],[animal animal],'k');
            end
            xlim([-5 5])
        subplot(2,2,3); hold on;
            histogram(AP5_day2day,[-20:5:20],'DisplayStyle','Stairs')
            ylim([0 12])
        subplot(2,2,4); hold on;
            histogram(saline_day2day,[-20:5:20],'DisplayStyle','Stairs')
            ylim([0 12])            
        
%% FIG 1 E & F | sFIG 2 & 4 | LEARNING & AP5/SALINE COHORT | VELOCITY PROFILE CORRELATION (INDIVIDUAL ANIMAL)

    total_saline_change = [];
    total_AP5_change = [];
    total_learning_change = [];
    
    %%% compute total change in velocity profile correlation in AP5 cohort

        AP5_x = cell(1,10);
        AP5_y = cell(1,10);
        AP5_xvel = cell(1,10);
        AP5_yvel = cell(1,10);

        for animal = 1:6
            for day = 1:10

                files = dir([data_path '\AP' num2str(animal) '\behavioral_data\day_' num2str(day) '\*csv']);

                %%% FIND PELLET POSITION
                pellet_position = [];
                for file = 1:length(files)
                    trial_kin = readmatrix([data_path '\AP' num2str(animal) '\behavioral_data\day_' num2str(day) '\' files(file).name]);
                    tmp_pellet_x = trial_kin(1:25,5);
                    tmp_pellet_y = trial_kin(1:25,6);
                    tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
                    tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
                    pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                    pellet_position = [pellet_position; pellet_location];
                end
                pellet_position = nanmean(pellet_position);

                %%% GET TRAJECTORIES AND VELOCITY PROFILES
                all_diffx = [];
                all_diffy = [];
                all_x = [];
                all_y = [];
                for file = 1:length(files)
                    trial_kin = readmatrix([data_path '\AP' num2str(animal) '\behavioral_data\day_' num2str(day) '\' files(file).name]);
                    %%% get paw trajectory and velocity
                    tmp_x = trial_kin(:,2);
                    tmp_y = trial_kin(:,3);
                    tmp_x(trial_kin(:,4)<.99) = NaN;
                    tmp_y(trial_kin(:,4)<.99) = NaN;
                    diff_x = diff(tmp_x);
                    diff_y = diff(tmp_y);
                    %%% find frame with paw closest to pellet
                    [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));
                    %%% skip bad trials
                    if ~(i<31 || i>size(trial_kin,1)-31)
                        all_diffx = [all_diffx; diff_x(i-30:i+30)'];
                        all_diffy = [all_diffy; diff_y(i-30:i+30)'];
                        all_x = [all_x; tmp_x(i-30:i+30)'-pellet_position(1)];
                        all_y = [all_y; tmp_y(i-30:i+30)'-pellet_position(2)];
                    end
                end
                AP5_x{1,day} = all_x;
                AP5_y{1,day} = all_y;
                AP5_xvel{1,day} = all_diffx;
                AP5_yvel{1,day} = all_diffy;

            end

            AP5_x_corr = [];
            AP5_y_corr = [];
            for n = 1:length(AP5_xvel)
                AP5_x_corr = [AP5_x_corr corr(smooth(nanmean(AP5_xvel{1,n}(:,:)),1),smooth(nanmean(AP5_xvel{1,end}(:,:)),1))];
                AP5_y_corr = [AP5_y_corr corr(smooth(nanmean(AP5_yvel{1,n}(:,:)),1),smooth(nanmean(AP5_yvel{1,end}(:,:)),1))];
            end
            AP5_corr = mean([AP5_x_corr; AP5_y_corr]);

            figure;
            hold on;
            plot(AP5_corr,'k')
            for day = 1:10
                scatter(day,AP5_corr(day),50,[0 0 0],'filled')
            end
            xlabel('days')
            ylabel('correlation')

            if (animal == 1 | animal==3 | animal==6)
                title(['Saline day 1-5, AP5 day 6-10 - velocity profile correlation'])
                total_saline_change = [total_saline_change sum(diff(AP5_corr(1:5)))];
                total_AP5_change = [total_AP5_change sum(diff(AP5_corr(6:9)))];
            else
                title(['AP5 day 1-5, Saline day 6-10 - velocity profile correlation'])                   
                total_AP5_change = [total_AP5_change sum(diff(AP5_corr(1:5)))];
                total_saline_change = [total_saline_change sum(diff(AP5_corr(6:9)))];   
            end
    %         
    %         if animal==1
    %             ylim([.7 1]) 
    %         elseif animal == 2
    %             ylim([.5 1])   
    %         elseif animal == 3
    %             ylim([.6 1])     
    %         elseif animal == 4
    %             ylim([0 1])  
    %         elseif animal == 5
    %             ylim([.6 1])     
    %         elseif animal == 6
    %             ylim([.6 1])        
    %         end

        end
    
    %%% compute total change in velocity profile correlation in learning cohort

        animals = {'LC3','LC4'};
        days = {1:8,2:6};
        for animal = 1:length(animals)

            animal_days = days{animal};
            xvel = cell(1,length(animal_days));
            yvel = cell(1,length(animal_days));

            day_count = 1;
            for day = animal_days

                files = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.csv']);

                %%% find mean pellet position across trials
                pellet_position = [];
                for file = 1:length(files)
                    trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);
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

                    trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);

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
            total_learning_change = [total_learning_change sum(diff(animal_corr(1:end-1)))];

            figure;
            hold on;
            plot(animal_corr,'k')
            for day = 1:length(animal_corr)
                scatter(day,animal_corr(day),50,[0 0 0],'filled')
            end
            xlabel('days')
            ylabel('correlation')


        end

        animals = {'LC2','LC1'};
        days = {4:13,1:8};
        for animal = 1:length(animals)

            animal_days = days{animal};
            xvel = cell(1,length(animal_days));
            yvel = cell(1,length(animal_days));

            day_count = 1;
            for day = animal_days

                files = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.csv']);

                %%% find mean pellet position across trials
                pellet_position = [];
                for file = 1:length(files)
                    trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);
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

                    trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);

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
            total_learning_change = [total_learning_change sum(diff(animal_corr(1:end-1)))];

            figure;
            hold on;
            plot(animal_corr,'k')
            for day = 1:length(animal_corr)
                scatter(day,animal_corr(day),50,[0 0 0],'filled')
            end        
            xlabel('days')
            ylabel('correlation')

        end    

        animals = {'LC6','LC5'};
        days = {2:10,4:17};
        for animal = 1:length(animals)

            animal_days = days{animal};
            xvel = cell(1,length(animal_days));
            yvel = cell(1,length(animal_days));

            day_count = 1;
            for day = animal_days

                files = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.csv']);
                files = sort_nat({files.name});
                vid_timestamp = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.mat']); 
                vid_timestamp = sort_nat({vid_timestamp.name});

                all_diffx = [];
                all_diffy = [];

                for block = 1:length(files)

                    %%% load all trial kinematics
                    kin_all = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files{block}]);
                    load([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' vid_timestamp{block}]);

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
            total_learning_change = [total_learning_change sum(diff(animal_corr(1:end-1)))];

            figure;
            hold on;
            plot(animal_corr,'k')
            for day = 1:length(animal_corr)
                scatter(day,animal_corr(day),50,[0 0 0],'filled')
            end                 
            xlabel('days')
            ylabel('correlation')

        end
    
    %%% plot comparison

        fig = figure;
        hold on;

            scatter(ones(1,length(total_saline_change)),total_saline_change,30,[0.5 0.5 0.5],'filled')
            errorbar(1.5,mean(total_saline_change),std(total_saline_change)/sqrt(length(total_saline_change)),'color','k')

            scatter(2*ones(1,length(total_AP5_change)),total_AP5_change,30,[0.5 0.5 0.5],'filled')
            errorbar(2.5,mean(total_AP5_change),std(total_AP5_change)/sqrt(length(total_AP5_change)),'color','k')

            scatter(3*ones(1,length(total_learning_change)),total_learning_change,30,[0.5 0.5 0.5],'filled')
            errorbar(3.5,mean(total_learning_change),std(total_learning_change)/sqrt(length(total_learning_change)),'color','k')

            xlim([0 4.5])

            [p1, h, stats] = ranksum(total_saline_change,total_AP5_change)
            [p2, h, stats] = ranksum(total_learning_change,total_saline_change)
            [p3, h, stats] = ranksum(total_learning_change,total_AP5_change)

            title(['SUM -> saline vs. ap5 p=' num2str(p1) ' - saline vs. learning p=' num2str(p2) ' - ap5 vs. learning p=' num2str(p3)])

            mean(total_AP5_change)
            std(total_AP5_change)/sqrt(6)

            mean(total_saline_change)
            std(total_saline_change)/sqrt(6)

            mean(total_learning_change)
            std(total_learning_change)/sqrt(6)

%% sFIG 3 | LEARNING COHORT | PEAK VELOCITY TIMECOURSE OVER LEARNING
    
%%% plot individual animal single trial peak velocity over learning

    animals = {'LC3','LC4'};
    days = {[1:8],[2:6]};

    for animal = 1:length(animals)

        animal_days = days{animal};
        day_velocity = cell(1,length(animal_days));
        day_count = 1;

        for day = animal_days

            files = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.csv']);

            %%% find mean pellet position across trials
            pellet_position = [];
            for file = 1:length(files)
                trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);
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
                trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);

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
            %%% 0.0711 cm/pixel
            %%% 30 frames/second
            day_velocity{day_count} = (max(((abs(all_diffx)+abs(all_diffy))/2)')*0.0711)*30;
            day_count = day_count + 1;
        end

        figure;
        hold on;
        tmp_velocity = [];
        for day = 1:length(day_velocity)
            errorbar(day,nanmean(day_velocity{day}),nanstd(day_velocity{day})/sqrt(length(day_velocity{day})),'color','k')
            tmp_velocity = [tmp_velocity nanmean(day_velocity{day})];
        end
        plot(tmp_velocity,'k')
        xlabel('days')
        ylabel('velocity')

    end

    animals = {'LC2'};
    days = {[4:13]};

    for animal = 1:length(animals)

        animal_days = days{animal};
        day_velocity = cell(1,length(animal_days));
        day_count = 1;

        for day = animal_days

            files = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.csv']);

            %%% find mean pellet position across trials
            pellet_position = [];
            for file = 1:length(files)
                trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);
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
                trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);
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

            %%% 0.0237 cm/pixel
            %%% 75 frames/second        
            day_velocity{day_count} = (max(((abs(all_diffx)+abs(all_diffy))/2)')*0.0237)*75;
            day_count = day_count + 1;
        end

        figure;
        hold on;
        tmp_velocity = [];
        for day = 1:length(day_velocity)
            errorbar(day,nanmean(day_velocity{day}),nanstd(day_velocity{day})/sqrt(length(day_velocity{day})),'color','k')
            tmp_velocity = [tmp_velocity nanmean(day_velocity{day})];
        end
        plot(tmp_velocity,'k')
        xlabel('days')
        ylabel('velocity')

    end

    animals = {'LC1'};
    days = {[1:8]};

    for animal = 1:length(animals)

        animal_days = days{animal};
        day_velocity = cell(1,length(animal_days));
        day_count = 1;

        for day = animal_days

            files = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.csv']);

            %%% find mean pellet position across trials
            pellet_position = [];
            for file = 1:length(files)
                trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);
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
                trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);
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

            %%% 0.0127 cm/pixel
            %%% 100 frames/second
            day_velocity{day_count} = (max(((abs(all_diffx)+abs(all_diffy))/2)')*0.0127)*100;

            day_count = day_count + 1;
        end

        figure;
        hold on;
        tmp_velocity = [];
        for day = 1:length(day_velocity)
            errorbar(day,nanmean(day_velocity{day}),nanstd(day_velocity{day})/sqrt(length(day_velocity{day})),'color','k')
            tmp_velocity = [tmp_velocity nanmean(day_velocity{day})];
        end
        plot(tmp_velocity,'k')
        xlabel('days')
        ylabel('velocity')

    end

    animals = {'LC6'};
    days = {[2:10]};

    for animal = 1:length(animals)

        animal_days = days{animal};
        day_velocity = cell(1,length(animal_days));
        day_count = 1;

         for day = animal_days

            files = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.csv']);
            files = sort_nat({files.name});
            vid_timestamp = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.mat']);
            vid_timestamp = sort_nat({vid_timestamp.name});

            all_diffx = [];
            all_diffy = [];
            for block = 1:length(files)

                %%% load all trial kinematics
                kin_all = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files{block}]);
                load([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' vid_timestamp{block}]);

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

            %%% 0.0222 cm/pixel
            %%% 60 frames/second
            day_velocity{day_count} = (max(((abs(all_diffx)+abs(all_diffy))/2)')*0.0222)*60;

            day_count = day_count + 1;
        end

        figure;
        hold on;
        tmp_velocity = [];
        for day = 1:length(day_velocity)
            errorbar(day,nanmean(day_velocity{day}),nanstd(day_velocity{day})/sqrt(length(day_velocity{day})),'color','k')
            tmp_velocity = [tmp_velocity nanmean(day_velocity{day})];
        end
        plot(tmp_velocity,'k')
        xlabel('days')
        ylabel('velocity')

    end

    animals = {'LC5'};
    days = {[4:17]};

    for animal = 1:length(animals)

        animal_days = days{animal};
        day_velocity = cell(1,length(animal_days));
        day_count = 1;

        for day = animal_days

            files = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.csv']);
            files = sort_nat({files.name});
            vid_timestamp = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.mat']);
            vid_timestamp = sort_nat({vid_timestamp.name});

            all_diffx = [];
            all_diffy = [];
            for block = 1:length(files)

                %%% load all trial kinematics
                kin_all = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files{block}]);
                load([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' vid_timestamp{block}]);

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

            %%% 0.0254 cm/pixel
            %%% 60 frames/second
            day_velocity{day_count} = (max(((abs(all_diffx)+abs(all_diffy))/2)')*0.0254)*60;

            day_count = day_count + 1;
        end

        figure;
        hold on;
        tmp_velocity = [];
        for day = 1:length(day_velocity)
            errorbar(day,nanmean(day_velocity{day}),nanstd(day_velocity{day})/sqrt(length(day_velocity{day})),'color','k')
            tmp_velocity = [tmp_velocity nanmean(day_velocity{day})];
        end
        plot(tmp_velocity,'k')
        xlabel('days')
        ylabel('velocity')    

    end

%% sFIG 3 | LEARNING COHORT | VELOCITY PROFILE CORRELATION VS. SINGLE TRIAL PEAK VELOCITY

learning_corr_day = [];
learning_vel_day = [];

animals = {'LC3','LC4'};
days = {1:8,2:6};
for animal = 1:length(animals)
    
    animal_days = days{animal};
    xvel = cell(1,length(animal_days));
    yvel = cell(1,length(animal_days));
    day_velocity = cell(1,length(animal_days));

    day_count = 1;
    for day = animal_days
        
        files = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.csv']);
        
        %%% find mean pellet position across trials
        pellet_position = [];
        for file = 1:length(files)
            trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);
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
            trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);
            
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
        day_velocity{day_count} = (max(((abs(all_diffx)+abs(all_diffy))/2)')*0.0711)*30;
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
    
    tmp_velocity = [];
    for day = 1:length(day_velocity)
        tmp_velocity = [tmp_velocity nanmean(day_velocity{day})];
    end
    
    learning_corr_day = [learning_corr_day animal_corr];
    learning_vel_day = [learning_vel_day tmp_velocity];

end

animals = {'LC2'};
days = {4:13};
for animal = 1:length(animals)
    
    animal_days = days{animal};
    xvel = cell(1,length(animal_days));
    yvel = cell(1,length(animal_days));
    day_velocity = cell(1,length(animal_days));

    day_count = 1;
    for day = animal_days
        
        files = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.csv']);
        
        %%% find mean pellet position across trials
        pellet_position = [];
        for file = 1:length(files)
            trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);
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
            
            trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);
            
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
        day_velocity{day_count} = (max(((abs(all_diffx)+abs(all_diffy))/2)')*0.0237)*75;
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
    
    tmp_velocity = [];
    for day = 1:length(day_velocity)
        tmp_velocity = [tmp_velocity nanmean(day_velocity{day})];
    end
    
    learning_corr_day = [learning_corr_day animal_corr];
    learning_vel_day = [learning_vel_day tmp_velocity];

end

animals = {'LC1'};
days = {[1:8]};
for animal = 1:length(animals)
    
    animal_days = days{animal};
    xvel = cell(1,length(animal_days));
    yvel = cell(1,length(animal_days));
    day_velocity = cell(1,length(animal_days));
    
    day_count = 1;
    for day = animal_days
        
        files = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.csv']);
        
        %%% find mean pellet position across trials
        pellet_position = [];
        for file = 1:length(files)
            trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);
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
            
            trial_kin = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files(file).name]);
            
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
        day_velocity{day_count} = (max(((abs(all_diffx)+abs(all_diffy))/2)')*0.0127)*100;
        
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
    
    tmp_velocity = [];
    for day = 1:length(day_velocity)
        tmp_velocity = [tmp_velocity nanmean(day_velocity{day})];
    end
    
    learning_corr_day = [learning_corr_day animal_corr];
    learning_vel_day = [learning_vel_day tmp_velocity];

end

animals = {'LC6'};
days = {2:10};
for animal = 1:length(animals)
    
    animal_days = days{animal};
    xvel = cell(1,length(animal_days));
    yvel = cell(1,length(animal_days));
    day_velocity = cell(1,length(animal_days));

    day_count = 1;
    for day = animal_days
        
        files = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.csv']);
        files = sort_nat({files.name});
        vid_timestamp = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.mat']);
        vid_timestamp = sort_nat({vid_timestamp.name});
        
        all_diffx = [];
        all_diffy = [];
        
        for block = 1:length(files)
            
            %%% load all trial kinematics
            kin_all = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files{block}]);
            load([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' vid_timestamp{block}]);
            
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
        day_velocity{day_count} = (max(((abs(all_diffx)+abs(all_diffy))/2)')*0.0222)*60;
        
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
    
    tmp_velocity = [];
    for day = 1:length(day_velocity)
        tmp_velocity = [tmp_velocity nanmean(day_velocity{day})];
    end
    
    learning_corr_day = [learning_corr_day animal_corr];
    learning_vel_day = [learning_vel_day tmp_velocity];

end

animals = {'LC5'};
days = {4:17};
for animal = 1:length(animals)
    
    animal_days = days{animal};
    xvel = cell(1,length(animal_days));
    yvel = cell(1,length(animal_days));
    day_velocity = cell(1,length(animal_days));

    day_count = 1;
    for day = animal_days
        
        files = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.csv']);
        files = sort_nat({files.name});
        vid_timestamp = dir([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\*.mat']);
        vid_timestamp = sort_nat({vid_timestamp.name});
        
        all_diffx = [];
        all_diffy = [];
        
        for block = 1:length(files)
            
            %%% load all trial kinematics
            kin_all = readmatrix([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' files{block}]);
            load([data_path '\' animals{animal} '\reach_trajectories\day_' num2str(day) '\' vid_timestamp{block}]);
            
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
        day_velocity{day_count} = (max(((abs(all_diffx)+abs(all_diffy))/2)')*0.0254)*60;
        
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
    
    tmp_velocity = [];
    for day = 1:length(day_velocity)
        tmp_velocity = [tmp_velocity nanmean(day_velocity{day})];
    end
    
    learning_corr_day = [learning_corr_day animal_corr];
    learning_vel_day = [learning_vel_day tmp_velocity];

end

learning_vel_day = learning_vel_day(learning_corr_day<.99);
learning_corr_day = learning_corr_day(learning_corr_day<.99);

figure;
hold on;
scatter(learning_corr_day,learning_vel_day,30,[.5 .5 .5],'filled')
[R,P] = corrcoef(learning_corr_day', learning_vel_day')
[p] = polyfit(learning_corr_day', learning_vel_day',1)
x=[.3:.1:1]
y = x*p(1)+p(2)
plot(x,y,'color','k','LineStyle','--','LineWidth',3)
title(['stab vs. vel - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);

%% sFIG 5 | AP5/SALINE COHORTS | DIRECT COMPARISON / PERMUTATIONS

    all_animal_corr = cell(1,6);
    all_animal_vel = cell(1,6);

    for animal = 1:6

        for day = 1:10

            files = dir([data_path '\AP' num2str(animal) '\behavioral_data\day_' num2str(day) '\*csv']);

            %%% FIND PELLET POSITION
            pellet_position = [];
            for file = 1:length(files)
                trial_kin = readmatrix([data_path '\AP' num2str(animal) '\behavioral_data\day_' num2str(day) '\' files(file).name]);
                tmp_pellet_x = trial_kin(1:25,5);
                tmp_pellet_y = trial_kin(1:25,6);
                tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
                pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                pellet_position = [pellet_position; pellet_location];
            end
            pellet_position = nanmean(pellet_position);

            %%% GET TRAJECTORIES AND VELOCITY PROFILES
            all_diffx = [];
            all_diffy = [];
            all_x = [];
            all_y = [];
            for file = 1:length(files)
                trial_kin = readmatrix([data_path '\AP' num2str(animal) '\behavioral_data\day_' num2str(day) '\' files(file).name]);
                %%% get paw trajectory and velocity
                tmp_x = trial_kin(:,2);
                tmp_y = trial_kin(:,3);
                tmp_x(trial_kin(:,4)<.99) = NaN;
                tmp_y(trial_kin(:,4)<.99) = NaN;
                diff_x = diff(tmp_x);
                diff_y = diff(tmp_y);
                %%% find frame with paw closest to pellet
                [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));
                %%% skip bad trials
                if ~(i<31 || i>size(trial_kin,1)-31)
                    all_diffx = [all_diffx; diff_x(i-30:i+30)'];
                    all_diffy = [all_diffy; diff_y(i-30:i+30)'];
                    all_x = [all_x; tmp_x(i-30:i+30)'-pellet_position(1)];
                    all_y = [all_y; tmp_y(i-30:i+30)'-pellet_position(2)];
                end
            end
            AP5_x{1,day} = all_x;
            AP5_y{1,day} = all_y;
            AP5_xvel{1,day} = all_diffx;
            AP5_yvel{1,day} = all_diffy;
            AP5_maxVel{1,day} = (max(((abs(all_diffx)+abs(all_diffy))/2)')*0.0237)*75;

        end

        AP5_x_corr = [];
        AP5_y_corr = [];
        for n = 1:length(AP5_xvel)
            AP5_x_corr = [AP5_x_corr corr(smooth(nanmean(AP5_xvel{1,n}(:,:)),1),smooth(nanmean(AP5_xvel{1,end}(:,:)),1))];
            AP5_y_corr = [AP5_y_corr corr(smooth(nanmean(AP5_yvel{1,n}(:,:)),1),smooth(nanmean(AP5_yvel{1,end}(:,:)),1))];
        end
        AP5_corr = mean([AP5_x_corr; AP5_y_corr]);

        tmp_velocity = [];
        for day = 1:length(AP5_maxVel)
            tmp_velocity = [tmp_velocity nanmean(AP5_maxVel{day})];
        end

        tmp_velocity_norm = [];
        for day = 1:length(AP5_maxVel)
            tmp_velocity_norm = [tmp_velocity_norm tmp_velocity(day)/tmp_velocity(end)];
        end

        all_animal_corr{animal} = AP5_corr;
        all_animal_vel{animal} = tmp_velocity;

    end

    saline_corr_diff = [];
    saline_vel_diff = [];
    AP5_corr_diff = [];
    AP5_vel_diff = [];

    for animal = [1 3 6]

        all_animal_vel_norm = all_animal_vel{animal}/all_animal_vel{animal}(end);

        saline_corr_diff = [saline_corr_diff diff(all_animal_corr{animal}(1:5))];    
        saline_vel_diff = [saline_vel_diff diff(all_animal_vel_norm(1:5))];

        AP5_corr_diff = [AP5_corr_diff diff(all_animal_corr{animal}(6:9))];    
        AP5_vel_diff = [AP5_vel_diff diff(all_animal_vel_norm(6:10))];

    end

    for animal = [2 4 5]

        all_animal_vel_norm = all_animal_vel{animal}/all_animal_vel{animal}(end);

        saline_corr_diff = [saline_corr_diff diff(all_animal_corr{animal}(6:9))];
        saline_vel_diff = [saline_vel_diff diff(all_animal_vel_norm(6:10))];

        AP5_corr_diff = [AP5_corr_diff diff(all_animal_corr{animal}(1:5))];    
        AP5_vel_diff = [AP5_vel_diff diff(all_animal_vel_norm(1:5))];

    end

%%% COMPUTE CORR: AP5 vs. SALINE

    tmp_vals = [saline_corr_diff AP5_corr_diff];

    Sh_saline = [];
    Sh_AP5 = [];
    for shuffles = 1:10000
        tmpIdx = randperm(length([saline_corr_diff AP5_corr_diff]));
        tmp_saline_corr_diff = tmp_vals(tmpIdx(1:21));
        tmp_AP5_corr_diff = tmp_vals(tmpIdx(22:42));

        Sh_saline = [Sh_saline; mean(tmp_saline_corr_diff)];
        Sh_AP5 = [Sh_AP5; mean(tmp_AP5_corr_diff)];
    end

    figure; hold on;
    histogram(Sh_saline-Sh_AP5,'DisplayStyle','Stairs')
    plot([mean(AP5_corr_diff)-mean(saline_corr_diff) mean(AP5_corr_diff)-mean(saline_corr_diff)],[ylim])
    title(num2str(sum((Sh_saline-Sh_AP5)<(mean(AP5_corr_diff)-mean(saline_corr_diff)))/10000))

%%% COMPUTE VEL: AP5 vs. SALINE

    tmp_vals = [saline_vel_diff AP5_vel_diff];

    Sh_saline = [];
    Sh_AP5 = [];
    for shuffles = 1:10000
        tmpIdx = randperm(length([saline_vel_diff AP5_vel_diff]));
        tmp_saline_vel_diff = tmp_vals(tmpIdx(1:24));
        tmp_AP5_vel_diff = tmp_vals(tmpIdx(25:48));

        Sh_saline = [Sh_saline; mean(tmp_saline_vel_diff)];
        Sh_AP5 = [Sh_AP5; mean(tmp_AP5_vel_diff)];
    end

    figure; hold on;
    histogram(Sh_saline-Sh_AP5,'DisplayStyle','Stairs')
    plot([mean(AP5_vel_diff)-mean(saline_vel_diff) mean(AP5_vel_diff)-mean(saline_vel_diff)],[ylim])
    title(num2str(sum((Sh_saline-Sh_AP5)<(mean(AP5_vel_diff)-mean(saline_vel_diff)))/10000))

