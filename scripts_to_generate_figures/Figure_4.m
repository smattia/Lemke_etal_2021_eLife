%% SCRIPT TO GENERATE FIGURE 4 & FIGURE 4 FIGURE SUPPLEMENTS
% 
% to generate figures, download data from https://doi.org/10.7272/Q6KK9927

clc; clear all; close all;
data_path = 'D:\Lemke_etal_2021\Lemke_eLife_2021_data';

%% FIG 4 A | PLOT INDIVIDUAL ANIMAL/DAY DLS PCA TRAJECTORIES AND M1 PREDICTION OF DLS TRAJECTORIES

    animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
    animal_days = {[1:8],[2:6],[2:10],[4:17],[4:13],[3:5 7:9]};

    for animal = 1:length(animal_id)

        animal_path = [data_path '\' animal_id{animal}];

        for day = animal_days{animal}

            load([animal_path '\reach_activity\day_' num2str(day) '_reach_mod.mat']);

            %%% DLS PCA

            all_dls_raster = [];
            for unit = 1:length(dls_reach)
                tmp_dls_raster = [];
                for trial = 1:size(dls_reach(unit).raster_100ms,1)
                    tmp_dls_raster = [tmp_dls_raster dls_reach(unit).raster_100ms(trial,:)];
                end
                all_dls_raster = [all_dls_raster; tmp_dls_raster];
            end
            [coeff,score,latent] = pca(all_dls_raster');
            dls_pc1 = reshape(score(:,1),[size(dls_reach(unit).raster_100ms,2) size(dls_reach(unit).raster_100ms,1)])';
            dls_pc2 = reshape(score(:,2),[size(dls_reach(unit).raster_100ms,2) size(dls_reach(unit).raster_100ms,1)])';

            %%%

            trial_ind = 1:size(m1_reach(1).raster_100ms,1);
            cv_sets = [1 ...
                round(size(m1_reach(1).raster_100ms,1)/5) ...
                round(size(m1_reach(1).raster_100ms,1)*2/5) ...
                round(size(m1_reach(1).raster_100ms,1)*3/5) ...
                round(size(m1_reach(1).raster_100ms,1)*4/5) ...
                size(m1_reach(1).raster_100ms,1)];

            tmp_task1 = [];
            tmp_pred1 = [];

            tmp_task2 = [];
            tmp_pred2 = [];

            for cv_n = 1:5

                test_trials = zeros(1,size(m1_reach(1).raster_100ms,1));
                training_trials = ones(1,size(m1_reach(1).raster_100ms,1));
                test_trials(trial_ind(cv_sets(cv_n):cv_sets(cv_n+1)))=1;
                training_trials(trial_ind(cv_sets(cv_n):cv_sets(cv_n+1)))=0;

                %%% TRAINING DATA

                task_related_dls1 = [];
                task_related_dls2 = [];
                for trial = find(training_trials==1)
                    task_related_dls1 = [task_related_dls1; dls_pc1(trial,41:60)];
                    task_related_dls2 = [task_related_dls2; dls_pc2(trial,41:60)];
                end
                task_related_dls1 = mean(task_related_dls1)';
                task_related_dls2 = mean(task_related_dls2)';

                task_related_m1 = [];
                for time_step = 41:60
                    tmp_task_related_m1 = [];
                    for unit=1:length(m1_reach)
                        tmp_raster = [];
                        for trial = find(training_trials==1)
                            tmp_raster = [tmp_raster; m1_reach(unit).raster_100ms(trial,time_step-15:time_step)];
                        end
                        tmp_task_related_m1 = [tmp_task_related_m1 mean(tmp_raster)];
                    end
                    task_related_m1 = [task_related_m1; tmp_task_related_m1];
                end

                tbl = array2table([task_related_m1 task_related_dls1]);
                mdl1 = fitlm(tbl);

                tbl = array2table([task_related_m1 task_related_dls2]);
                mdl2 = fitlm(tbl);

                %%% TESTING DATA

                task_related_dls1 = [];
                task_related_dls2 = [];
                for trial = find(test_trials==1)
                    task_related_dls1 = [task_related_dls1; dls_pc1(trial,31:60)];
                    task_related_dls2 = [task_related_dls2; dls_pc2(trial,31:60)];
                end
                task_related_dls1 = mean(task_related_dls1)';
                task_related_dls2 = mean(task_related_dls2)';

                task_related_m1 = [];
                for time_step = 31:60
                    tmp_task_related_m1 = [];
                    for unit=1:length(m1_reach)
                        tmp_raster = [];
                        for trial = find(test_trials==1)
                            tmp_raster = [tmp_raster; m1_reach(unit).raster_100ms(trial,time_step-15:time_step)];
                        end
                        tmp_task_related_m1 = [tmp_task_related_m1 mean(tmp_raster)];
                    end
                    task_related_m1 = [task_related_m1; tmp_task_related_m1];
                end

                ypred1 = predict(mdl1,task_related_m1);
                ypred2 = predict(mdl2,task_related_m1);

                tmp_task1 = [tmp_task1; task_related_dls1'];
                tmp_task2 = [tmp_task2; task_related_dls2'];

                tmp_pred1 = [tmp_pred1; ypred1'];
                tmp_pred2 = [tmp_pred2; ypred2'];

            end

            fig = figure;

            subplot(1,2,1)
            hold on

            plot(smooth(mean(tmp_task1(:,1:end)),5),smooth(mean(tmp_task2(:,1:end)),5),'color','k')

            tmp1 = smooth(mean(tmp_task1(:,1:end)),5);
            tmp2 = smooth(mean(tmp_task2(:,1:end)),5);
            scatter(tmp1(1),tmp2(1),30,[0 1 0])
            scatter(tmp1(11),tmp2(11),30,[0 0 0])
            scatter(tmp1(end),tmp2(end),30,[1 0 0])

            subplot(1,2,2)
            hold on

            plot(smooth(mean(tmp_pred1(:,1:end)),5),smooth(mean(tmp_pred2(:,1:end)),5),'color','k')

            tmp1 = smooth(mean(tmp_pred1(:,1:end)),5);
            tmp2 = smooth(mean(tmp_pred2(:,1:end)),5);
            scatter(tmp1(1),tmp2(1),30,[0 1 0])
            scatter(tmp1(11),tmp2(11),30,[0 0 0])
            scatter(tmp1(end),tmp2(end),30,[1 0 0])

        end

    end
    
%% FIG 4 B | PLOT M1 PREDICTION OF DLS TRAJECTORIES OVER THE COURSE OF LEARNING

animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
animal_days = {[1:8],[2:6],[2:10],[4:17],[4:13],[3:5 7:9]};
animal_curves = cell(1,length(animal_id));
baseline_curves = cell(1,length(animal_id));

for animal = 1:length(animal_id)
    
    animal_path = [data_path '\' animal_id{animal}];
    
    %%% REACH ACTIVITY

        tmp_corr1 = [];
        tmp_corr2 = [];
        tmp_corr3 = [];
        num_units = [];

        for day = animal_days{animal}

            load([animal_path '\reach_activity\day_' num2str(day) '_reach_mod.mat']);

            if length(dls_reach)>6 && length(m1_reach)>6

                %%% GENERATE PCs
                all_dls_raster = [];
                for unit = 1:length(dls_reach)
                    tmp_dls_raster = [];
                    for trial = 1:size(dls_reach(unit).raster_100ms,1)
                        tmp_dls_raster = [tmp_dls_raster dls_reach(unit).raster_100ms(trial,:)];
                    end
                    all_dls_raster = [all_dls_raster; tmp_dls_raster];
                end
                [coeff,score,latent] = pca(all_dls_raster');
                dls_pc1 = reshape(score(:,1),[size(dls_reach(unit).raster_100ms,2) size(dls_reach(unit).raster_100ms,1)])';
                dls_pc2 = reshape(score(:,2),[size(dls_reach(unit).raster_100ms,2) size(dls_reach(unit).raster_100ms,1)])';
                dls_pc3 = reshape(score(:,3),[size(dls_reach(unit).raster_100ms,2) size(dls_reach(unit).raster_100ms,1)])';

                %%% 5 FOLD CROSS VALIDATION
                trial_ind = 1:size(m1_reach(1).raster_100ms,1);
                cv_sets = [1 ...
                    round(size(m1_reach(1).raster_100ms,1)/5) ...
                    round(size(m1_reach(1).raster_100ms,1)*2/5) ...
                    round(size(m1_reach(1).raster_100ms,1)*3/5) ...
                    round(size(m1_reach(1).raster_100ms,1)*4/5) ...
                    size(m1_reach(1).raster_100ms,1)];

                cv_corr1 = [];
                cv_corr2 = [];
                cv_corr3 = [];

                for cv_n = 1:5

                    test_trials = zeros(1,size(m1_reach(1).raster_100ms,1));
                    training_trials = ones(1,size(m1_reach(1).raster_100ms,1));
                    test_trials(trial_ind(cv_sets(cv_n):cv_sets(cv_n+1)))=1;
                    training_trials(trial_ind(cv_sets(cv_n):cv_sets(cv_n+1)))=0;

                    %%% TRAINING DATA

                    task_related_dls1 = [];
                    task_related_dls2 = [];
                    task_related_dls3 = [];
                    for trial = find(training_trials==1)
                        task_related_dls1 = [task_related_dls1; dls_pc1(trial,41:60)];
                        task_related_dls2 = [task_related_dls2; dls_pc2(trial,41:60)];
                        task_related_dls3 = [task_related_dls3; dls_pc3(trial,41:60)];
                    end
                    task_related_dls1 = mean(task_related_dls1)';
                    task_related_dls2 = mean(task_related_dls2)';
                    task_related_dls3 = mean(task_related_dls3)';

                    task_related_m1 = [];
                    for time_step = 41:60
                        tmp_task_related_m1 = [];
                        for unit=1:length(m1_reach)
                            tmp_raster = [];
                            for trial = find(training_trials==1)
                                tmp_raster = [tmp_raster; m1_reach(unit).raster_100ms(trial,time_step-15:time_step)];
                            end
                            tmp_task_related_m1 = [tmp_task_related_m1 mean(tmp_raster)];
                        end
                        task_related_m1 = [task_related_m1; tmp_task_related_m1];
                    end
                    tbl = array2table([task_related_m1 task_related_dls1]);
                    mdl1 = fitlm(tbl);
                    tbl = array2table([task_related_m1 task_related_dls2]);
                    mdl2 = fitlm(tbl);
                    tbl = array2table([task_related_m1 task_related_dls3]);
                    mdl3 = fitlm(tbl);

                    %%% TESTING DATA

                    task_related_dls1 = [];
                    task_related_dls2 = [];
                    task_related_dls3 = [];
                    for trial = find(test_trials==1)
                        task_related_dls1 = [task_related_dls1; dls_pc1(trial,41:60)];
                        task_related_dls2 = [task_related_dls2; dls_pc2(trial,41:60)];
                        task_related_dls3 = [task_related_dls3; dls_pc3(trial,41:60)];
                    end
                    task_related_dls1 = mean(task_related_dls1)';
                    task_related_dls2 = mean(task_related_dls2)';
                    task_related_dls3 = mean(task_related_dls3)';

                    task_related_m1 = [];
                    for time_step = 41:60
                        tmp_task_related_m1 = [];
                        for unit=1:length(m1_reach)
                            tmp_raster = [];
                            for trial = find(test_trials==1)
                                tmp_raster = [tmp_raster; m1_reach(unit).raster_100ms(trial,time_step-15:time_step)];
                            end
                            tmp_task_related_m1 = [tmp_task_related_m1 mean(tmp_raster)];
                        end
                        task_related_m1 = [task_related_m1; tmp_task_related_m1];
                    end

                    ypred1 = predict(mdl1,task_related_m1);
                    ypred2 = predict(mdl2,task_related_m1);
                    ypred3 = predict(mdl3,task_related_m1);

                    [r p] = corrcoef(ypred1,task_related_dls1);
                    cv_corr1 = [cv_corr1; r(1,2)];

                    [r p] = corrcoef(ypred2,task_related_dls2);
                    cv_corr2 = [cv_corr2; r(1,2)];

                    [r p] = corrcoef(ypred3,task_related_dls3);
                    cv_corr3 = [cv_corr3; r(1,2)];

                end

                tmp_corr1 = [tmp_corr1; mean(cv_corr1)];
                tmp_corr2 = [tmp_corr2; mean(cv_corr2)];
                tmp_corr3 = [tmp_corr3; mean(cv_corr3)];

            else

                tmp_corr1 = [tmp_corr1; NaN];
                tmp_corr2 = [tmp_corr2; NaN];
                tmp_corr3 = [tmp_corr3; NaN];

            end

            num_units = [num_units; length(dls_reach) length(m1_reach)];

        end

        animal_curves{animal} = [tmp_corr1'; tmp_corr2'; tmp_corr3'];

    %%% BASELINE

        tmp_corr1 = [];
        tmp_corr2 = [];
        tmp_corr3 = [];
        num_units = [];

        for day = animal_days{animal}

            load([animal_path '\reach_activity\day_' num2str(day) '_reach_mod.mat']);

            if length(dls_reach)>6 && length(m1_reach)>6

                %%% GENERATE PCs
                all_dls_raster = [];
                for unit = 1:length(dls_reach)
                    tmp_dls_raster = [];
                    for trial = 1:size(dls_reach(unit).raster_100ms,1)
                        tmp_dls_raster = [tmp_dls_raster dls_reach(unit).raster_100ms(trial,:)];
                    end
                    all_dls_raster = [all_dls_raster; tmp_dls_raster];
                end
                [coeff,score,latent] = pca(all_dls_raster');
                dls_pc1 = reshape(score(:,1),[size(dls_reach(unit).raster_100ms,2) size(dls_reach(unit).raster_100ms,1)])';
                dls_pc2 = reshape(score(:,2),[size(dls_reach(unit).raster_100ms,2) size(dls_reach(unit).raster_100ms,1)])';
                dls_pc3 = reshape(score(:,3),[size(dls_reach(unit).raster_100ms,2) size(dls_reach(unit).raster_100ms,1)])';

                %%% 5 FOLD CROSS VALIDATION
                trial_ind = 1:size(m1_reach(1).raster_100ms,1);
                cv_sets = [1 ...
                    round(size(m1_reach(1).raster_100ms,1)/5) ...
                    round(size(m1_reach(1).raster_100ms,1)*2/5) ...
                    round(size(m1_reach(1).raster_100ms,1)*3/5) ...
                    round(size(m1_reach(1).raster_100ms,1)*4/5) ...
                    size(m1_reach(1).raster_100ms,1)];

                cv_corr1 = [];
                cv_corr2 = [];
                cv_corr3 = [];

                for cv_n = 1:5

                    test_trials = zeros(1,size(m1_reach(1).raster_100ms,1));
                    training_trials = ones(1,size(m1_reach(1).raster_100ms,1));
                    test_trials(trial_ind(cv_sets(cv_n):cv_sets(cv_n+1)))=1;
                    training_trials(trial_ind(cv_sets(cv_n):cv_sets(cv_n+1)))=0;

                    %%% TRAINING DATA

                    task_related_dls1 = [];
                    task_related_dls2 = [];
                    task_related_dls3 = [];
                    for trial = find(training_trials==1)
                        task_related_dls1 = [task_related_dls1; dls_pc1(trial,16:36)];
                        task_related_dls2 = [task_related_dls2; dls_pc2(trial,16:36)];
                        task_related_dls3 = [task_related_dls3; dls_pc3(trial,16:36)];
                    end
                    task_related_dls1 = mean(task_related_dls1)';
                    task_related_dls2 = mean(task_related_dls2)';
                    task_related_dls3 = mean(task_related_dls3)';

                    task_related_m1 = [];
                    for time_step = 16:36
                        tmp_task_related_m1 = [];
                        for unit=1:length(m1_reach)
                            tmp_raster = [];
                            for trial = find(training_trials==1)
                                tmp_raster = [tmp_raster; m1_reach(unit).raster_100ms(trial,time_step-15:time_step)];
                            end
                            tmp_task_related_m1 = [tmp_task_related_m1 mean(tmp_raster)];
                        end
                        task_related_m1 = [task_related_m1; tmp_task_related_m1];
                    end
                    tbl = array2table([task_related_m1 task_related_dls1]);
                    mdl1 = fitlm(tbl);
                    tbl = array2table([task_related_m1 task_related_dls2]);
                    mdl2 = fitlm(tbl);
                    tbl = array2table([task_related_m1 task_related_dls3]);
                    mdl3 = fitlm(tbl);

                    %%% TESTING DATA

                    task_related_dls1 = [];
                    task_related_dls2 = [];
                    task_related_dls3 = [];
                    for trial = find(test_trials==1)
                        task_related_dls1 = [task_related_dls1; dls_pc1(trial,16:36)];
                        task_related_dls2 = [task_related_dls2; dls_pc2(trial,16:36)];
                        task_related_dls3 = [task_related_dls3; dls_pc3(trial,16:36)];
                    end
                    task_related_dls1 = mean(task_related_dls1)';
                    task_related_dls2 = mean(task_related_dls2)';
                    task_related_dls3 = mean(task_related_dls3)';

                    task_related_m1 = [];
                    for time_step = 16:36
                        tmp_task_related_m1 = [];
                        for unit=1:length(m1_reach)
                            tmp_raster = [];
                            for trial = find(test_trials==1)
                                tmp_raster = [tmp_raster; m1_reach(unit).raster_100ms(trial,time_step-15:time_step)];
                            end
                            tmp_task_related_m1 = [tmp_task_related_m1 mean(tmp_raster)];
                        end
                        task_related_m1 = [task_related_m1; tmp_task_related_m1];
                    end

                    ypred1 = predict(mdl1,task_related_m1);
                    ypred2 = predict(mdl2,task_related_m1);
                    ypred3 = predict(mdl3,task_related_m1);

                    [r p] = corrcoef(ypred1,task_related_dls1);
                    cv_corr1 = [cv_corr1; r(1,2)];

                    [r p] = corrcoef(ypred2,task_related_dls2);
                    cv_corr2 = [cv_corr2; r(1,2)];

                    [r p] = corrcoef(ypred3,task_related_dls3);
                    cv_corr3 = [cv_corr3; r(1,2)];

                end

                tmp_corr1 = [tmp_corr1; mean(cv_corr1)];
                tmp_corr2 = [tmp_corr2; mean(cv_corr2)];
                tmp_corr3 = [tmp_corr3; mean(cv_corr3)];

            else

                tmp_corr1 = [tmp_corr1; NaN];
                tmp_corr2 = [tmp_corr2; NaN];
                tmp_corr3 = [tmp_corr3; NaN];

            end

            num_units = [num_units; length(dls_reach) length(m1_reach)];

        end

        baseline_curves{animal} = [tmp_corr1'; tmp_corr2'; tmp_corr3'];
end

fig = figure; 

    subplot(2,1,1); hold on;
        tmp_line = [];
        for day = 1:8
            mean_day = [];
            for animal = 1:4
                if size(baseline_curves{animal},2)>=day
                    mean_day = [mean_day; mean(baseline_curves{animal}(:,day))];
                end
            end
            scatter(day*ones(1,sum(~isnan(mean_day))),mean_day(~isnan(mean_day)),30,[.5 .5 .5],'filled')
            errorbar(day,mean(mean_day(~isnan(mean_day))),std(mean_day(~isnan(mean_day)))/sqrt(length(mean_day(~isnan(mean_day)))),'color','k','LineWidth',2)
            tmp_line = [tmp_line mean(mean_day(~isnan(mean_day)))];
        end
        plot(tmp_line,'color','k','LineWidth',2)
        xlim([0 9])
        ylim([-0.2 0.7])

    subplot(2,1,2); hold on;
        tmp_line = [];
        for day = 1:8
            mean_day = [];
            for animal = 1:4
                if size(animal_curves{animal},2)>=day
                    mean_day = [mean_day; mean(animal_curves{animal}(:,day))];
                end
            end
            scatter(day*ones(1,sum(~isnan(mean_day))),mean_day(~isnan(mean_day)),30,[.5 .5 .5],'filled')
            errorbar(day,mean(mean_day(~isnan(mean_day))),std(mean_day(~isnan(mean_day)))/sqrt(length(mean_day(~isnan(mean_day)))),'color','k','LineWidth',2)
            tmp_line = [tmp_line mean(mean_day(~isnan(mean_day)))];
        end
        plot(tmp_line,'color','k','LineWidth',2)
        xlim([0 9])
        ylim([-0.2 0.7])

    %%% STATS

        early_reach = [];
        for n = 1:4
            early_reach = [early_reach mean(animal_curves{n}(:,1:3))];
        end
        early_reach = early_reach(~isnan(early_reach));

        late_reach = [];
        for n = 1:4
            late_reach = [late_reach mean(animal_curves{n}(:,end-2:end))];
        end
        late_reach = late_reach(~isnan(late_reach));

        early_baseline = [];
        for n = 1:4
            early_baseline = [early_baseline mean(baseline_curves{n}(:,1:3))];
        end
        early_baseline = early_baseline(~isnan(early_baseline));

        late_baseline = [];
        for n = 1:4
            late_baseline = [late_baseline mean(baseline_curves{n}(:,end-2:end))];
        end
        late_baseline = late_baseline(~isnan(late_baseline));

        length(early_reach)
        length(late_reach)

        mean(early_reach)
        std(early_reach)/sqrt(length(early_reach))

        mean(late_reach)
        std(late_reach)/sqrt(length(late_reach))

        [p h] = ranksum(early_reach(~isnan(early_reach)), late_reach(~isnan(late_reach)))

        length(early_baseline)
        length(late_baseline)

        mean(early_baseline)
        std(early_baseline)/sqrt(length(early_baseline))

        mean(late_baseline)
        std(late_baseline)/sqrt(length(late_baseline))

        [p h] = ranksum(early_baseline(~isnan(early_baseline)), late_baseline(~isnan(late_baseline)))  
        
%% FIG 4 C | CORRELATION BETWEEN ABILITY TO PREDICT DLS PCA TRAJ FROM M1 AND EACH DAY'S 4-8HZ LFP COHERENCE

%%% CALCULATE PREDICTION VALUES

    animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
    animal_days = {[1:8],[2:6],[2:10],[4:17],[4:13],[3:5 7:9]};

    all_animal_reach_pred = [];

    for animal = 1:length(animal_id)

        animal_path = [data_path '\' animal_id{animal}];

        tmp_corr1 = [];
        tmp_corr2 = [];
        tmp_corr3 = [];

        for day = animal_days{animal}

            load([animal_path '\reach_activity\day_' num2str(day) '_reach_mod.mat']);

            if length(dls_reach)>6 && length(m1_reach)>6

                %%% GENERATE PCs
                all_dls_raster = [];
                for unit = 1:length(dls_reach)
                    tmp_dls_raster = [];
                    for trial = 1:size(dls_reach(unit).raster_100ms,1)
                        tmp_dls_raster = [tmp_dls_raster dls_reach(unit).raster_100ms(trial,:)];
                    end
                    all_dls_raster = [all_dls_raster; tmp_dls_raster];
                end
                [coeff,score,latent] = pca(all_dls_raster');
                dls_pc1 = reshape(score(:,1),[size(dls_reach(unit).raster_100ms,2) size(dls_reach(unit).raster_100ms,1)])';
                dls_pc2 = reshape(score(:,2),[size(dls_reach(unit).raster_100ms,2) size(dls_reach(unit).raster_100ms,1)])';
                dls_pc3 = reshape(score(:,3),[size(dls_reach(unit).raster_100ms,2) size(dls_reach(unit).raster_100ms,1)])';

                %%% 5 FOLD CROSS VALIDATION
                trial_ind = 1:size(m1_reach(1).raster_100ms,1);
                cv_sets = [1 ...
                    round(size(m1_reach(1).raster_100ms,1)/5) ...
                    round(size(m1_reach(1).raster_100ms,1)*2/5) ...
                    round(size(m1_reach(1).raster_100ms,1)*3/5) ...
                    round(size(m1_reach(1).raster_100ms,1)*4/5) ...
                    size(m1_reach(1).raster_100ms,1)];

                cv_corr1 = [];
                cv_corr2 = [];
                cv_corr3 = [];

                for cv_n = 1:5

                    test_trials = zeros(1,size(m1_reach(1).raster_100ms,1));
                    training_trials = ones(1,size(m1_reach(1).raster_100ms,1));
                    test_trials(trial_ind(cv_sets(cv_n):cv_sets(cv_n+1)))=1;
                    training_trials(trial_ind(cv_sets(cv_n):cv_sets(cv_n+1)))=0;

                    %%% TRAINING DATA

                    task_related_dls1 = [];
                    task_related_dls2 = [];
                    task_related_dls3 = [];
                    for trial = find(training_trials==1)
                        task_related_dls1 = [task_related_dls1; dls_pc1(trial,41:60)];
                        task_related_dls2 = [task_related_dls2; dls_pc2(trial,41:60)];
                        task_related_dls3 = [task_related_dls3; dls_pc3(trial,41:60)];
                    end
                    task_related_dls1 = mean(task_related_dls1)';
                    task_related_dls2 = mean(task_related_dls2)';
                    task_related_dls3 = mean(task_related_dls3)';

                    task_related_m1 = [];
                    for time_step = 41:60
                        tmp_task_related_m1 = [];
                        for unit=1:length(m1_reach)
                            tmp_raster = [];
                            for trial = find(training_trials==1)
                                tmp_raster = [tmp_raster; m1_reach(unit).raster_100ms(trial,time_step-15:time_step)];
                            end
                            tmp_task_related_m1 = [tmp_task_related_m1 mean(tmp_raster)];
                        end
                        task_related_m1 = [task_related_m1; tmp_task_related_m1];
                    end
                    tbl = array2table([task_related_m1 task_related_dls1]);
                    mdl1 = fitlm(tbl);
                    tbl = array2table([task_related_m1 task_related_dls2]);
                    mdl2 = fitlm(tbl);
                    tbl = array2table([task_related_m1 task_related_dls3]);
                    mdl3 = fitlm(tbl);

                    %%% TESTING DATA

                    task_related_dls1 = [];
                    task_related_dls2 = [];
                    task_related_dls3 = [];
                    for trial = find(test_trials==1)
                        task_related_dls1 = [task_related_dls1; dls_pc1(trial,41:60)];
                        task_related_dls2 = [task_related_dls2; dls_pc2(trial,41:60)];
                        task_related_dls3 = [task_related_dls3; dls_pc3(trial,41:60)];
                    end
                    task_related_dls1 = mean(task_related_dls1)';
                    task_related_dls2 = mean(task_related_dls2)';
                    task_related_dls3 = mean(task_related_dls3)';

                    task_related_m1 = [];
                    for time_step = 41:60
                        tmp_task_related_m1 = [];
                        for unit=1:length(m1_reach)
                            tmp_raster = [];
                            for trial = find(test_trials==1)
                                tmp_raster = [tmp_raster; m1_reach(unit).raster_100ms(trial,time_step-15:time_step)];
                            end
                            tmp_task_related_m1 = [tmp_task_related_m1 mean(tmp_raster)];
                        end
                        task_related_m1 = [task_related_m1; tmp_task_related_m1];
                    end

                    ypred1 = predict(mdl1,task_related_m1);
                    ypred2 = predict(mdl2,task_related_m1);
                    ypred3 = predict(mdl3,task_related_m1);

                    [r p] = corrcoef(ypred1,task_related_dls1);
                    cv_corr1 = [cv_corr1; r(1,2)];

                    [r p] = corrcoef(ypred2,task_related_dls2);
                    cv_corr2 = [cv_corr2; r(1,2)];

                    [r p] = corrcoef(ypred3,task_related_dls3);
                    cv_corr3 = [cv_corr3; r(1,2)];

                end

                tmp_corr1 = [tmp_corr1; mean(cv_corr1)];
                tmp_corr2 = [tmp_corr2; mean(cv_corr2)];
                tmp_corr3 = [tmp_corr3; mean(cv_corr3)];

            else

                tmp_corr1 = [tmp_corr1; NaN];
                tmp_corr2 = [tmp_corr2; NaN];
                tmp_corr3 = [tmp_corr3; NaN];

            end
        end

        tmp_corr = mean([tmp_corr1'; tmp_corr1'; tmp_corr1']);
        non_nan_index = find(~isnan(tmp_corr));
        tmp_corr(non_nan_index) = zscore(tmp_corr(~isnan(tmp_corr)));
        all_animal_reach_pred = [all_animal_reach_pred tmp_corr];

    end

    all_animal_baseline_pred = [];

    for animal = 1:length(animal_id)

        animal_path = [data_path '\' animal_id{animal}];

        tmp_corr1 = [];
        tmp_corr2 = [];
        tmp_corr3 = [];

        for day = animal_days{animal}

            load([animal_path '\reach_activity\day_' num2str(day) '_reach_mod.mat']);

            if length(dls_reach)>6 && length(m1_reach)>6

                %%% GENERATE PCs
                all_dls_raster = [];
                for unit = 1:length(dls_reach)
                    tmp_dls_raster = [];
                    for trial = 1:size(dls_reach(unit).raster_100ms,1)
                        tmp_dls_raster = [tmp_dls_raster dls_reach(unit).raster_100ms(trial,:)];
                    end
                    all_dls_raster = [all_dls_raster; tmp_dls_raster];
                end
                [coeff,score,latent] = pca(all_dls_raster');
                dls_pc1 = reshape(score(:,1),[size(dls_reach(unit).raster_100ms,2) size(dls_reach(unit).raster_100ms,1)])';
                dls_pc2 = reshape(score(:,2),[size(dls_reach(unit).raster_100ms,2) size(dls_reach(unit).raster_100ms,1)])';
                dls_pc3 = reshape(score(:,3),[size(dls_reach(unit).raster_100ms,2) size(dls_reach(unit).raster_100ms,1)])';

                %%% 5 FOLD CROSS VALIDATION
                trial_ind = 1:size(m1_reach(1).raster_100ms,1);
                cv_sets = [1 ...
                    round(size(m1_reach(1).raster_100ms,1)/5) ...
                    round(size(m1_reach(1).raster_100ms,1)*2/5) ...
                    round(size(m1_reach(1).raster_100ms,1)*3/5) ...
                    round(size(m1_reach(1).raster_100ms,1)*4/5) ...
                    size(m1_reach(1).raster_100ms,1)];

                cv_corr1 = [];
                cv_corr2 = [];
                cv_corr3 = [];

                for cv_n = 1:5

                    test_trials = zeros(1,size(m1_reach(1).raster_100ms,1));
                    training_trials = ones(1,size(m1_reach(1).raster_100ms,1));
                    test_trials(trial_ind(cv_sets(cv_n):cv_sets(cv_n+1)))=1;
                    training_trials(trial_ind(cv_sets(cv_n):cv_sets(cv_n+1)))=0;

                    %%% TRAINING DATA

                    task_related_dls1 = [];
                    task_related_dls2 = [];
                    task_related_dls3 = [];
                    for trial = find(training_trials==1)
                        task_related_dls1 = [task_related_dls1; dls_pc1(trial,16:36)];
                        task_related_dls2 = [task_related_dls2; dls_pc2(trial,16:36)];
                        task_related_dls3 = [task_related_dls3; dls_pc3(trial,16:36)];
                    end
                    task_related_dls1 = mean(task_related_dls1)';
                    task_related_dls2 = mean(task_related_dls2)';
                    task_related_dls3 = mean(task_related_dls3)';

                    task_related_m1 = [];
                    for time_step = 16:36
                        tmp_task_related_m1 = [];
                        for unit=1:length(m1_reach)
                            tmp_raster = [];
                            for trial = find(training_trials==1)
                                tmp_raster = [tmp_raster; m1_reach(unit).raster_100ms(trial,time_step-15:time_step)];
                            end
                            tmp_task_related_m1 = [tmp_task_related_m1 mean(tmp_raster)];
                        end
                        task_related_m1 = [task_related_m1; tmp_task_related_m1];
                    end
                    tbl = array2table([task_related_m1 task_related_dls1]);
                    mdl1 = fitlm(tbl);
                    tbl = array2table([task_related_m1 task_related_dls2]);
                    mdl2 = fitlm(tbl);
                    tbl = array2table([task_related_m1 task_related_dls3]);
                    mdl3 = fitlm(tbl);

                    %%% TESTING DATA

                    task_related_dls1 = [];
                    task_related_dls2 = [];
                    task_related_dls3 = [];
                    for trial = find(test_trials==1)
                        task_related_dls1 = [task_related_dls1; dls_pc1(trial,16:36)];
                        task_related_dls2 = [task_related_dls2; dls_pc2(trial,16:36)];
                        task_related_dls3 = [task_related_dls3; dls_pc3(trial,16:36)];
                    end
                    task_related_dls1 = mean(task_related_dls1)';
                    task_related_dls2 = mean(task_related_dls2)';
                    task_related_dls3 = mean(task_related_dls3)';

                    task_related_m1 = [];
                    for time_step = 16:36
                        tmp_task_related_m1 = [];
                        for unit=1:length(m1_reach)
                            tmp_raster = [];
                            for trial = find(test_trials==1)
                                tmp_raster = [tmp_raster; m1_reach(unit).raster_100ms(trial,time_step-15:time_step)];
                            end
                            tmp_task_related_m1 = [tmp_task_related_m1 mean(tmp_raster)];
                        end
                        task_related_m1 = [task_related_m1; tmp_task_related_m1];
                    end

                    ypred1 = predict(mdl1,task_related_m1);
                    ypred2 = predict(mdl2,task_related_m1);
                    ypred3 = predict(mdl3,task_related_m1);

                    [r p] = corrcoef(ypred1,task_related_dls1);
                    cv_corr1 = [cv_corr1; r(1,2)];

                    [r p] = corrcoef(ypred2,task_related_dls2);
                    cv_corr2 = [cv_corr2; r(1,2)];

                    [r p] = corrcoef(ypred3,task_related_dls3);
                    cv_corr3 = [cv_corr3; r(1,2)];

                end

                tmp_corr1 = [tmp_corr1; mean(cv_corr1)];
                tmp_corr2 = [tmp_corr2; mean(cv_corr2)];
                tmp_corr3 = [tmp_corr3; mean(cv_corr3)];

            else

                tmp_corr1 = [tmp_corr1; NaN];
                tmp_corr2 = [tmp_corr2; NaN];
                tmp_corr3 = [tmp_corr3; NaN];

            end
        end

        tmp_corr = mean([tmp_corr1'; tmp_corr1'; tmp_corr1']);
        non_nan_index = find(~isnan(tmp_corr));
        tmp_corr(non_nan_index) = zscore(tmp_corr(~isnan(tmp_corr)));
        all_animal_baseline_pred = [all_animal_baseline_pred tmp_corr];

    end

%%% LFP COHERENCE

    all_animal_coh = [];

    %%% LC3

        load([data_path '\LC3\sleep_activity\LFP_coherence.mat'])
        increase_channels = [];
        for pairs=1:256
            if mean(mean([pre_LFP_coh{pairs}(8,48:113); post_LFP_coh{pairs}(8,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)])))
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
        all_animal_coh = [all_animal_coh zscore(mean(day_pre))];

    %%% LC4

        load([data_path '\LC4\sleep_activity\LFP_coherence.mat'])
        increase_channels = [];
        for pairs=1:256
            if mean(mean([pre_LFP_coh{pairs}(6,48:113); post_LFP_coh{pairs}(6,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)])))
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
        all_animal_coh = [all_animal_coh zscore(mean(day_pre))];

    %%% LC6

        load([data_path '\LC6\sleep_activity\LFP_coherence.mat'])
        increase_channels = [];
        for pairs=1:256
            if mean(mean([pre_LFP_coh{pairs}(9,48:113); post_LFP_coh{pairs}(9,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)])))
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
        all_animal_coh = [all_animal_coh zscore(mean(day_pre))];

    %%% LC5

        load([data_path '\LC5\sleep_activity\LFP_coherence.mat'])
        increase_channels = [];
        for pairs=1:512
            if mean(mean([pre_LFP_coh{pairs}(15,48:113); post_LFP_coh{pairs}(15,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(2,48:113); post_LFP_coh{pairs}(2,48:113)])))
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
            for day=2:15
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,48:113))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,48:113))];
            end
            day_pre = [day_pre; tmp_pre];
            day_post = [day_post; tmp_post];
        end
        all_animal_coh = [all_animal_coh zscore(mean(day_pre))];

    %%% LC2

        load([data_path '\LC2\sleep_activity\LFP_coherence.mat'])
        increase_channels = [];
        for pairs=1:1024
            if mean(mean([pre_LFP_coh{pairs}(10,48:113); post_LFP_coh{pairs}(10,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)])))
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
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,48:113))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,48:113))];
            end
            day_pre = [day_pre; tmp_pre];
            day_post = [day_post; tmp_post];
        end
        all_animal_coh = [all_animal_coh zscore(mean([mean(day_pre); mean(day_post)]))];

    %%% LC1

        load([data_path '\LC1\sleep_activity\LFP_coherence.mat'])
        increase_channels = [];
        for pairs=1:1406
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
            for day=3:8
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,48:113))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,48:113))];
            end
            day_pre = [day_pre; tmp_pre];
            day_post = [day_post; tmp_post];
        end
        all_animal_coh = [all_animal_coh zscore(mean([mean(day_pre); mean(day_post)]))];

%%% PLOT PRED/COH SCATTER

fig = figure;
    hold on;
    scatter(all_animal_coh(~isnan(all_animal_baseline_pred)),all_animal_baseline_pred(~isnan(all_animal_baseline_pred)),30,[0 0 0],'filled')
    [R,P] = corrcoef(all_animal_coh(~isnan(all_animal_baseline_pred)),all_animal_baseline_pred(~isnan(all_animal_baseline_pred)))
    [p] = polyfit(all_animal_coh(~isnan(all_animal_baseline_pred)),all_animal_baseline_pred(~isnan(all_animal_baseline_pred)),1)
    x=[-3:.1:3]
    xlim([-2 3])
    ylim([-3 3])
    y = x*p(1)+p(2)
    plot(x,y,'color','r','LineStyle','--','LineWidth',3)
    title(['baseline pred vs. lfp coh - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);

fig = figure;
    hold on;
    scatter(all_animal_coh(~isnan(all_animal_reach_pred)),all_animal_reach_pred(~isnan(all_animal_reach_pred)),30,[0 0 0],'filled')
    [R,P] = corrcoef(all_animal_coh(~isnan(all_animal_reach_pred)),all_animal_reach_pred(~isnan(all_animal_reach_pred)))
    [p] = polyfit(all_animal_coh(~isnan(all_animal_reach_pred)),all_animal_reach_pred(~isnan(all_animal_reach_pred)),1)
    x=[-3:.1:3]
    xlim([-2 3])
    ylim([-3 3])
    y = x*p(1)+p(2)
    plot(x,y,'color','r','LineStyle','--','LineWidth',3)
    title(['reach pred vs. lfp coh - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2)) ' days: ' num2str(length(all_animal_coh(~isnan(all_animal_reach_pred))))]);
            
%% sFIG 1 | REACH MODULATION EARLY VS. LATE

    %%% FIRST THREE DAYS

        animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
        animal_days = {[1:3],[2:4],[2:4],[4:6],[4:6],[3:5]};

        all_early_m1 = [];
        all_early_m1_max = [];
        all_early_dls = [];
        all_early_dls_max = [];
        
        all_early_m1_mod = [];
        all_early_dls_mod = [];

        for animal = 1:length(animal_id)

            animal_path = [data_path '\' animal_id{animal}];

            for day = animal_days{animal}

                load([animal_path '\reach_activity\day_' num2str(day) '_reach_mod.mat']);

                for n = 1:length(m1_reach)
                    all_early_m1 = [all_early_m1; smooth(zscore(mean(m1_reach(n).raster)),3)'];
                    [~, i] = max(smooth(zscore(mean(m1_reach(n).raster)),3)');
                    all_early_m1_max = [all_early_m1_max i];
                    
%                     tmp_mod = abs(mean(m1_reach(n).raster)-mean(mean(m1_reach(n).raster)));
                    tmp_mod = abs(zscore(mean(m1_reach(n).raster)));
                    all_early_m1_mod = [all_early_m1_mod sum(tmp_mod(160:220))/sum(tmp_mod(60:120))];
                end

                for n = 1:length(dls_reach)
                    all_early_dls = [all_early_dls; smooth(zscore(mean(dls_reach(n).raster)),3)'];
                    [~, i] = max(smooth(zscore(mean(dls_reach(n).raster)),3)');
                    all_early_dls_max = [all_early_dls_max i];
                    
%                     tmp_mod = abs(mean(dls_reach(n).raster)-mean(mean(dls_reach(n).raster)));
                    tmp_mod = abs(zscore(mean(dls_reach(n).raster)));                    
                    all_early_dls_mod = [all_early_dls_mod sum(tmp_mod(160:220))/sum(tmp_mod(60:120))];
                end

            end
        end

        [~,all_early_m1_max] = sort(all_early_m1_max);

        fig = figure;
        hold on;
            imagesc(all_early_m1(all_early_m1_max,:),[-3 8])
            xlim([100 300])
            ylim([.5 .5+size(all_early_m1_max,2)])
            title('EARLY M1 REACH MODULATION')
            colorbar
            
        [~,all_early_dls_max] = sort(all_early_dls_max);

        fig = figure;
        hold on;
            imagesc(all_early_dls(all_early_dls_max,:),[-3 8])
            xlim([100 300])
            ylim([.5 .5+size(all_early_dls_max,2)])
            title('EARLY DLS REACH MODULATION')
            colorbar

    %%% LAST THREE DAYS

        animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
        animal_days = {[6:8],[4:6],[8:10],[15:17],[11:13],[7:9]};

        all_late_m1 = [];
        all_late_m1_max = [];
        all_late_dls = [];
        all_late_dls_max = [];
        
        all_late_m1_mod = [];
        all_late_dls_mod = [];
        
        for animal = 1:length(animal_id)

            animal_path = [data_path '\' animal_id{animal}];

            for day = animal_days{animal}

                load([animal_path '\reach_activity\day_' num2str(day) '_reach_mod.mat']);

                for n = 1:length(m1_reach)
                    all_late_m1 = [all_late_m1; smooth(zscore(mean(m1_reach(n).raster)),3)'];
                    [~, i] = max(smooth(zscore(mean(m1_reach(n).raster)),3)');
                    all_late_m1_max = [all_late_m1_max i];
                    
%                     tmp_mod = abs(mean(m1_reach(n).raster)-mean(mean(m1_reach(n).raster)));
                    tmp_mod = abs(zscore(mean(m1_reach(n).raster)));
                    all_late_m1_mod = [all_late_m1_mod sum(tmp_mod(160:220))/sum(tmp_mod(60:120))];
                end

                for n = 1:length(dls_reach)
                    all_late_dls = [all_late_dls; smooth(zscore(mean(dls_reach(n).raster)),3)'];
                    [~, i] = max(smooth(zscore(mean(dls_reach(n).raster)),3)');
                    all_late_dls_max = [all_late_dls_max i];
                    
%                     tmp_mod = abs(mean(dls_reach(n).raster)-mean(mean(dls_reach(n).raster)));
                    tmp_mod = abs(zscore(mean(dls_reach(n).raster)));  
                    all_late_dls_mod = [all_late_dls_mod sum(tmp_mod(160:220))/sum(tmp_mod(60:120))];
                end

            end
        end

        [~,all_late_m1_max] = sort(all_late_m1_max);

        fig = figure;
        hold on;
            imagesc(all_late_m1(all_late_m1_max,:),[-3 8])
            xlim([100 300])
            ylim([.5 .5+size(all_late_m1_max,2)])
            title('LATE M1 REACH MODULATION')
            colorbar

        [~,all_late_dls_max] = sort(all_late_dls_max);

        fig = figure;
        hold on;
            imagesc(all_late_dls(all_late_dls_max,:),[-3 8])
            xlim([100 300])
            ylim([.5 .5+size(all_late_dls_max,2)])
            title('LATE DLS REACH MODULATION')
            colorbar

        [h,p,ci,stats] = ttest2(all_early_m1_mod,all_late_m1_mod)

        mean(all_early_m1_mod)
        std(all_early_m1_mod)/sqrt(length(all_early_m1_mod))
        mean(all_late_m1_mod)
        std(all_late_m1_mod)/sqrt(length(all_late_m1_mod))

        [h,p,ci,stats] = ttest2(all_early_dls_mod,all_late_dls_mod)
    
        mean(all_early_dls_mod)
        std(all_early_dls_mod)/sqrt(length(all_early_dls_mod))
        mean(all_late_dls_mod)
        std(all_late_dls_mod)/sqrt(length(all_late_dls_mod))

%% COMPUTE VARIANCE EXPLAINED BY FIRST TWO PC ON FIRST/LAST THREE DAYS

animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
animal_days = {[1:8],[2:6],[2:10],[4:17],[4:13],[3:5 7:9]};

early_VE = [];
late_VE = [];
early_num_units = [];
late_num_units = [];

for animal = 1:length(animal_id)
    
    animal_path = [data_path '\' animal_id{animal}];
    
    animal_explained = [];
    num_DLS_units =[];
    
    for day = animal_days{animal}
        load([animal_path '\reach_activity\day_' num2str(day) '_reach_mod.mat']);
        
        all_dls_raster = [];
        for unit = 1:length(dls_reach)
            tmp_dls_raster = [];
            for trial = 1:size(dls_reach(unit).raster_100ms,1)
                tmp_dls_raster = [tmp_dls_raster dls_reach(unit).raster_100ms(trial,:)];
            end
            all_dls_raster = [all_dls_raster; tmp_dls_raster];
        end
        
        if length(dls_reach)>6 && length(m1_reach)>6
            [coeff,score,latent,~,explained] = pca(all_dls_raster');
            animal_explained = [animal_explained sum(explained(1:3))];   
            num_DLS_units = [num_DLS_units length(dls_reach)];            
        else
            animal_explained = [animal_explained NaN];   
            num_DLS_units = [num_DLS_units NaN];       
        end
            
    end
    
    early_VE = [early_VE animal_explained(1:3)];
    late_VE = [late_VE animal_explained(end-2:end)];
    early_num_units = [early_num_units num_DLS_units(1:3)];
    late_num_units = [late_num_units num_DLS_units(end-2:end)];

end

length(early_VE(~isnan(early_VE)))
mean(early_VE(~isnan(early_VE)))
std(early_VE(~isnan(early_VE)))/sqrt(length(early_VE(~isnan(early_VE))))

length(early_VE(~isnan(late_VE)))
mean(late_VE(~isnan(late_VE)))
std(late_VE(~isnan(late_VE)))/sqrt(length(late_VE(~isnan(late_VE))))

[h p ci stats] = ttest2(early_VE(~isnan(early_VE)),late_VE(~isnan(late_VE)))

mean(early_num_units(~isnan(early_num_units)))
std(early_num_units(~isnan(early_num_units)))/sqrt(length(early_num_units(~isnan(early_num_units))))
mean(late_num_units(~isnan(late_num_units)))
std(late_num_units(~isnan(late_num_units)))/sqrt(length(late_num_units(~isnan(late_num_units))))
[h p ci stats] = ttest2(early_num_units(~isnan(early_num_units)),late_num_units(~isnan(late_num_units)))






    


