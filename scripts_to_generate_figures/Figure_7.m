%% SCRIPT TO GENERATE FIGURE 7 & FIGURE 7 FIGURE SUPPLEMENTS
% 
% to generate figures, download data from https://doi.org/10.7272/Q6KK9927

clc; clear all; close all;
data_path = 'D:\Lemke_etal_2021\Lemke_eLife_2021_data';

%% FIG 7 A | PRE/POST SPINDLE NESTING

    animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
    animal_days = {[1:8],[1:6],[2:10],[4:17],[4:13],[1:5 7:9]};

    pre_SO_nesting = [];
    pre_delta_nesting = [];
    post_SO_nesting = [];
    post_delta_nesting = [];

    for animal = 1:length(animal_id)

        CC_path = [data_path '\' animal_id{animal} '\cross_correlations'];

        for day = animal_days{animal}

            load([data_path '\' animal_id{animal} '\sleep_activity\day_'  num2str(day) '_block_1_sleep_activity.mat']);

            for ind = 1:length(sleep_rhythms.spindles{1,1}.pks)
                tmp_diff_SO = sleep_rhythms.spindles{1,1}.pks(ind)-sleep_rhythms.so_delta.SO_zc;
                tmp_diff_delta = sleep_rhythms.spindles{1,1}.pks(ind)-sleep_rhythms.so_delta.delta_zc;
                if ~isempty(tmp_diff_SO(tmp_diff_SO>0)) & ~isempty(tmp_diff_delta(tmp_diff_delta>0))
                    pre_SO_nesting = [pre_SO_nesting min(tmp_diff_SO(tmp_diff_SO>0))];
                    pre_delta_nesting = [pre_delta_nesting min(tmp_diff_delta(tmp_diff_delta>0))];
                end
            end

            if animal == 1 && day == 6
                load([data_path '\' animal_id{animal} '\sleep_activity\day_'  num2str(day) '_block_3_sleep_activity.mat']);
            else
                load([data_path '\' animal_id{animal} '\sleep_activity\day_'  num2str(day) '_block_2_sleep_activity.mat']);
            end

            for ind = 1:length(sleep_rhythms.spindles{1,1}.pks)
                tmp_diff_SO = sleep_rhythms.spindles{1,1}.pks(ind)-sleep_rhythms.so_delta.SO_zc;
                tmp_diff_delta = sleep_rhythms.spindles{1,1}.pks(ind)-sleep_rhythms.so_delta.delta_zc;
                if ~isempty(tmp_diff_SO(tmp_diff_SO>0)) & ~isempty(tmp_diff_delta(tmp_diff_delta>0))
                    post_SO_nesting = [post_SO_nesting min(tmp_diff_SO(tmp_diff_SO>0))];
                    post_delta_nesting = [post_delta_nesting min(tmp_diff_delta(tmp_diff_delta>0))];
                end
            end

        end
    end

    figure; hold on
        histogram(pre_SO_nesting,[0:0.5:10],'normalization','probability','DisplayStyle','stairs')
        histogram(post_SO_nesting,[0:0.5:10],'normalization','probability','DisplayStyle','stairs')
        [h p] = kstest2(pre_SO_nesting,post_SO_nesting)
        title(num2str(p));
        ylim([0 0.1])

    figure; hold on
        histogram(pre_delta_nesting,[0:0.5:10],'normalization','probability','DisplayStyle','stairs')
        histogram(post_delta_nesting,[0:0.5:10],'normalization','probability','DisplayStyle','stairs')
        [h p] = kstest2(pre_delta_nesting,post_delta_nesting)
        title(num2str(p));
   
%% FIG 7 B | sFIG 1 | CC CHANGE BEFORE/AFTER NESTED OR NON-NESTED SPINDLES AND M1/DLS FIRING RATES BEFORE/AFTER NESTED OR NON-NESTED SPINDLES

    animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
    animal_days = {[1:8],[1:6],[2:10],[4:17],[4:13],[1:5 7:9]};

    all_pre_spindle_cc = [];
    all_post_spindle_cc = [];

    pre_SO_nested_spindle_cc = [];
    pre_SO_not_nested_spindle_cc = [];
    pre_SO_nested_spindle_raster_m1 = [];
    pre_SO_not_nested_spindle_raster_m1 = [];
    pre_SO_nested_spindle_raster_dls = []; 
    pre_SO_not_nested_spindle_raster_dls = [];

    pre_delta_nested_spindle_cc = [];
    pre_delta_not_nested_spindle_cc = [];
    pre_delta_nested_spindle_raster_m1 = [];
    pre_delta_not_nested_spindle_raster_m1 = [];
    pre_delta_nested_spindle_raster_dls = [];
    pre_delta_not_nested_spindle_raster_dls = [];

    pre_not_nested_spindle_cc = [];

    post_SO_nested_spindle_cc = [];
    post_SO_not_nested_spindle_cc = [];
    post_SO_nested_spindle_raster_m1 = [];
    post_SO_not_nested_spindle_raster_m1 = [];
    post_SO_nested_spindle_raster_dls = [];
    post_SO_not_nested_spindle_raster_dls = [];

    post_delta_nested_spindle_cc = [];
    post_delta_not_nested_spindle_cc = [];
    post_delta_nested_spindle_raster_m1 = [];
    post_delta_not_nested_spindle_raster_m1 = [];
    post_delta_nested_spindle_raster_dls = [];
    post_delta_not_nested_spindle_raster_dls = [];

    post_not_nested_spindle_cc = [];

    for animal = 1:length(animal_id)

        CC_path = [data_path '\' animal_id{animal} '\cross_correlations'];

        for day = animal_days{animal}

            load([CC_path '\day_' num2str(day) '\spindle_cc_before_after.mat']);

            if isfield(all_pairs,'pre_spindle_cc')

                load([data_path '\' animal_id{animal} '\sleep_activity\day_'  num2str(day) '_block_1_sleep_activity.mat']);

                for pair = 1:length(all_pairs)
                    if ~isempty(all_pairs(pair).pre_spindle_cc)
                        all_pre_spindle_cc = [all_pre_spindle_cc; mean(all_pairs(pair).pre_spindle_cc,1)];
                        all_post_spindle_cc = [all_post_spindle_cc; mean(all_pairs(pair).post_spindle_cc,1)];
                    end
                end

                SO_nesting = [];
                delta_nesting = [];
                for ind = 1:length(sleep_rhythms.spindles{1,1}.pks)
                    tmp_diff_SO = sleep_rhythms.spindles{1,1}.pks(ind)-sleep_rhythms.so_delta.SO_zc;
                    tmp_diff_delta = sleep_rhythms.spindles{1,1}.pks(ind)-sleep_rhythms.so_delta.delta_zc;
                    if ~isempty(tmp_diff_SO(tmp_diff_SO>0)) & ~isempty(tmp_diff_delta(tmp_diff_delta>0))
                        SO_nesting = [SO_nesting min(tmp_diff_SO(tmp_diff_SO>0))];
                        delta_nesting = [delta_nesting min(tmp_diff_delta(tmp_diff_delta>0))];
                    end
                end

                for pair = 1:length(all_pairs)
                    if ~isempty(all_pairs(pair).pre_spindle_cc)

                        pre_SO_nested_spindle_cc = [pre_SO_nested_spindle_cc; mean(all_pairs(pair).pre_spindle_cc(SO_nesting<1,:),1)];                    
                        pre_delta_nested_spindle_cc = [pre_delta_nested_spindle_cc; mean(all_pairs(pair).pre_spindle_cc(delta_nesting<1,:),1)];

                        pre_SO_not_nested_spindle_cc = [pre_SO_not_nested_spindle_cc; mean(all_pairs(pair).pre_spindle_cc(SO_nesting>5,:),1)];
                        pre_delta_not_nested_spindle_cc = [pre_delta_not_nested_spindle_cc; mean(all_pairs(pair).pre_spindle_cc(delta_nesting>5,:),1)];

                        pre_not_nested_spindle_cc = [pre_not_nested_spindle_cc; mean(all_pairs(pair).pre_spindle_cc(SO_nesting>5 & delta_nesting>5,:),1)];

                        pre_SO_nested_spindle_raster_m1 = [pre_SO_nested_spindle_raster_m1; mean(all_pairs(pair).pre_spindle_raster_m1(SO_nesting<1,:),1)];
                        pre_SO_nested_spindle_raster_dls = [pre_SO_nested_spindle_raster_dls; mean(all_pairs(pair).pre_spindle_raster_dls(SO_nesting<1,:),1)];
                        pre_SO_not_nested_spindle_raster_m1 = [pre_SO_not_nested_spindle_raster_m1; mean(all_pairs(pair).pre_spindle_raster_m1(SO_nesting>5,:),1)];
                        pre_SO_not_nested_spindle_raster_dls = [pre_SO_not_nested_spindle_raster_dls; mean(all_pairs(pair).pre_spindle_raster_dls(SO_nesting>5,:),1)];

                        pre_delta_nested_spindle_raster_m1 = [pre_delta_nested_spindle_raster_m1; mean(all_pairs(pair).pre_spindle_raster_m1(delta_nesting<1,:),1)];
                        pre_delta_nested_spindle_raster_dls = [pre_delta_nested_spindle_raster_dls; mean(all_pairs(pair).pre_spindle_raster_dls(delta_nesting<1,:),1)];
                        pre_delta_not_nested_spindle_raster_m1 = [pre_delta_not_nested_spindle_raster_m1; mean(all_pairs(pair).pre_spindle_raster_m1(delta_nesting>5,:),1)];
                        pre_delta_not_nested_spindle_raster_dls = [pre_delta_not_nested_spindle_raster_dls; mean(all_pairs(pair).pre_spindle_raster_dls(delta_nesting>5,:),1)];

                    end
                end

                if animal == 1 && day == 6
                    load([data_path '\' animal_id{animal} '\sleep_activity\day_'  num2str(day) '_block_3_sleep_activity.mat']);
                else
                    load([data_path '\' animal_id{animal} '\sleep_activity\day_'  num2str(day) '_block_2_sleep_activity.mat']);
                end

                SO_nesting = [];
                delta_nesting = [];
                for ind = 1:length(sleep_rhythms.spindles{1,1}.pks)
                    tmp_diff_SO = sleep_rhythms.spindles{1,1}.pks(ind)-sleep_rhythms.so_delta.SO_zc;
                    tmp_diff_delta = sleep_rhythms.spindles{1,1}.pks(ind)-sleep_rhythms.so_delta.delta_zc;
                    if ~isempty(tmp_diff_SO(tmp_diff_SO>0)) & ~isempty(tmp_diff_delta(tmp_diff_delta>0))
                        SO_nesting = [SO_nesting min(tmp_diff_SO(tmp_diff_SO>0))];
                        delta_nesting = [delta_nesting min(tmp_diff_delta(tmp_diff_delta>0))];
                    end
                end

                for pair = 1:length(all_pairs)
                    if ~isempty(all_pairs(pair).post_spindle_cc)

                        post_SO_nested_spindle_cc = [post_SO_nested_spindle_cc; mean(all_pairs(pair).post_spindle_cc(SO_nesting<1,:),1)];                    
                        post_delta_nested_spindle_cc = [post_delta_nested_spindle_cc; mean(all_pairs(pair).post_spindle_cc(delta_nesting<1,:),1)];

                        post_SO_not_nested_spindle_cc = [post_SO_not_nested_spindle_cc; mean(all_pairs(pair).post_spindle_cc(SO_nesting>5,:),1)];
                        post_delta_not_nested_spindle_cc = [post_delta_not_nested_spindle_cc; mean(all_pairs(pair).post_spindle_cc(delta_nesting>5,:),1)];

                        post_not_nested_spindle_cc = [post_not_nested_spindle_cc; mean(all_pairs(pair).post_spindle_cc(SO_nesting>5 & delta_nesting>5,:),1)];

                        post_SO_nested_spindle_raster_m1 = [post_SO_nested_spindle_raster_m1; mean(all_pairs(pair).post_spindle_raster_m1(SO_nesting<1,:),1)];
                        post_SO_nested_spindle_raster_dls = [post_SO_nested_spindle_raster_dls; mean(all_pairs(pair).post_spindle_raster_dls(SO_nesting<1,:),1)];
                        post_SO_not_nested_spindle_raster_m1 = [post_SO_not_nested_spindle_raster_m1; mean(all_pairs(pair).post_spindle_raster_m1(SO_nesting>5,:),1)];
                        post_SO_not_nested_spindle_raster_dls = [post_SO_not_nested_spindle_raster_dls; mean(all_pairs(pair).post_spindle_raster_dls(SO_nesting>5,:),1)];

                        post_delta_nested_spindle_raster_m1 = [post_delta_nested_spindle_raster_m1; mean(all_pairs(pair).post_spindle_raster_m1(delta_nesting<1,:),1)];
                        post_delta_nested_spindle_raster_dls = [post_delta_nested_spindle_raster_dls; mean(all_pairs(pair).post_spindle_raster_dls(delta_nesting<1,:),1)];
                        post_delta_not_nested_spindle_raster_m1 = [post_delta_not_nested_spindle_raster_m1; mean(all_pairs(pair).post_spindle_raster_m1(delta_nesting>5,:),1)];
                        post_delta_not_nested_spindle_raster_dls = [post_delta_not_nested_spindle_raster_dls; mean(all_pairs(pair).post_spindle_raster_dls(delta_nesting>5,:),1)];

                    end
                end
            end
        end

    end

    all_SO_nested_spindle = [pre_SO_nested_spindle_cc; post_SO_nested_spindle_cc];
    all_not_nested_spindle = [pre_not_nested_spindle_cc; post_SO_not_nested_spindle_cc];

    idx = ~isnan(mean(all_SO_nested_spindle')) & ~isnan(mean(all_not_nested_spindle'));
    norm_all_SO_nested_spindle_cc = all_SO_nested_spindle(idx,:)-repmat(mean(all_SO_nested_spindle(idx,1),2),1,6);
    norm_all_not_nested_spindle_cc = all_not_nested_spindle(idx,:)-repmat(mean(all_not_nested_spindle(idx,1),2),1,6);

    figure;

        subplot(2,1,1)
            hold on;
            for bin = 1:6
                errorbar(bin,mean(norm_all_not_nested_spindle_cc(:,bin)),std(norm_all_not_nested_spindle_cc(:,bin))/sqrt(size(norm_all_not_nested_spindle_cc,1)),'color','k')
                errorbar(bin,mean(norm_all_SO_nested_spindle_cc(:,bin)),std(norm_all_SO_nested_spindle_cc(:,bin))/sqrt(size(norm_all_SO_nested_spindle_cc,1)),'color','r')
                [h p] = ttest(norm_all_SO_nested_spindle_cc(:,bin),norm_all_not_nested_spindle_cc(:,bin))
                if p<0.05
                    scatter(bin,3e-4,30,[0 0 0],'filled')
                end
            end
            xlim([0 7])
            ylim([-2e-4 3.5e-4])

        subplot(2,1,2)
            hold on;
            for bin = 1:6
                errorbar(bin,mean(norm_all_SO_nested_spindle_cc(:,bin)-norm_all_not_nested_spindle_cc(:,bin)),std(norm_all_SO_nested_spindle_cc(:,bin)-norm_all_not_nested_spindle_cc(:,bin))/sqrt(size(norm_all_not_nested_spindle_cc,1)),'color','k')
                [h p] = ttest(norm_all_SO_nested_spindle_cc(:,bin)-norm_all_not_nested_spindle_cc(:,bin))
                if p<0.05
                    scatter(bin,2.5e-4,30,[0 0 0],'filled')
                end
            end
            plot([0 7],[0 0],'k')
            xlim([0 7])
            ylim([-1e-4 3.5e-4])

    tmp_SO_nested_m1_fr = [];
    tmp_SO_not_nested_m1_fr = [];

        all_SO_nested_spindle = [pre_SO_nested_spindle_raster_m1; post_SO_nested_spindle_raster_m1];
        all_not_nested_spindle = [pre_SO_not_nested_spindle_raster_m1; post_SO_not_nested_spindle_raster_m1];
        for n = 1:size(all_SO_nested_spindle,1)
            tmp_SO_nested_m1_fr = [tmp_SO_nested_m1_fr; sum(reshape(all_SO_nested_spindle(n,:),100,300))];
            tmp_SO_not_nested_m1_fr = [tmp_SO_not_nested_m1_fr; sum(reshape(all_not_nested_spindle(n,:),100,300))];    
        end
        tmp_SO_nested_m1_fr = tmp_SO_nested_m1_fr(:,59:240);
        tmp_SO_not_nested_m1_fr = tmp_SO_not_nested_m1_fr(:,59:240);

        tmp_SO_nested_m1_fr = tmp_SO_nested_m1_fr-repmat(nanmean(tmp_SO_nested_m1_fr(:,1:30)'),182,1)';
        tmp_SO_not_nested_m1_fr = tmp_SO_not_nested_m1_fr-repmat(nanmean(tmp_SO_not_nested_m1_fr(:,1:30)'),182,1)';

    tmp_SO_nested_dls_fr = [];
    tmp_SO_not_nested_dls_fr = [];

        all_SO_nested_spindle = [pre_SO_nested_spindle_raster_dls; post_SO_nested_spindle_raster_dls];
        all_not_nested_spindle = [pre_SO_not_nested_spindle_raster_dls; post_SO_not_nested_spindle_raster_dls];
        for n = 1:size(all_SO_nested_spindle,1)
            tmp_SO_nested_dls_fr = [tmp_SO_nested_dls_fr; sum(reshape(all_SO_nested_spindle(n,:),100,300))];
            tmp_SO_not_nested_dls_fr = [tmp_SO_not_nested_dls_fr; sum(reshape(all_not_nested_spindle(n,:),100,300))];    
        end

        tmp_SO_nested_dls_fr = tmp_SO_nested_dls_fr(:,59:240);
        tmp_SO_not_nested_dls_fr = tmp_SO_not_nested_dls_fr(:,59:240);

        tmp_SO_nested_dls_fr = tmp_SO_nested_dls_fr-repmat(nanmean(tmp_SO_nested_dls_fr(:,1:30)'),182,1)';
        tmp_SO_not_nested_dls_fr = tmp_SO_not_nested_dls_fr-repmat(nanmean(tmp_SO_not_nested_dls_fr(:,1:30)'),182,1)';

    figure;

        subplot(1,2,1); hold on;
        plot(nanmean(tmp_SO_nested_m1_fr)+nanstd(tmp_SO_nested_m1_fr)/sqrt(size(tmp_SO_nested_m1_fr,1)))
        plot(nanmean(tmp_SO_nested_m1_fr)-nanstd(tmp_SO_nested_m1_fr)/sqrt(size(tmp_SO_nested_m1_fr,1)))
        ylim([-1 1])
        xlim([31 182])

        subplot(1,2,2); hold on;
        plot(nanmean(tmp_SO_not_nested_m1_fr)+nanstd(tmp_SO_not_nested_m1_fr)/sqrt(size(tmp_SO_not_nested_m1_fr,1)))
        plot(nanmean(tmp_SO_not_nested_m1_fr)-nanstd(tmp_SO_not_nested_m1_fr)/sqrt(size(tmp_SO_not_nested_m1_fr,1)))
        ylim([-1 1])
        xlim([31 182])

    figure;

        subplot(1,2,1); hold on;
        plot(nanmean(tmp_SO_nested_dls_fr)+nanstd(tmp_SO_nested_dls_fr)/sqrt(size(tmp_SO_nested_dls_fr,1)))
        plot(nanmean(tmp_SO_nested_dls_fr)-nanstd(tmp_SO_nested_dls_fr)/sqrt(size(tmp_SO_nested_dls_fr,1)))
        ylim([-.5 .5])
        xlim([31 182])

        subplot(1,2,2); hold on;
        plot(nanmean(tmp_SO_not_nested_dls_fr)+nanstd(tmp_SO_not_nested_dls_fr)/sqrt(size(tmp_SO_not_nested_dls_fr,1)))
        plot(nanmean(tmp_SO_not_nested_dls_fr)-nanstd(tmp_SO_not_nested_dls_fr)/sqrt(size(tmp_SO_not_nested_dls_fr,1)))
        ylim([-.5 .5])
        xlim([31 182])

%% sFIG 1 | SPINDLE PROBABILITY BEFORE/AFTER NESTED OR NON-NESTED SPINDLES

    animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
    animal_days = {[1:8],[1:6],[2:10],[4:17],[4:13],[1:5 7:9]};

    spindles_around_nested = [];
    spindles_around_iso = [];
    SO_around_nested = [];
    SO_around_iso = [];

    for animal = 1:length(animal_id)

        CC_path = [data_path '\' animal_id{animal} '\cross_correlations'];

        for day = animal_days{animal}

            load([data_path '\' animal_id{animal} '\sleep_activity\day_'  num2str(day) '_block_1_sleep_activity.mat']);

            pre_SO_nesting = [];
            for ind = 1:length(sleep_rhythms.spindles{1,1}.pks)
                tmp_diff_SO = sleep_rhythms.spindles{1,1}.pks(ind)-sleep_rhythms.so_delta.SO_zc;
                if ~isempty(tmp_diff_SO(tmp_diff_SO>0))
                    pre_SO_nesting = [pre_SO_nesting min(tmp_diff_SO(tmp_diff_SO>0))];
                end
            end
            
            for ind = find(pre_SO_nesting<1)
                if length(sleep_rhythms.spindles{1,1}.pks(sleep_rhythms.spindles{1,1}.pks>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.spindles{1,1}.pks<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind))>1
                    spindles_around_nested = [spindles_around_nested; histc(sleep_rhythms.spindles{1,1}.pks(sleep_rhythms.spindles{1,1}.pks>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.spindles{1,1}.pks<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind),[-91:1:91])'];
                elseif length(sleep_rhythms.spindles{1,1}.pks(sleep_rhythms.spindles{1,1}.pks>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.spindles{1,1}.pks<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind))==1
                    spindles_around_nested = [spindles_around_nested; histc(sleep_rhythms.spindles{1,1}.pks(sleep_rhythms.spindles{1,1}.pks>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.spindles{1,1}.pks<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind),[-91:1:91])];
                end
            end
            for ind = find(pre_SO_nesting>5)
                if length(sleep_rhythms.spindles{1,1}.pks(sleep_rhythms.spindles{1,1}.pks>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.spindles{1,1}.pks<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind))>1
                    spindles_around_iso = [spindles_around_iso; histc(sleep_rhythms.spindles{1,1}.pks(sleep_rhythms.spindles{1,1}.pks>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.spindles{1,1}.pks<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind),[-91:1:91])'];
                elseif length(sleep_rhythms.spindles{1,1}.pks(sleep_rhythms.spindles{1,1}.pks>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.spindles{1,1}.pks<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind))==1
                    spindles_around_iso = [spindles_around_iso; histc(sleep_rhythms.spindles{1,1}.pks(sleep_rhythms.spindles{1,1}.pks>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.spindles{1,1}.pks<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind),[-91:1:91])];
                end
            end

            for ind = find(pre_SO_nesting<1)
                if length(sleep_rhythms.so_delta.SO_zc(sleep_rhythms.so_delta.SO_zc>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.so_delta.SO_zc<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind))>1
                    SO_around_nested = [SO_around_nested; histc(sleep_rhythms.so_delta.SO_zc(sleep_rhythms.so_delta.SO_zc>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.so_delta.SO_zc<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind),[-91:1:91])'];
                elseif length(sleep_rhythms.so_delta.SO_zc(sleep_rhythms.so_delta.SO_zc>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.so_delta.SO_zc<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind))==1
                    SO_around_nested = [SO_around_nested; histc(sleep_rhythms.so_delta.SO_zc(sleep_rhythms.so_delta.SO_zc>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.so_delta.SO_zc<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind),[-91:1:91])];
                end
            end
            for ind = find(pre_SO_nesting>5)
                if length(sleep_rhythms.so_delta.SO_zc(sleep_rhythms.so_delta.SO_zc>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.so_delta.SO_zc<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind))>1
                    SO_around_iso = [SO_around_iso; histc(sleep_rhythms.so_delta.SO_zc(sleep_rhythms.so_delta.SO_zc>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.so_delta.SO_zc<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind),[-91:1:91])'];
                elseif length(sleep_rhythms.so_delta.SO_zc(sleep_rhythms.so_delta.SO_zc>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.so_delta.SO_zc<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind))==1
                    SO_around_iso = [SO_around_iso; histc(sleep_rhythms.so_delta.SO_zc(sleep_rhythms.so_delta.SO_zc>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.so_delta.SO_zc<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind),[-91:1:91])];
                end
            end
            
            if animal == 1 && day == 6
                load([data_path '\' animal_id{animal} '\sleep_activity\day_'  num2str(day) '_block_3_sleep_activity.mat']);
            else
                load([data_path '\' animal_id{animal} '\sleep_activity\day_'  num2str(day) '_block_2_sleep_activity.mat']);
            end
            
            post_SO_nesting = [];
            for ind = 1:length(sleep_rhythms.spindles{1,1}.pks)
                tmp_diff_SO = sleep_rhythms.spindles{1,1}.pks(ind)-sleep_rhythms.so_delta.SO_zc;
                if ~isempty(tmp_diff_SO(tmp_diff_SO>0))
                    post_SO_nesting = [post_SO_nesting min(tmp_diff_SO(tmp_diff_SO>0))];
                end
            end

            for ind = find(post_SO_nesting<1)
                if length(sleep_rhythms.spindles{1,1}.pks(sleep_rhythms.spindles{1,1}.pks>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.spindles{1,1}.pks<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind))>1
                    spindles_around_nested = [spindles_around_nested; histc(sleep_rhythms.spindles{1,1}.pks(sleep_rhythms.spindles{1,1}.pks>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.spindles{1,1}.pks<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind),[-91:1:91])'];
                elseif length(sleep_rhythms.spindles{1,1}.pks(sleep_rhythms.spindles{1,1}.pks>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.spindles{1,1}.pks<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind))==1
                    spindles_around_nested = [spindles_around_nested; histc(sleep_rhythms.spindles{1,1}.pks(sleep_rhythms.spindles{1,1}.pks>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.spindles{1,1}.pks<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind),[-91:1:91])];
                end
            end
            for ind = find(post_SO_nesting>5)
                if length(sleep_rhythms.spindles{1,1}.pks(sleep_rhythms.spindles{1,1}.pks>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.spindles{1,1}.pks<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind))>1
                    spindles_around_iso = [spindles_around_iso; histc(sleep_rhythms.spindles{1,1}.pks(sleep_rhythms.spindles{1,1}.pks>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.spindles{1,1}.pks<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind),[-91:1:91])'];
                elseif length(sleep_rhythms.spindles{1,1}.pks(sleep_rhythms.spindles{1,1}.pks>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.spindles{1,1}.pks<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind))==1
                    spindles_around_iso = [spindles_around_iso; histc(sleep_rhythms.spindles{1,1}.pks(sleep_rhythms.spindles{1,1}.pks>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.spindles{1,1}.pks<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind),[-91:1:91])];
                end
            end

            for ind = find(post_SO_nesting<1)
                if length(sleep_rhythms.so_delta.SO_zc(sleep_rhythms.so_delta.SO_zc>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.so_delta.SO_zc<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind))>1
                    SO_around_nested = [SO_around_nested; histc(sleep_rhythms.so_delta.SO_zc(sleep_rhythms.so_delta.SO_zc>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.so_delta.SO_zc<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind),[-91:1:91])'];
                elseif length(sleep_rhythms.so_delta.SO_zc(sleep_rhythms.so_delta.SO_zc>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.so_delta.SO_zc<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind))==1
                    SO_around_nested = [SO_around_nested; histc(sleep_rhythms.so_delta.SO_zc(sleep_rhythms.so_delta.SO_zc>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.so_delta.SO_zc<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind),[-91:1:91])];
                end
            end
            for ind = find(post_SO_nesting>5)
                if length(sleep_rhythms.so_delta.SO_zc(sleep_rhythms.so_delta.SO_zc>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.so_delta.SO_zc<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind))>1
                    SO_around_iso = [SO_around_iso; histc(sleep_rhythms.so_delta.SO_zc(sleep_rhythms.so_delta.SO_zc>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.so_delta.SO_zc<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind),[-91:1:91])'];
                elseif length(sleep_rhythms.so_delta.SO_zc(sleep_rhythms.so_delta.SO_zc>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.so_delta.SO_zc<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind))==1
                    SO_around_iso = [SO_around_iso; histc(sleep_rhythms.so_delta.SO_zc(sleep_rhythms.so_delta.SO_zc>sleep_rhythms.spindles{1,1}.pks(ind)-91 & sleep_rhythms.so_delta.SO_zc<sleep_rhythms.spindles{1,1}.pks(ind)+91)-sleep_rhythms.spindles{1,1}.pks(ind),[-91:1:91])];
                end
            end
            
        end
    end

    spindles_around_nested = spindles_around_nested-repmat(nanmean(spindles_around_nested(:,1:30)'),183,1)';
    spindles_around_iso = spindles_around_iso-repmat(nanmean(spindles_around_iso(:,1:30)'),183,1)';
    
    figure;
        
        subplot(1,2,1); hold on;
        plot(nanmean(spindles_around_nested)+nanstd(spindles_around_nested)/sqrt(size(spindles_around_nested,1)))
        plot(nanmean(spindles_around_nested)-nanstd(spindles_around_nested)/sqrt(size(spindles_around_nested,1)))
        xlim([31 182])
        
        subplot(1,2,2); hold on;        
        plot(nanmean(spindles_around_iso)+nanstd(spindles_around_iso)/sqrt(size(spindles_around_iso,1)))
        plot(nanmean(spindles_around_iso)-nanstd(spindles_around_iso)/sqrt(size(spindles_around_iso,1)))        
        xlim([31 182])
        