%% SCRIPT TO GENERATE FIGURE 6 & FIGURE 6 FIGURE SUPPLEMENTS
% 
% to generate figures, download data from https://doi.org/10.7272/Q6KK9927

clc; clear all; close all;
data_path = 'D:\Lemke_etal_2021\Lemke_eLife_2021_data';

%% FIG 6 A, B, C, & D | PRE POST CC CHANGE 

    animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
    animal_days = {[1:8],[1:6],[2:10],[4:17],[4:13],[1:5 7:9]};

    all_animal_pre_cc_1 = cell(1,6);
    all_animal_pre_cc_2 = cell(1,6);
    all_animal_post_cc_1 = cell(1,6);
    all_animal_post_cc_2 = cell(1,6);

    all_animal_pre_dls_spindle = cell(1,6);
    all_animal_post_dls_spindle = cell(1,6);
    all_animal_pre_m1_spindle = cell(1,6);
    all_animal_post_m1_spindle = cell(1,6);

    all_animal_dls_reach_mod = cell(1,6);
    all_animal_m1_reach_mod = cell(1,6);

    all_animal_coh = cell(1,6);

    theta_band = [48 113; 48 113; 48 113; 48 113; 40 95; 40 95];

    for animal = 1:length(animal_id)

        %%% LFP COHERENCE

            load([data_path '\' animal_id{animal} '\sleep_activity\LFP_coherence.mat'])

            increase_channels = [];
            for pairs=1:length(pre_LFP_coh)
                if mean(mean([pre_LFP_coh{pairs}(end,theta_band(animal,1):theta_band(animal,2)); post_LFP_coh{pairs}(end,theta_band(animal,1):theta_band(animal,2))])) > (mean(mean([pre_LFP_coh{pairs}(1,theta_band(animal,1):theta_band(animal,2)); post_LFP_coh{pairs}(1,theta_band(animal,1):theta_band(animal,2))])))+0.025
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
                for day=1:length(animal_days{animal})
                    tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,theta_band(animal,1):theta_band(animal,2)))];
                    tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,theta_band(animal,1):theta_band(animal,2)))];
                end
                day_pre = [day_pre; tmp_pre];
                day_post = [day_post; tmp_post];
            end

            all_animal_coh{animal} = [mean(day_pre); mean(day_post)];

        %%% CROSS CORRELATION CHANGE AND REACH MOD

            if animal == 6
                pre_cc_1 = cell(1,max(animal_days{animal})-min(animal_days{animal}));
                pre_cc_2 = cell(1,max(animal_days{animal})-min(animal_days{animal}));
                post_cc_1 = cell(1,max(animal_days{animal})-min(animal_days{animal}));
                post_cc_2 = cell(1,max(animal_days{animal})-min(animal_days{animal}));

                pre_dls_spindle_raster = cell(1,max(animal_days{animal})-min(animal_days{animal}));
                post_dls_spindle_raster = cell(1,max(animal_days{animal})-min(animal_days{animal}));
                pre_m1_spindle_raster = cell(1,max(animal_days{animal})-min(animal_days{animal}));
                post_m1_spindle_raster = cell(1,max(animal_days{animal})-min(animal_days{animal}));

                dls_reach_mod = cell(1,max(animal_days{animal})-min(animal_days{animal}));
                m1_reach_mod = cell(1,max(animal_days{animal})-min(animal_days{animal}));
            else
                pre_cc_1 = cell(1,1+max(animal_days{animal})-min(animal_days{animal}));
                pre_cc_2 = cell(1,1+max(animal_days{animal})-min(animal_days{animal}));
                post_cc_1 = cell(1,1+max(animal_days{animal})-min(animal_days{animal}));
                post_cc_2 = cell(1,1+max(animal_days{animal})-min(animal_days{animal}));

                pre_dls_spindle_raster = cell(1,1+max(animal_days{animal})-min(animal_days{animal}));
                post_dls_spindle_raster = cell(1,1+max(animal_days{animal})-min(animal_days{animal}));
                pre_m1_spindle_raster = cell(1,1+max(animal_days{animal})-min(animal_days{animal}));
                post_m1_spindle_raster = cell(1,1+max(animal_days{animal})-min(animal_days{animal}));

                dls_reach_mod = cell(1,1+max(animal_days{animal})-min(animal_days{animal}));
                m1_reach_mod = cell(1,1+max(animal_days{animal})-min(animal_days{animal}));
            end

            day_count = 1;
            for day = animal_days{animal}

                load([data_path '\' animal_id{animal} '\cross_correlations\day_' num2str(day) '\CC_by_sleep_rhythm.mat']);
                load([data_path '\' animal_id{animal} '\cross_correlations\day_' num2str(day) '\behavioral_state_cc.mat']);
                sig = load([data_path '\' animal_id{animal} '\cross_correlations\day_' num2str(day) '\jitter_test.mat']);

                %%% REACH MOD 

                    if ~isfile([data_path '\' animal_id{animal} '\reach_activity\day_' num2str(day) '_reach_mod.mat'])

                        tmp_m1_sig = zeros(1,length(all_m1_spindle));
                        tmp_dls_sig = zeros(1,length(all_dls_spindle));
                        
                    else
                        
                        load([data_path '\' animal_id{animal} '\reach_activity\day_' num2str(day) '_reach_mod.mat']);

                        shuffled_ind = [];
                        for shuffle = 1:1000
                            shuffled_ind = [shuffled_ind; randperm(400)];
                        end

                        tmp_m1_sig = [];
                        for unit = 1:length(all_m1_spindle)
                            tmp_raster = (mean(m1_reach(unit).raster)-mean(mean(m1_reach(unit).raster(:,1:100))))/std(mean(m1_reach(unit).raster(:,1:100)));
                            shuffled_mod = sum(abs(tmp_raster(shuffled_ind(:,161:240)))');
                            if sum(abs(tmp_raster(161:240)))>prctile(shuffled_mod,95)
                                tmp_m1_sig = [tmp_m1_sig 1];
                            else
                                tmp_m1_sig = [tmp_m1_sig 0];
                            end
                        end

                        tmp_dls_sig = [];        
                        for unit = 1:length(all_dls_spindle)
                            tmp_raster = (mean(dls_reach(unit).raster)-mean(mean(dls_reach(unit).raster(:,1:100))))/std(mean(dls_reach(unit).raster(:,1:100)));
                            shuffled_mod = sum(abs(tmp_raster(shuffled_ind(:,161:240)))');
                            if sum(abs(tmp_raster(161:240)))>prctile(shuffled_mod,95)
                                tmp_dls_sig = [tmp_dls_sig 1];
                            else
                                tmp_dls_sig = [tmp_dls_sig 0];
                            end
                        end
                        
                    end

                %%% SPINDLE MOD

                    tmp_pre_m1_spindle_raster = [];
                    tmp_post_m1_spindle_raster = [];
                    for unit = 1:length(all_m1_spindle)
                        tmp_pre_m1_spindle_raster = [tmp_pre_m1_spindle_raster; mean(all_m1_spindle(unit).pre_m1_spindle_raster)];
                        tmp_post_m1_spindle_raster = [tmp_post_m1_spindle_raster; mean(all_m1_spindle(unit).post_m1_spindle_raster)];
                    end

                    tmp_pre_dls_spindle_raster = [];
                    tmp_post_dls_spindle_raster = [];        
                    for unit = 1:length(all_dls_spindle)
                        tmp_pre_dls_spindle_raster = [tmp_pre_dls_spindle_raster; mean(all_dls_spindle(unit).pre_dls_spindle_raster)];
                        tmp_post_dls_spindle_raster = [tmp_post_dls_spindle_raster; mean(all_dls_spindle(unit).post_dls_spindle_raster)];
                    end

                %%% CROSS CORRELATION

                    day_pre_dls_spindle_raster = [];
                    day_post_dls_spindle_raster = [];
                    day_pre_m1_spindle_raster = [];
                    day_post_m1_spindle_raster = [];            

                    day_pre_cc_1 = [];
                    day_pre_cc_2 = [];
                    day_post_cc_1 = [];
                    day_post_cc_2 = [];

                    day_reach_mod_m1 = [];
                    day_reach_mod_dls = [];

                    pair_count = 1;
                    for m1_unit = 1:length(all_m1_spindle)
                        for dls_unit = 1:length(all_dls_spindle)
                            if sig.all_pairs(pair_count).p_val>0.99

                                day_pre_dls_spindle_raster = [day_pre_dls_spindle_raster; tmp_pre_dls_spindle_raster(dls_unit,:)];
                                day_post_dls_spindle_raster = [day_post_dls_spindle_raster; tmp_post_dls_spindle_raster(dls_unit,:)];
                                day_pre_m1_spindle_raster = [day_pre_m1_spindle_raster; tmp_pre_m1_spindle_raster(m1_unit,:)];
                                day_post_m1_spindle_raster = [day_post_m1_spindle_raster; tmp_post_m1_spindle_raster(m1_unit,:)];

                                day_pre_cc_1 = [day_pre_cc_1; mean(repmat(all_pairs(pair_count).pre_CC_nrem_1,100,1)-all_pairs(pair_count).shuffle_pre_CC_nrem_1)];
                                day_pre_cc_2 = [day_pre_cc_2; mean(repmat(all_pairs(pair_count).pre_CC_nrem_2,100,1)-all_pairs(pair_count).shuffle_pre_CC_nrem_2)];
                                day_post_cc_1 = [day_post_cc_1; mean(repmat(all_pairs(pair_count).post_CC_nrem_1,100,1)-all_pairs(pair_count).shuffle_post_CC_nrem_1)];
                                day_post_cc_2 = [day_post_cc_2; mean(repmat(all_pairs(pair_count).post_CC_nrem_2,100,1)-all_pairs(pair_count).shuffle_post_CC_nrem_2)];

                                day_reach_mod_m1 = [day_reach_mod_m1 tmp_m1_sig(m1_unit)];
                                day_reach_mod_dls = [day_reach_mod_dls tmp_dls_sig(dls_unit)];                        

                            end
                            pair_count = pair_count + 1;
                        end
                    end

                    pre_cc_1{day_count} = day_pre_cc_1;
                    pre_cc_2{day_count} = day_pre_cc_2;
                    post_cc_1{day_count} = day_post_cc_1;
                    post_cc_2{day_count} = day_post_cc_2;

                    pre_dls_spindle_raster{day_count} = day_pre_dls_spindle_raster;
                    post_dls_spindle_raster{day_count} = day_post_dls_spindle_raster;
                    pre_m1_spindle_raster{day_count} = day_pre_m1_spindle_raster;
                    post_m1_spindle_raster{day_count} = day_post_m1_spindle_raster;

                    m1_reach_mod{day_count} = day_reach_mod_m1;
                    dls_reach_mod{day_count} = day_reach_mod_dls;

                    day_count = day_count + 1;

            end

            all_animal_pre_cc_1{animal} = pre_cc_1;
            all_animal_pre_cc_2{animal} = pre_cc_2;
            all_animal_post_cc_1{animal} = post_cc_1;
            all_animal_post_cc_2{animal} = post_cc_2;

            all_animal_pre_dls_spindle{animal} = pre_dls_spindle_raster;
            all_animal_post_dls_spindle{animal} = post_dls_spindle_raster;
            all_animal_pre_m1_spindle{animal} = pre_m1_spindle_raster;
            all_animal_post_m1_spindle{animal} = post_m1_spindle_raster;

            all_animal_m1_reach_mod{animal} = m1_reach_mod;
            all_animal_dls_reach_mod{animal} = dls_reach_mod;

    end
 
    %%% PLOT

        pre_cc_1 = [];
        pre_cc_2 = [];
        post_cc_1 = [];
        post_cc_2 = [];

        pre_cc_1_m1_dls_spindlemod = [];
        pre_cc_2_m1_dls_spindlemod = [];
        post_cc_1_m1_dls_spindlemod = [];
        post_cc_2_m1_dls_spindlemod = [];

        pre_cc_1_m1_spindlemod = [];
        pre_cc_2_m1_spindlemod = [];
        post_cc_1_m1_spindlemod = [];
        post_cc_2_m1_spindlemod = [];

        pre_cc_1_dls_spindlemod = [];
        pre_cc_2_dls_spindlemod = [];
        post_cc_1_dls_spindlemod = [];
        post_cc_2_dls_spindlemod = [];

        pre_cc_1_no_spindlemod = [];
        pre_cc_2_no_spindlemod = [];
        post_cc_1_no_spindlemod = [];
        post_cc_2_no_spindlemod = [];

        for animal = 1:length(animal_id)
            for day = 1:length(animal_days{animal})

                shuffled_ind = [];
                for shuffle = 1:1000
                    shuffled_ind = [shuffled_ind; randperm(200)];
                end

                pre_m1_spindle_sig = [];
                post_m1_spindle_sig = [];
                pre_dls_spindle_sig = [];
                post_dls_spindle_sig = [];

                for pair = 1:size(all_animal_pre_dls_spindle{animal}{day},1)    

                    tmp_spindle = (all_animal_pre_m1_spindle{animal}{day}(pair,:)-mean(all_animal_pre_m1_spindle{animal}{day}(pair,1:60)))/std(all_animal_pre_m1_spindle{animal}{day}(pair,1:60));
                    shuffled_mod = sum(abs(tmp_spindle(shuffled_ind(:,61:140)))');
                    if sum(abs(tmp_spindle(61:140)))>prctile(shuffled_mod,95)
                        pre_m1_spindle_sig = [pre_m1_spindle_sig 1];
                    else
                        pre_m1_spindle_sig = [pre_m1_spindle_sig 0];
                    end

                    tmp_spindle = (all_animal_post_m1_spindle{animal}{day}(pair,:)-mean(all_animal_post_m1_spindle{animal}{day}(pair,1:60)))/std(all_animal_post_m1_spindle{animal}{day}(pair,1:60));
                    shuffled_mod = sum(abs(tmp_spindle(shuffled_ind(:,61:140)))');
                    if sum(abs(tmp_spindle(61:140)))>prctile(shuffled_mod,95)
                        post_m1_spindle_sig = [post_m1_spindle_sig 1];
                    else
                        post_m1_spindle_sig = [post_m1_spindle_sig 0];
                    end   

                    tmp_spindle = (all_animal_pre_dls_spindle{animal}{day}(pair,:)-mean(all_animal_pre_dls_spindle{animal}{day}(pair,1:60)))/std(all_animal_pre_dls_spindle{animal}{day}(pair,1:60));
                    shuffled_mod = sum(abs(tmp_spindle(shuffled_ind(:,61:140)))');
                    if sum(abs(tmp_spindle(61:140)))>prctile(shuffled_mod,95)
                        pre_dls_spindle_sig = [pre_dls_spindle_sig 1];
                    else
                        pre_dls_spindle_sig = [pre_dls_spindle_sig 0];
                    end

                    tmp_spindle = (all_animal_post_dls_spindle{animal}{day}(pair,:)-mean(all_animal_post_dls_spindle{animal}{day}(pair,1:60)))/std(all_animal_post_dls_spindle{animal}{day}(pair,1:60));
                    shuffled_mod = sum(abs(tmp_spindle(shuffled_ind(:,61:140)))');
                    if sum(abs(tmp_spindle(61:140)))>prctile(shuffled_mod,95)
                        post_dls_spindle_sig = [post_dls_spindle_sig 1];
                    else
                        post_dls_spindle_sig = [post_dls_spindle_sig 0];
                    end      

                end

                pre_cc_1 = [pre_cc_1; all_animal_pre_cc_1{animal}{day}];
                pre_cc_2 = [pre_cc_2; all_animal_pre_cc_2{animal}{day}];
                post_cc_1 = [post_cc_1; all_animal_post_cc_1{animal}{day}];
                post_cc_2 = [post_cc_2; all_animal_post_cc_2{animal}{day}];

                pre_cc_1_m1_dls_spindlemod = [pre_cc_1_m1_dls_spindlemod; all_animal_pre_cc_1{animal}{day}(pre_m1_spindle_sig==1 & pre_dls_spindle_sig==1,:)];
                pre_cc_2_m1_dls_spindlemod = [pre_cc_2_m1_dls_spindlemod; all_animal_pre_cc_2{animal}{day}(pre_m1_spindle_sig==1 & pre_dls_spindle_sig==1,:)];
                post_cc_1_m1_dls_spindlemod = [post_cc_1_m1_dls_spindlemod; all_animal_post_cc_1{animal}{day}(post_m1_spindle_sig==1 & post_dls_spindle_sig==1,:)];
                post_cc_2_m1_dls_spindlemod = [post_cc_2_m1_dls_spindlemod; all_animal_post_cc_2{animal}{day}(post_m1_spindle_sig==1 & post_dls_spindle_sig==1,:)];

                pre_cc_1_m1_spindlemod = [pre_cc_1_m1_spindlemod; all_animal_pre_cc_1{animal}{day}(pre_m1_spindle_sig==1 & pre_dls_spindle_sig==0,:)];
                pre_cc_2_m1_spindlemod = [pre_cc_2_m1_spindlemod; all_animal_pre_cc_2{animal}{day}(pre_m1_spindle_sig==1 & pre_dls_spindle_sig==0,:)];
                post_cc_1_m1_spindlemod = [post_cc_1_m1_spindlemod; all_animal_post_cc_1{animal}{day}(post_m1_spindle_sig==1 & post_dls_spindle_sig==0,:)];
                post_cc_2_m1_spindlemod = [post_cc_2_m1_spindlemod; all_animal_post_cc_2{animal}{day}(post_m1_spindle_sig==1 & post_dls_spindle_sig==0,:)];

                pre_cc_1_dls_spindlemod = [pre_cc_1_dls_spindlemod; all_animal_pre_cc_1{animal}{day}(pre_m1_spindle_sig==0 & pre_dls_spindle_sig==1,:)];
                pre_cc_2_dls_spindlemod = [pre_cc_2_dls_spindlemod; all_animal_pre_cc_2{animal}{day}(pre_m1_spindle_sig==0 & pre_dls_spindle_sig==1,:)];
                post_cc_1_dls_spindlemod = [post_cc_1_dls_spindlemod; all_animal_post_cc_1{animal}{day}(post_m1_spindle_sig==0 & post_dls_spindle_sig==1,:)];
                post_cc_2_dls_spindlemod = [post_cc_2_dls_spindlemod; all_animal_post_cc_2{animal}{day}(post_m1_spindle_sig==0 & post_dls_spindle_sig==1,:)];

                pre_cc_1_no_spindlemod = [pre_cc_1_no_spindlemod; all_animal_pre_cc_1{animal}{day}(pre_m1_spindle_sig==0 | pre_dls_spindle_sig==0,:)];
                pre_cc_2_no_spindlemod = [pre_cc_2_no_spindlemod; all_animal_pre_cc_2{animal}{day}(pre_m1_spindle_sig==0 | pre_dls_spindle_sig==0,:)];
                post_cc_1_no_spindlemod = [post_cc_1_no_spindlemod; all_animal_post_cc_1{animal}{day}(post_m1_spindle_sig==0 | post_dls_spindle_sig==0,:)];
                post_cc_2_no_spindlemod = [post_cc_2_no_spindlemod; all_animal_post_cc_2{animal}{day}(post_m1_spindle_sig==0 | post_dls_spindle_sig==0,:)];

            end
        end

        fig = figure;
            subplot(2,2,1);
                hold on;
                plot(nanmean(pre_cc_1_m1_dls_spindlemod)-(nanstd(pre_cc_1_m1_dls_spindlemod)/sqrt(size(pre_cc_1_m1_dls_spindlemod,1))),'k')
                plot(nanmean(pre_cc_1_m1_dls_spindlemod)+(nanstd(pre_cc_1_m1_dls_spindlemod)/sqrt(size(pre_cc_1_m1_dls_spindlemod,1))),'k')
                plot([101 101],[-8e-4 15e-4],'k');
                title(['PRE CC 1 SPINDLE MOD - ' num2str(size(pre_cc_1_m1_dls_spindlemod,1)) ' pairs'])
                xlim([51 150]);  
                ylim([-12e-4 20e-4])
            subplot(2,2,2);
                hold on;
                plot(nanmean(pre_cc_2_m1_dls_spindlemod)-(nanstd(pre_cc_2_m1_dls_spindlemod)/sqrt(size(pre_cc_2_m1_dls_spindlemod,1))),'k')
                plot(nanmean(pre_cc_2_m1_dls_spindlemod)+(nanstd(pre_cc_2_m1_dls_spindlemod)/sqrt(size(pre_cc_2_m1_dls_spindlemod,1))),'k')
                plot([101 101],[-8e-4 15e-4],'k');
                title(['PRE CC 2 SPINDLE MOD - ' num2str(size(pre_cc_2_m1_dls_spindlemod,1)) ' pairs'])
                xlim([51 150]);  
                ylim([-12e-4 20e-4])
            subplot(2,2,3);
                hold on;
                plot(nanmean(post_cc_1_m1_dls_spindlemod)-(nanstd(post_cc_1_m1_dls_spindlemod)/sqrt(size(post_cc_1_m1_dls_spindlemod,1))),'k')
                plot(nanmean(post_cc_1_m1_dls_spindlemod)+(nanstd(post_cc_1_m1_dls_spindlemod)/sqrt(size(post_cc_1_m1_dls_spindlemod,1))),'k')
                plot([101 101],[-8e-4 15e-4],'k');
                title(['POST CC 1 SPINDLE MOD - ' num2str(size(post_cc_1_m1_dls_spindlemod,1)) ' pairs'])
                xlim([51 150]);  
                ylim([-12e-4 20e-4])
            subplot(2,2,4);
                hold on;
                plot(nanmean(post_cc_2_m1_dls_spindlemod)-(nanstd(post_cc_2_m1_dls_spindlemod)/sqrt(size(post_cc_2_m1_dls_spindlemod,1))),'k')
                plot(nanmean(post_cc_2_m1_dls_spindlemod)+(nanstd(post_cc_2_m1_dls_spindlemod)/sqrt(size(post_cc_2_m1_dls_spindlemod,1))),'k')
                plot([101 101],[-8e-4 15e-4],'k');
                title(['POST CC 2 SPINDLE MOD - ' num2str(size(post_cc_2_m1_dls_spindlemod,1)) ' pairs'])
                xlim([51 150]);  
                ylim([-12e-4 20e-4])

        fig = figure;
            subplot(2,2,1);
                hold on;
                plot(nanmean(pre_cc_1_no_spindlemod)-(nanstd(pre_cc_1_no_spindlemod)/sqrt(size(pre_cc_1_no_spindlemod,1))),'k')
                plot(nanmean(pre_cc_1_no_spindlemod)+(nanstd(pre_cc_1_no_spindlemod)/sqrt(size(pre_cc_1_no_spindlemod,1))),'k')
                plot([101 101],[-8e-4 15e-4],'k');
                title(['PRE CC 1 NO SPINDLE MOD - ' num2str(size(pre_cc_1_no_spindlemod,1)) ' pairs'])
                xlim([51 150]);  
                ylim([-8e-4 15e-4])
            subplot(2,2,2);
                hold on;
                plot(nanmean(pre_cc_2_no_spindlemod)-(nanstd(pre_cc_2_no_spindlemod)/sqrt(size(pre_cc_2_no_spindlemod,1))),'k')
                plot(nanmean(pre_cc_2_no_spindlemod)+(nanstd(pre_cc_2_no_spindlemod)/sqrt(size(pre_cc_2_no_spindlemod,1))),'k')
                plot([101 101],[-8e-4 15e-4],'k');
                title(['PRE CC 2 NO SPINDLE MOD - ' num2str(size(pre_cc_2_no_spindlemod,1)) ' pairs'])
                xlim([51 150]);  
                ylim([-8e-4 15e-4])
            subplot(2,2,3);
                hold on;
                plot(nanmean(post_cc_1_no_spindlemod)-(nanstd(post_cc_1_no_spindlemod)/sqrt(size(post_cc_1_no_spindlemod,1))),'k')
                plot(nanmean(post_cc_1_no_spindlemod)+(nanstd(post_cc_1_no_spindlemod)/sqrt(size(post_cc_1_no_spindlemod,1))),'k')
                plot([101 101],[-8e-4 15e-4],'k');
                title(['POST CC 1 NO SPINDLE MOD - ' num2str(size(post_cc_1_no_spindlemod,1)) ' pairs'])
                xlim([51 150]);  
                ylim([-8e-4 15e-4])
            subplot(2,2,4);
                hold on;
                plot(nanmean(post_cc_2_no_spindlemod)-(nanstd(post_cc_2_no_spindlemod)/sqrt(size(post_cc_2_no_spindlemod,1))),'k')
                plot(nanmean(post_cc_2_no_spindlemod)+(nanstd(post_cc_2_no_spindlemod)/sqrt(size(post_cc_2_no_spindlemod,1))),'k')
                plot([101 101],[-8e-4 15e-4],'k');
                title(['POST CC 2 NO SPINDLE MOD - ' num2str(size(post_cc_2_no_spindlemod,1)) ' pairs'])
                xlim([51 150]);  
                ylim([-8e-4 15e-4])

        fig = figure;
            hold on;
            errorbar(1,nanmean(nanmean(pre_cc_2_m1_dls_spindlemod(:,86:100)')-nanmean(pre_cc_1_m1_dls_spindlemod(:,86:100)')), ...
                       nanstd(nanmean(pre_cc_2_m1_dls_spindlemod(:,86:100)')-nanmean(pre_cc_1_m1_dls_spindlemod(:,86:100)'))/sqrt(size(pre_cc_2_m1_dls_spindlemod,1)),'color','k')
            errorbar(2,nanmean(nanmean(post_cc_2_m1_dls_spindlemod(:,86:100)')-nanmean(post_cc_1_m1_dls_spindlemod(:,86:100)')), ...
                       nanstd(nanmean(post_cc_2_m1_dls_spindlemod(:,86:100)')-nanmean(post_cc_1_m1_dls_spindlemod(:,86:100)'))/sqrt(size(post_cc_2_m1_dls_spindlemod,1)),'color','r')
            errorbar(3,nanmean(nanmean(pre_cc_2_no_spindlemod(:,86:100)')-nanmean(pre_cc_1_no_spindlemod(:,86:100)')), ...
                       nanstd(nanmean(pre_cc_2_no_spindlemod(:,86:100)')-nanmean(pre_cc_1_no_spindlemod(:,86:100)'))/sqrt(size(pre_cc_2_no_spindlemod,1)),'color','k')
            errorbar(4,nanmean(nanmean(post_cc_2_no_spindlemod(:,86:100)')-nanmean(post_cc_1_no_spindlemod(:,86:100)')), ...
                       nanstd(nanmean(post_cc_2_no_spindlemod(:,86:100)')-nanmean(post_cc_1_no_spindlemod(:,86:100)'))/sqrt(size(post_cc_2_no_spindlemod,1)),'color','r')
            plot([0 5],[0 0],'k','LineStyle','--')
            xlim([0 5])   
            ylim([-4e-4 2e-4])   

            [h p] = ttest2(nanmean(pre_cc_2_no_spindlemod(:,86:100)')-nanmean(pre_cc_1_no_spindlemod(:,86:100)'),nanmean(pre_cc_2_m1_dls_spindlemod(:,86:100)')-nanmean(pre_cc_1_m1_dls_spindlemod(:,86:100)'))
            [h p] = ttest2(nanmean(pre_cc_2_m1_dls_spindlemod(:,86:100)')-nanmean(pre_cc_1_m1_dls_spindlemod(:,86:100)'),nanmean(post_cc_2_no_spindlemod(:,86:100)')-nanmean(post_cc_1_no_spindlemod(:,86:100)'))
            [h p] = ttest2(nanmean(pre_cc_2_no_spindlemod(:,86:100)')-nanmean(pre_cc_1_no_spindlemod(:,86:100)'),nanmean(post_cc_2_no_spindlemod(:,86:100)')-nanmean(post_cc_1_no_spindlemod(:,86:100)'))

            [h p] = ttest2(nanmean(post_cc_2_m1_dls_spindlemod(:,86:100)')-nanmean(post_cc_1_m1_dls_spindlemod(:,86:100)'),nanmean(post_cc_2_no_spindlemod(:,86:100)')-nanmean(post_cc_1_no_spindlemod(:,86:100)'))
            [h p] = ttest2(nanmean(post_cc_2_m1_dls_spindlemod(:,86:100)')-nanmean(post_cc_1_m1_dls_spindlemod(:,86:100)'),nanmean(pre_cc_2_m1_dls_spindlemod(:,86:100)')-nanmean(pre_cc_1_m1_dls_spindlemod(:,86:100)'))
            [h p] = ttest2(nanmean(post_cc_2_m1_dls_spindlemod(:,86:100)')-nanmean(post_cc_1_m1_dls_spindlemod(:,86:100)'),nanmean(pre_cc_2_no_spindlemod(:,86:100)')-nanmean(pre_cc_1_no_spindlemod(:,86:100)'))

    %%% PLOT CC VS. COH

        pre_cc_change = [];
        post_cc_change = [];

        coh_change = [];

        for animal = 1:length(animal_id)

            tmp_pre_coh = [];
            tmp_post_coh = [];

            tmp_pre_cc_change = [];
            tmp_post_cc_change = [];  

            for day = 1:length(animal_days{animal})

                tmp_pre_coh = [tmp_pre_coh all_animal_coh{animal}(1,day)];
                tmp_post_coh = [tmp_post_coh all_animal_coh{animal}(2,day)];

                if ~isempty(all_animal_pre_cc_2{animal}{day})

                    tmp_pre_cc_change = [tmp_pre_cc_change mean(mean(all_animal_pre_cc_2{animal}{day}(:,91:100)))-mean(mean(all_animal_pre_cc_1{animal}{day}(:,91:100)))];
                    tmp_post_cc_change = [tmp_post_cc_change mean(mean(all_animal_post_cc_2{animal}{day}(:,91:100)))-mean(mean(all_animal_post_cc_1{animal}{day}(:,91:100)))];

                else

                    tmp_pre_cc_change = [tmp_pre_cc_change NaN];
                    tmp_post_cc_change = [tmp_post_cc_change NaN];

                end
            end

    %         if animal == 6 
    %             pre_cc_change = [pre_cc_change tmp_pre_cc_change([1:4 6 7])];
    %             post_cc_change = [post_cc_change tmp_post_cc_change([1:4 6 7])];
    % 
    %             coh_change = [coh_change tmp_pre_coh([2:5 7 8])-tmp_post_coh(([1:4 6 7]))];
    %         else
                pre_cc_change = [pre_cc_change tmp_pre_cc_change(1:end-1)];
                post_cc_change = [post_cc_change tmp_post_cc_change(1:end-1)];

                coh_change = [coh_change tmp_pre_coh(2:end)-tmp_post_coh(1:end-1)];
    %         end

        end

        fig = figure;
            hold on;
            y = coh_change(~isnan(pre_cc_change) & ~isnan(post_cc_change));
            x = post_cc_change(~isnan(pre_cc_change) & ~isnan(post_cc_change)) - pre_cc_change(~isnan(pre_cc_change) & ~isnan(post_cc_change));     
            scatter(x,y,30,[0 0 0],'filled')
            [R,P] = corrcoef(x,y)
            [p] = polyfit(x,y,1)
            x_2 = linspace(min(x)-std(x),max(x)+std(x),100);
            y_2 = x_2*p(1)+p(2)
            plot(x_2,y_2,'color','r','LineStyle','--','LineWidth',3)
            title(['R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            xlim([min(x)-std(x) max(x)+std(x)])
       
%% sFIG 1 | LENGTH OF EACH BEHAVIORAL STATE

    animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
    animal_days = {[1:8],[1:6],[2:10],[4:17],[4:13],[1:5 7:9]};

    pre_nrem_length = [];
    post_nrem_length = [];

    pre_rem_length = [];
    post_rem_length = [];

    pre_wake_length = [];
    post_wake_length = [];

    pre_total_length = [];
    post_total_length = [];

    for animal = 1:length(animal_id)

        for day = animal_days{animal}

            load([data_path '\' animal_id{animal} '\sleep_activity\day_'  num2str(day) '_block_1_sleep_activity.mat']);
            pre_nrem_length = [pre_nrem_length (sum(beh_state.nrem==1)*10)/60];
            pre_rem_length = [pre_rem_length (sum(beh_state.rem==1)*10)/60];
            pre_wake_length = [pre_wake_length (sum(beh_state.wake==1)*10)/60];

            if animal == 1 && day == 6
                load([data_path '\' animal_id{animal} '\sleep_activity\day_'  num2str(day) '_block_3_sleep_activity.mat']);
            else
                load([data_path '\' animal_id{animal} '\sleep_activity\day_'  num2str(day) '_block_2_sleep_activity.mat']);
            end
            post_nrem_length = [post_nrem_length (sum(beh_state.nrem==1)*10)/60];
            post_rem_length = [post_rem_length (sum(beh_state.rem==1)*10)/60];
            post_wake_length = [post_wake_length (sum(beh_state.wake==1)*10)/60];


        end
    end

    fig = figure;
    hold on;

        scatter(ones(1,length(pre_nrem_length))-((rand(1,55)-0.5)*.15),pre_nrem_length,30,[.5 .5 .5],'filled');
        errorbar(1,mean(pre_nrem_length),std(pre_nrem_length)/sqrt(length(pre_nrem_length)),'color','r','LineWidth',2)

        scatter(2*ones(1,length(post_nrem_length))-((rand(1,55)-0.5)*.15),post_nrem_length,30,[.5 .5 .5],'filled');
        errorbar(2,mean(post_nrem_length),std(post_nrem_length)/sqrt(length(post_nrem_length)),'color','r','LineWidth',2)

        scatter(4*ones(1,length(pre_rem_length))-((rand(1,55)-0.5)*.15),pre_rem_length,30,[.5 .5 .5],'filled');
        errorbar(4,mean(pre_rem_length),std(pre_rem_length)/sqrt(length(pre_rem_length)),'color','r','LineWidth',2)

        scatter(5*ones(1,length(post_rem_length))-((rand(1,55)-0.5)*.15),post_rem_length,30,[.5 .5 .5],'filled');
        errorbar(5,mean(post_rem_length),std(post_rem_length)/sqrt(length(post_rem_length)),'color','r','LineWidth',2)

        scatter(7*ones(1,length(pre_wake_length))-((rand(1,55)-0.5)*.15),pre_wake_length,30,[.5 .5 .5],'filled');
        errorbar(7,mean(pre_wake_length),std(pre_wake_length)/sqrt(length(pre_wake_length)),'color','r','LineWidth',2)

        scatter(8*ones(1,length(post_wake_length))-((rand(1,55)-0.5)*.15),post_wake_length,30,[.5 .5 .5],'filled');
        errorbar(8,mean(post_wake_length),std(post_wake_length)/sqrt(length(post_wake_length)),'color','r','LineWidth',2)

        xlim([0 9])
        ylim([0 200])

    [h p ci stats] = ttest(pre_nrem_length,post_nrem_length)
        mean(pre_nrem_length)
        std(pre_nrem_length)/sqrt(length(pre_nrem_length))
        mean(post_nrem_length)
        std(post_nrem_length)/sqrt(length(post_nrem_length))

    [h p ci stats] = ttest(pre_rem_length,post_rem_length)
        mean(pre_rem_length)
        std(pre_rem_length)/sqrt(length(pre_rem_length))
        mean(post_rem_length)
        std(post_rem_length)/sqrt(length(post_rem_length))

    [h p ci stats] = ttest(pre_wake_length,post_wake_length)
        mean(pre_wake_length)
        std(pre_wake_length)/sqrt(length(pre_wake_length))
        mean(post_wake_length)
        std(post_wake_length)/sqrt(length(post_wake_length))

    mean(pre_nrem_length + pre_rem_length + pre_wake_length)
    std(pre_nrem_length + pre_rem_length + pre_wake_length)/sqrt(length(pre_nrem_length))
    mean(post_nrem_length + post_rem_length + post_wake_length)
    std(post_nrem_length + post_rem_length + post_wake_length)/sqrt(length(post_nrem_length))

%% sFIG 2 | COMPARE DELTA POWER FIRST HALF SECOND HALF 

    animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
    animal_days = {[1:8],[1:6],[2:10],[4:17],[4:13],[1:5 7:9]};

    pre_first_half_delta_interp = [];
    pre_second_half_delta_interp = [];
    post_first_half_delta_interp = [];
    post_second_half_delta_interp = [];

    for animal = 1:length(animal_id)

        for day = animal_days{animal}

            load([data_path '\' animal_id{animal} '\sleep_activity\day_'  num2str(day) '_block_1_sleep_activity.mat']);
            tmp_delta = beh_state.delta_power(beh_state.nrem==1);
            pre_first_half_delta_interp = [pre_first_half_delta_interp; interp1(1:length(tmp_delta(1:floor(size(tmp_delta,1)/2))),tmp_delta(1:floor(size(tmp_delta,1)/2)),linspace(1,length(tmp_delta(1:floor(size(tmp_delta,1)/2))),1000))];
            pre_second_half_delta_interp = [pre_second_half_delta_interp; interp1(1:length(tmp_delta(1+floor(size(tmp_delta,1)/2):end)),tmp_delta(1+floor(size(tmp_delta,1)/2):end),linspace(1,length(tmp_delta(1+floor(size(tmp_delta,1)/2):end)),1000))];


            if animal == 1 && day == 6
                load([data_path '\' animal_id{animal} '\sleep_activity\day_'  num2str(day) '_block_3_sleep_activity.mat']);
            else
                load([data_path '\' animal_id{animal} '\sleep_activity\day_'  num2str(day) '_block_2_sleep_activity.mat']);
            end
            tmp_delta = beh_state.delta_power(beh_state.nrem==1);
            post_first_half_delta_interp = [post_first_half_delta_interp; interp1(1:length(tmp_delta(1:floor(size(tmp_delta,1)/2))),tmp_delta(1:floor(size(tmp_delta,1)/2)),linspace(1,length(tmp_delta(1:floor(size(tmp_delta,1)/2))),1000))];
            post_second_half_delta_interp = [post_second_half_delta_interp; interp1(1:length(tmp_delta(1+floor(size(tmp_delta,1)/2):end)),tmp_delta(1+floor(size(tmp_delta,1)/2):end),linspace(1,length(tmp_delta(1+floor(size(tmp_delta,1)/2):end)),1000))];

        end
    end

    fig = figure;
        subplot(1,4,1); hold on;
            plot(mean(pre_first_half_delta_interp)+std(pre_first_half_delta_interp)/sqrt(size(pre_first_half_delta_interp,1)),'k')
            plot(mean(pre_first_half_delta_interp)-std(pre_first_half_delta_interp)/sqrt(size(pre_first_half_delta_interp,1)),'k')
            ylim([.5 1.8])
        subplot(1,4,2); hold on;
            plot(mean(pre_second_half_delta_interp)+std(pre_second_half_delta_interp)/sqrt(size(pre_second_half_delta_interp,1)),'k')
            plot(mean(pre_second_half_delta_interp)-std(pre_second_half_delta_interp)/sqrt(size(pre_second_half_delta_interp,1)),'k')
            ylim([.5 1.8])
        subplot(1,4,3); hold on;
            plot(mean(post_first_half_delta_interp)+std(post_first_half_delta_interp)/sqrt(size(post_first_half_delta_interp,1)),'k')
            plot(mean(post_first_half_delta_interp)-std(post_first_half_delta_interp)/sqrt(size(post_first_half_delta_interp,1)),'k')
            ylim([.5 1.8])
        subplot(1,4,4); hold on;
            plot(mean(post_second_half_delta_interp)+std(post_second_half_delta_interp)/sqrt(size(post_second_half_delta_interp,1)),'k')
            plot(mean(post_second_half_delta_interp)-std(post_second_half_delta_interp)/sqrt(size(post_second_half_delta_interp,1)),'k')
            ylim([.5 1.8])

    fig = figure;
        hold on;
        errorbar(1,mean(mean(pre_first_half_delta_interp')),std(mean(pre_first_half_delta_interp'))/sqrt(length(mean(pre_first_half_delta_interp'))),'color','k')
        errorbar(2,mean(mean(pre_second_half_delta_interp')),std(mean(pre_second_half_delta_interp'))/sqrt(length(mean(pre_second_half_delta_interp'))),'color','k')
        errorbar(3,mean(mean(post_first_half_delta_interp')),std(mean(post_first_half_delta_interp'))/sqrt(length(mean(post_first_half_delta_interp'))),'color','k')
        errorbar(4,mean(mean(post_second_half_delta_interp')),std(mean(post_second_half_delta_interp'))/sqrt(length(mean(post_second_half_delta_interp'))),'color','k')
        ylim([.8 1.4])
        xlim([0 5])

        [h p] = ttest2(mean(pre_first_half_delta_interp'),mean(pre_second_half_delta_interp'))
        [h p] = ttest2(mean(pre_first_half_delta_interp'),mean(post_first_half_delta_interp'))
        [h p] = ttest2(mean(pre_first_half_delta_interp'),mean(post_second_half_delta_interp'))

        [h p] = ttest2(mean(pre_second_half_delta_interp'),mean(post_first_half_delta_interp'))
        [h p] = ttest2(mean(pre_second_half_delta_interp'),mean(post_second_half_delta_interp'))
        [h p] = ttest2(mean(post_first_half_delta_interp'),mean(post_second_half_delta_interp'))
