%% SCRIPT TO GENERATE FIGURE 2 & FIGURE 2 FIGURE SUPPLEMENTS
% 
% to generate figures, download data from https://doi.org/10.7272/Q6KK9927
%
% HOWEVER, several plots for figure 2 were generated from large raw data 
% files that are not uploaded to Dryad - I hope to generate and upload 
% intermediate processed data to allow generation of all panels of this 
% figure in the future. Code can be viewed below but only sFig 2 will be 
% generated (SL - 8/23/21)

clc; clear all; close all;
data_path = 'D:\Lemke_etal_2021\Lemke_eLife_2021_data';

%% NEED RAW DATA NOT ON DRYAD | FIG 2 A, B, C, & D | EXAMPLE DLS UNIT-M1 4-8Hz LFP PHASE LOCKING & EXAMPLE M1-DLS 4-8HZ LFP COH

%%% EXAMPLE

    %%% T200 

    day_blocks = {[1],[2 7 8],[9 11 12 13 14],[15 16 17 18 19],[20 21 22 23 24],[25 26 27],[28 29 30],[32 33 34 35 36],[37 38 39 40 41],[42 43 44 45 46 47]};
    pre_sleep_blocks = {[1],[1],[2],[1],[1],[1],[1],[1],[1],[1]};
    post_sleep_blocks = {[],[3],[5],[5],[4],[3],[3],[5],[5],[6]};
    plx_files = {'T200_blocks_1-01','T200_blocks_2_7_8-01','T200_blocks_9_11_12_13_14-01','T200_blocks_15_16_17_18_19-01','T200_blocks_20_21_22_23_24-01','T200_blocks_25_26_27-01','T200_blocks_28_29_30-01','T200_blocks_32_33_34_35_36-01','T200_blocks_37_38_39_40_41-01','T200_blocks_42_43_44_45_46_47-01'};

    load([data_path '\data\T200\sleep_activity\LFP_coherence.mat'])

    for day = 2%:10

        tmp_day_blocks = day_blocks{day};

        %%% LOAD BEH STATE, LFP, AND START/STOP

            load([data_path '\data\T200\sleep_activity\day_' num2str(day) '_block_1_sleep_activity.mat'])
            pre_beh_state = beh_state;
            data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\Block-' num2str(tmp_day_blocks(pre_sleep_blocks{day}))],'Type',4);
            pre_lfp = data.streams.LFPs;
            pre_video_start_stop_lfp = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
            pre_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];
            pre_lfp = pre_lfp.data;

            pre_lfp_norm = [];
            for n = 1:32
                tmp_lfp = pre_lfp(n,:);
                tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                b = find(abs(tmp_lfp)>6);
                tmp_lfp = pre_lfp(n,:);
                tmp_lfp(b)=NaN;
                tmp_mean=nanmean(tmp_lfp);
                tmp_sd=nanstd(tmp_lfp);
                tmp_lfp(b)=tmp_mean;
                pre_lfp_norm =  [pre_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
            end
            pre_lfp_norm(1:16,:) = pre_lfp_norm(1:16,:)-repmat(mean(pre_lfp_norm(1:16,:),1),[size(pre_lfp_norm(1:16,:),1) 1]);
            pre_lfp_norm(17:32,:) = pre_lfp_norm(17:32,:)-repmat(mean(pre_lfp_norm(17:32,:),1),[size(pre_lfp_norm(17:32,:),1) 1]);

            load(['D:\Lemke_etal_2021\data\T200\sleep_activity\day_' num2str(day) '_block_2_sleep_activity.mat'])
            post_beh_state = beh_state;
            data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\Block-' num2str(tmp_day_blocks(post_sleep_blocks{day}))],'Type',4);
            post_lfp = data.streams.LFPs;
            post_video_start_stop_lfp = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
            post_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];
            post_lfp = post_lfp.data;

            post_lfp_norm = [];
            for n = 1:32
                tmp_lfp = post_lfp(n,:);
                tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                b = find(abs(tmp_lfp)>6);
                tmp_lfp = post_lfp(n,:);
                tmp_lfp(b)=NaN;
                tmp_mean=nanmean(tmp_lfp);
                tmp_sd=nanstd(tmp_lfp);
                tmp_lfp(b)=tmp_mean;
                post_lfp_norm =  [post_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
            end
            post_lfp_norm(1:16,:) = post_lfp_norm(1:16,:)-repmat(mean(post_lfp_norm(1:16,:),1),[size(post_lfp_norm(1:16,:),1) 1]);
            post_lfp_norm(17:32,:) = post_lfp_norm(17:32,:)-repmat(mean(post_lfp_norm(17:32,:),1),[size(post_lfp_norm(17:32,:),1) 1]);

        %%% LOAD SPIKING DATA

            plx_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\' plx_files{day} '.txt']);
            plx_spiking = table2array(plx_spiking);

            tmp_ts = cell(32,length(tmp_day_blocks));
            block_count = 1;
            for block = tmp_day_blocks
                tmp_num = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\Block-' num2str(block)],'Type',3);
                for chan = 1:32
                    tmp_ts{chan,block_count} = tmp_num.snips.eNeu.ts(tmp_num.snips.eNeu.chan==chan);
                end
                block_count = block_count + 1;
            end

            tmp_length = zeros(32,length(tmp_day_blocks));
            for chan = 1:32
                for block = 1:length(tmp_day_blocks)
                    tmp_length(chan,block) = length(tmp_ts{chan,block});
                end
            end
            tmp_length = [zeros(32,1) tmp_length];

            tmp_sort = cell(32,length(tmp_day_blocks));
            for chan = 1:32
                chan_sorts = plx_spiking(plx_spiking(:,1)==chan,:);
                for block = 1:length(tmp_day_blocks)
                    tmp_sort{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),2);
                end
            end

            pre_dls_spiking = cell(16,max(plx_spiking(:,2)));
            post_dls_spiking = cell(16,max(plx_spiking(:,2)));
            count = 1;
            for dls_chan = 1:16
                for dls_unit = 1:max(plx_spiking(:,2))
                    pre_dls_spiking{count,dls_unit} = tmp_ts{dls_chan,pre_sleep_blocks{day}}(tmp_sort{dls_chan,pre_sleep_blocks{day}}==dls_unit);
                    post_dls_spiking{count,dls_unit} = tmp_ts{dls_chan,post_sleep_blocks{day}}(tmp_sort{dls_chan,post_sleep_blocks{day}}==dls_unit);
                end
                count = count + 1;
            end

        %%% EXAMPLE PHASE LOCKING HISTOGRAM

            for m1_lfp = 18%:32

                tmp_pre_m1_lfp = pre_lfp_norm(m1_lfp,:);
                tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2));
                tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
                tmp_pre_m1_lfp = eegfilt(tmp_pre_m1_lfp,data.streams.LFPs.fs,4,8);
                tmp_pre_m1_lfp = angle(hilbert(tmp_pre_m1_lfp));

                tmp_post_m1_lfp = post_lfp_norm(m1_lfp,:);
                tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2));
                tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
                tmp_post_m1_lfp = eegfilt(tmp_post_m1_lfp,data.streams.LFPs.fs,4,8);
                tmp_post_m1_lfp = angle(hilbert(tmp_post_m1_lfp));

                for dls_chan = 2%1:size(pre_dls_spiking,1)
                    for dls_unit = 1%1:size(pre_dls_spiking,2)

                        if isempty(pre_dls_spiking{dls_chan,dls_unit}) || isempty(post_dls_spiking{dls_chan,dls_unit})
                            continue
                        end

                        pre_dls_unit = pre_dls_spiking{dls_chan,dls_unit};%(pre_dls_spiking{dls_chan,dls_unit}>pre_video_start_stop(1) & pre_dls_spiking{dls_chan,dls_unit}<pre_video_start_stop(2))-pre_video_start_stop(1);
                        pre_dls_unit = histc(pre_dls_unit,linspace(pre_video_start_stop(1),pre_video_start_stop(2),length(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2))));
                        pre_dls_unit = pre_dls_unit(pre_beh_state.nrem_interp==1);

                        post_dls_unit = post_dls_spiking{dls_chan,dls_unit};%(post_dls_spiking{dls_chan,dls_unit}>post_video_start_stop(1) & post_dls_spiking{dls_chan,dls_unit}<post_video_start_stop(2))-post_video_start_stop(1);
                        post_dls_unit = histc(post_dls_unit,linspace(post_video_start_stop(1),post_video_start_stop(2),length(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2))));
                        post_dls_unit = post_dls_unit(post_beh_state.nrem_interp==1);

                        figure;
                            histogram(tmp_post_m1_lfp(post_dls_unit==1),[-pi:pi/8:pi],'normalization','probability','DisplayStyle','Stairs')
                            ylim([0.04 0.08])
                            [pval, z] = circ_rtest(tmp_post_m1_lfp(post_dls_unit==1));
                            sd = circ_std(tmp_post_m1_lfp(post_dls_unit==1)');
                            title(['CHAN ' num2str(dls_chan) ' UNIT ' num2str(dls_unit) ' P=' num2str(pval) ' SD=' num2str(sd)])

                    end
                end        
            end   

        %%% EXAMPLE PHASE LOCKING TIME SERIES

            scale = 1:2048736;

            for m1_lfp = 18%

                tmp_pre_m1_lfp = pre_lfp_norm(m1_lfp,:);
                tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2));
                tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
                tmp_pre_m1_lfp_filt = eegfilt(tmp_pre_m1_lfp,data.streams.LFPs.fs,4,8);
                tmp_pre_m1_lfp_phase = angle(hilbert(tmp_pre_m1_lfp_filt));

                tmp_post_m1_lfp = post_lfp_norm(m1_lfp,:);
                tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2));
                tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
                tmp_post_m1_lfp_filt = eegfilt(tmp_post_m1_lfp,data.streams.LFPs.fs,4,8);
                tmp_post_m1_lfp_phase = angle(hilbert(tmp_post_m1_lfp_filt));

                figure;

                ax1 = subplot(2,1,1); hold on;
                plot(zscore(tmp_post_m1_lfp(scale)),'color',[.5 .5 .5])
                plot(zscore(tmp_post_m1_lfp_filt(scale)),'color',[0 0 0])

                ax2 = subplot(2,1,2); hold on;
                plot(tmp_post_m1_lfp_phase(scale),'color','k')

                for dls_chan = 2
                    for dls_unit = 1

                        pre_dls_unit = pre_dls_spiking{dls_chan,dls_unit};%(pre_dls_spiking{dls_chan,dls_unit}>pre_video_start_stop(1) & pre_dls_spiking{dls_chan,dls_unit}<pre_video_start_stop(2))-pre_video_start_stop(1);
                        pre_dls_unit = histc(pre_dls_unit,linspace(pre_video_start_stop(1),pre_video_start_stop(2),length(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2))));
                        pre_dls_unit = pre_dls_unit(pre_beh_state.nrem_interp==1);

                        post_dls_unit = post_dls_spiking{dls_chan,dls_unit};%(post_dls_spiking{dls_chan,dls_unit}>post_video_start_stop(1) & post_dls_spiking{dls_chan,dls_unit}<post_video_start_stop(2))-post_video_start_stop(1);
                        post_dls_unit = histc(post_dls_unit,linspace(post_video_start_stop(1),post_video_start_stop(2),length(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2))));
                        post_dls_unit = post_dls_unit(post_beh_state.nrem_interp==1);

                        post_dls_unit = find(post_dls_unit((scale)));

                        for spike = 1:length(post_dls_unit)
                            plot([post_dls_unit(spike) post_dls_unit(spike)],[-2 2],'color','r');
                        end

                        linkaxes([ax1 ax2],'x')

                    end
                end        
            end 

        %%% M1 CHAN LFP COH VS. DLS UNIT PLV

            pre_coh = [];
            post_coh = [];
            for pairs=1:256
                pre_coh = [pre_coh mean(pre_LFP_coh{pairs}(day-1,48:113))];
                post_coh = [post_coh mean(post_LFP_coh{pairs}(day-1,48:113))];
            end
            pre_coh = mean(reshape(pre_coh,16,16));
            post_coh = mean(reshape(post_coh,16,16));

            dls_plv = cell(1,16);
            for m1_lfp = 17:32

                tmp_pre_m1_lfp = pre_lfp_norm(m1_lfp,:);
                tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2));
                tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
                tmp_pre_m1_lfp = eegfilt(tmp_pre_m1_lfp,data.streams.LFPs.fs,4,8);
                tmp_pre_m1_lfp = angle(hilbert(tmp_pre_m1_lfp));

                tmp_post_m1_lfp = post_lfp_norm(m1_lfp,:);
                tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2));
                tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
                tmp_post_m1_lfp = eegfilt(tmp_post_m1_lfp,data.streams.LFPs.fs,4,8);
                tmp_post_m1_lfp = angle(hilbert(tmp_post_m1_lfp));

                for dls_chan = 1:size(pre_dls_spiking,1)
                    for dls_unit = 1:size(pre_dls_spiking,2)

                        if isempty(pre_dls_spiking{dls_chan,dls_unit}) || isempty(post_dls_spiking{dls_chan,dls_unit})
                            continue
                        end

                        pre_dls_unit = pre_dls_spiking{dls_chan,dls_unit};%(pre_dls_spiking{dls_chan,dls_unit}>pre_video_start_stop(1) & pre_dls_spiking{dls_chan,dls_unit}<pre_video_start_stop(2))-pre_video_start_stop(1);
                        pre_dls_unit = histc(pre_dls_unit,linspace(pre_video_start_stop(1),pre_video_start_stop(2),length(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2))));
                        pre_dls_unit = pre_dls_unit(pre_beh_state.nrem_interp==1);     

                        post_dls_unit = post_dls_spiking{dls_chan,dls_unit};%(post_dls_spiking{dls_chan,dls_unit}>post_video_start_stop(1) & post_dls_spiking{dls_chan,dls_unit}<post_video_start_stop(2))-post_video_start_stop(1);
                        post_dls_unit = histc(post_dls_unit,linspace(post_video_start_stop(1),post_video_start_stop(2),length(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2))));
                        post_dls_unit = post_dls_unit(post_beh_state.nrem_interp==1);

                        sd = circ_std([tmp_pre_m1_lfp(pre_dls_unit==1) tmp_post_m1_lfp(post_dls_unit==1)]');
                        dls_plv{m1_lfp-16} = [dls_plv{m1_lfp-16} sd];

                    end
                end
            end

            [m, i] = sort(pre_coh);
            figure;
                subplot(2,1,1)
                    plot(pre_coh(i),'color','k')
                subplot(2,1,2); hold on;
                    for lfp_chan = 1:16
                        scatter(lfp_chan*ones(1,length(dls_plv{i(lfp_chan)})),dls_plv{i(lfp_chan)},30,[0 0 0],'filled')
                        errorbar(lfp_chan,mean(dls_plv{i(lfp_chan)}),std(dls_plv{i(lfp_chan)})/sqrt(length(dls_plv{i(lfp_chan)})),'color','r');
                    end    

        %%% M1 & DLS CHAN LFP COH

            scale = 1:2048736;

            for m1_lfp = 18%:16
                for dls_lfp = 1%:16

                    tmp_post_m1_lfp = post_lfp_norm(m1_lfp,:);
                    tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2));
                    tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
                    tmp_post_m1_lfp_filt = eegfilt(tmp_post_m1_lfp,data.streams.LFPs.fs,4,8);
                    tmp_post_m1_lfp_phase = angle(hilbert(tmp_post_m1_lfp_filt));

                    tmp_post_dls_lfp = post_lfp_norm(dls_lfp,:);
                    tmp_post_dls_lfp = tmp_post_dls_lfp(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2));
                    tmp_post_dls_lfp = tmp_post_dls_lfp(post_beh_state.nrem_interp==1);
                    tmp_post_dls_lfp_filt = eegfilt(tmp_post_dls_lfp,data.streams.LFPs.fs,4,8);
                    tmp_post_dls_lfp_phase = angle(hilbert(tmp_post_dls_lfp_filt));

                    figure;
                        ax1 = subplot(3,1,1); hold on;
                            plot(tmp_post_m1_lfp(scale),'color',[1 0 0])
                            plot(tmp_post_dls_lfp(scale),'color',[0 0 1])
                            xlim([1 length(scale)])
                        ax2 = subplot(3,1,2); hold on;
                            plot(tmp_post_m1_lfp_filt(scale),'color',[1 0 0])
                            plot(tmp_post_dls_lfp_filt(scale),'color',[0 0 1])
                            xlim([1 length(scale)])                        
                        ax3 = subplot(3,1,3); hold on;
                            plot(tmp_post_m1_lfp_phase(scale),'color',[1 0 0])
                            plot(tmp_post_dls_lfp_phase(scale),'color',[0 0 1])
                        linkaxes([ax1 ax2 ax3],'x')

                    figure;
                        histogram(tmp_post_m1_lfp_phase-tmp_post_dls_lfp_phase,[-pi:pi/8:pi],'normalization','probability','DisplayStyle','Stairs')
                        ylim([0 0.1])
                end
            end

    end
    
%% NEED RAW DATA NOT ON DRYAD | FIG 2 E & F | ALL ANIMAL DLS UNIT M1 PHASE LOCKING DIFFERENCE 
   
all_unit_coh_lfp_coh = [];
all_unit_coh_lfp_coh_day_norm = [];

%%% T102

        animal_unit_lfp_coh = [];
        animal_unit_lfp_coh_norm = [];

        day_blocks = {[2 3 7 8 9 11 12 13], ...
            [14 16 18 19 20 21], ...
            [22 24 25 26 27 29], ...
            [30 31 32 33 34 35 36 37 38], ...
            [39 43 44 45 47], ...
            [48 49 51 53 54 56], ...
            [57 58 59 60 61], ...
            [62 63 64 65 66]};
        pre_sleep_blocks = {[2],[1],[1],[1],[1],[1],[1],[1]};
        reach_blocks = {[3 4 6 7], [2 3 5],[2 3 6],[2 3 5 6 7 9],[2 3 5],[2 3 5],[2 3 5]};
        post_sleep_blocks = {[8],[4],[5],[8],[4],[5],[4],[4]};
        plx_files = {'T102_blocks_2_3_7_8_9_11_12_13-01','T102_blocks_14_16_18_19_20_21-01','T102_blocks_22_24_25_26_27_29-01','T102_blocks_30_31_32_33_34_35_36_37_38-01','T102_blocks_39_43_44_45_47-01','T102_blocks_48_49_51_53_54_56-01','T102_blocks_57_58_59_60_61-01','T102_blocks_62_63_64_65_66-01'};

        load(['D:\Lemke_etal_2021\data\T102\sleep_activity\LFP_coherence.mat'])

        for day = 1:8

            tmp_day_blocks = day_blocks{day};

            %%% PRE AND POST SLEEP NREM COH

            pre_coh = [];
            post_coh = [];
            for pairs=1:256
                pre_coh = [pre_coh mean(pre_LFP_coh{pairs}(day,48:113))];
                post_coh = [post_coh mean(post_LFP_coh{pairs}(day,48:113))];
            end
            pre_coh = mean(reshape(pre_coh,16,16));
            post_coh = mean(reshape(post_coh,16,16));

            %%% LOAD BEH STATE, LFP, AND START/STOP

            load(['D:\Lemke_etal_2021\data\T102\sleep_activity\day_' num2str(day) '_block_1_sleep_activity.mat'])
            pre_beh_state = beh_state;
            data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\Block-' num2str(tmp_day_blocks(pre_sleep_blocks{day}))],'Type',4);
            pre_lfp = data.streams.LFPs;
            pre_video_start_stop_lfp = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
            pre_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];
            pre_lfp = pre_lfp.data;

            pre_lfp_norm = [];
            for n = 1:16
                tmp_lfp = pre_lfp(n,:);
                tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                b = find(abs(tmp_lfp)>6);
                tmp_lfp = pre_lfp(n,:);
                tmp_lfp(b)=NaN;
                tmp_mean=nanmean(tmp_lfp);
                tmp_sd=nanstd(tmp_lfp);
                tmp_lfp(b)=tmp_mean;
                pre_lfp_norm =  [pre_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
            end
            pre_lfp_norm(1:16,:) = pre_lfp_norm(1:16,:)-repmat(mean(pre_lfp_norm(1:16,:),1),[size(pre_lfp_norm(1:16,:),1) 1]);

            if day==1 | day==4
                load(['D:\Lemke_etal_2021\data\T102\sleep_activity\day_' num2str(day) '_block_3_sleep_activity.mat'])
            else
                load(['D:\Lemke_etal_2021\data\T102\sleep_activity\day_' num2str(day) '_block_2_sleep_activity.mat'])
            end
            post_beh_state = beh_state;
            data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\Block-' num2str(tmp_day_blocks(post_sleep_blocks{day}))],'Type',4);
            post_lfp = data.streams.LFPs;
            post_video_start_stop_lfp = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
            post_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];
            post_lfp = post_lfp.data;

            post_lfp_norm = [];
            for n = 1:16
                tmp_lfp = post_lfp(n,:);
                tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                b = find(abs(tmp_lfp)>6);
                tmp_lfp = post_lfp(n,:);
                tmp_lfp(b)=NaN;
                tmp_mean=nanmean(tmp_lfp);
                tmp_sd=nanstd(tmp_lfp);
                tmp_lfp(b)=tmp_mean;
                post_lfp_norm =  [post_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
            end
            post_lfp_norm(1:16,:) = post_lfp_norm(1:16,:)-repmat(mean(post_lfp_norm(1:16,:),1),[size(post_lfp_norm(1:16,:),1) 1]);

            %%% LOAD SPIKING DATA

            plx_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\' plx_files{day} '.txt']);
            plx_spiking = table2array(plx_spiking);

            tmp_ts = cell(32,length(tmp_day_blocks));
            block_count = 1;
            for block = tmp_day_blocks
                tmp_num = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\Block-' num2str(block)],'Type',3);
                for chan = 1:32
                    tmp_ts{chan,block_count} = tmp_num.snips.eNeu.ts(tmp_num.snips.eNeu.chan==chan);
                end
                block_count = block_count + 1;
            end

            tmp_length = zeros(32,length(tmp_day_blocks));
            for chan = 1:32
                for block = 1:length(tmp_day_blocks)
                    tmp_length(chan,block) = length(tmp_ts{chan,block});
                end
            end
            tmp_length = [zeros(32,1) tmp_length];

            tmp_sort = cell(32,length(tmp_day_blocks));
            tmp_width = cell(32,length(tmp_day_blocks));
            tmp_peak2valley = cell(32,length(tmp_day_blocks));
            for chan = 1:32
                chan_sorts = plx_spiking(plx_spiking(:,1)==chan,:);
                for block = 1:length(tmp_day_blocks)
                    tmp_sort{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),2);
                    tmp_width{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),6);
                    tmp_peak2valley{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),4);
                end
            end

            pre_dls_spiking = cell(16,max(plx_spiking(:,2)));
            post_dls_spiking = cell(16,max(plx_spiking(:,2)));
            count = 1;
            for dls_chan = 17:32
                for dls_unit = 1:max(plx_spiking(:,2))
                    pre_dls_spiking{count,dls_unit} = tmp_ts{dls_chan,pre_sleep_blocks{day}}(tmp_sort{dls_chan,pre_sleep_blocks{day}}==dls_unit);
                    post_dls_spiking{count,dls_unit} = tmp_ts{dls_chan,post_sleep_blocks{day}}(tmp_sort{dls_chan,post_sleep_blocks{day}}==dls_unit);
                end
                count = count + 1;
            end

            %%% CALCULATE SPIKE COH

            coh_coh = [];

            for m1_lfp = 1:16

                tmp_pre_m1_lfp = pre_lfp_norm(m1_lfp,:);
                tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2));
                tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
                tmp_pre_m1_lfp = eegfilt(tmp_pre_m1_lfp,data.streams.LFPs.fs,4,8);
                tmp_pre_m1_lfp = angle(hilbert(tmp_pre_m1_lfp));

                tmp_post_m1_lfp = post_lfp_norm(m1_lfp,:);
                tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2));
                tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
                tmp_post_m1_lfp = eegfilt(tmp_post_m1_lfp,data.streams.LFPs.fs,4,8);
                tmp_post_m1_lfp = angle(hilbert(tmp_post_m1_lfp));

                tmp_chan_plv = [];

                for dls_chan = 1:size(pre_dls_spiking,1)
                    for dls_unit = 1:size(pre_dls_spiking,2)

                        if isempty(pre_dls_spiking{dls_chan,dls_unit}) || isempty(post_dls_spiking{dls_chan,dls_unit})
                            continue
                        end

                        pre_dls_unit = pre_dls_spiking{dls_chan,dls_unit};%(pre_dls_spiking{dls_chan,dls_unit}>pre_video_start_stop(1) & pre_dls_spiking{dls_chan,dls_unit}<pre_video_start_stop(2))-pre_video_start_stop(1);
                        pre_dls_unit = histc(pre_dls_unit,linspace(pre_video_start_stop(1),pre_video_start_stop(2),length(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2))));
                        pre_dls_unit = pre_dls_unit(pre_beh_state.nrem_interp==1);

                        post_dls_unit = post_dls_spiking{dls_chan,dls_unit};%(post_dls_spiking{dls_chan,dls_unit}>post_video_start_stop(1) & post_dls_spiking{dls_chan,dls_unit}<post_video_start_stop(2))-post_video_start_stop(1);
                        post_dls_unit = histc(post_dls_unit,linspace(post_video_start_stop(1),post_video_start_stop(2),length(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2))));
                        post_dls_unit = post_dls_unit(post_beh_state.nrem_interp==1);

                        [pval, z] = circ_rtest([tmp_pre_m1_lfp(pre_dls_unit==1) tmp_post_m1_lfp(post_dls_unit==1)]);
                        sd = circ_std([tmp_pre_m1_lfp(pre_dls_unit==1) tmp_post_m1_lfp(post_dls_unit==1)]');                        
                        tmp_chan_plv = [tmp_chan_plv; pval z sd pre_coh(m1_lfp)];

                    end
                end

                coh_coh = [coh_coh; mean(tmp_chan_plv)];

            end

            animal_unit_lfp_coh = [animal_unit_lfp_coh; coh_coh];
            animal_unit_lfp_coh_norm = [animal_unit_lfp_coh_norm; coh_coh(:,1) zscore(coh_coh(:,2)) zscore(coh_coh(:,3)) zscore(coh_coh(:,4))]; 

        end

        all_unit_coh_lfp_coh = [all_unit_coh_lfp_coh; animal_unit_lfp_coh];
        all_unit_coh_lfp_coh_day_norm = [all_unit_coh_lfp_coh_day_norm; animal_unit_lfp_coh_norm];
        
        figure;
            subplot(2,2,1); hold on;
                scatter(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,2),30,[0 0 0],'filled')
                [R,P] = corrcoef(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,2));
                [p] = polyfit(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,2),1);
                x=linspace(min(animal_unit_lfp_coh(:,4)),max(animal_unit_lfp_coh(:,4)),10);
                xlim([min(animal_unit_lfp_coh(:,4)) max(animal_unit_lfp_coh(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T102 - lfp coh vs. dls unit m1 lfp z val - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,2); hold on;
                scatter(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,2),30,[0 0 0],'filled');
                [R,P] = corrcoef(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,2));
                [p] = polyfit(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,2),1);
                x=linspace(min(animal_unit_lfp_coh_norm(:,4)),max(animal_unit_lfp_coh_norm(:,4)),10);
                xlim([min(animal_unit_lfp_coh_norm(:,4)) max(animal_unit_lfp_coh_norm(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T102 - norm lfp coh vs. norm dls unit m1 lfp z val - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,3); hold on;
                scatter(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,3),30,[0 0 0],'filled')
                [R,P] = corrcoef(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,3));
                [p] = polyfit(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,3),1);
                x=linspace(min(animal_unit_lfp_coh(:,4)),max(animal_unit_lfp_coh(:,4)),10);
                xlim([min(animal_unit_lfp_coh(:,4)) max(animal_unit_lfp_coh(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T102 - lfp coh vs. dls unit m1 lfp SD - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,4); hold on;
                scatter(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,3),30,[0 0 0],'filled');
                [R,P] = corrcoef(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,3));
                [p] = polyfit(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,3),1);
                x=linspace(min(animal_unit_lfp_coh_norm(:,4)),max(animal_unit_lfp_coh_norm(:,4)),10);
                xlim([min(animal_unit_lfp_coh_norm(:,4)) max(animal_unit_lfp_coh_norm(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T102 - norm lfp coh vs. norm dls unit m1 lfp SD - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
                
%%% T107

        animal_unit_lfp_coh = [];
        animal_unit_lfp_coh_norm = [];

        day_blocks = {[1 2],[3 4 5 6 8],[9 10 11 12 13 14 15],[16 17 18 19 20 21 22 24 25],[26 27 30 31 36 39],[41 42 43 47 49 51 52]};
        pre_sleep_blocks = {[1],[1],[1],[3],[1],[1]};
        post_sleep_blocks = {[2],[4],[6],[8],[5],[4]};
        plx_files = {'T107_blocks_1_2-01','T107_blocks_3_4_5_6_8-01','T107_blocks_9_10_11_12_13_14_15-01','T107_blocks_16_17_18_19_20_21_22_24_25-01','T107_blocks_26_27_30_31_36_39-01','T107_blocks_41_42_43_47_49_51_52-01'};

        load(['D:\Lemke_etal_2021\data\T107\sleep_activity\LFP_coherence.mat'])

        for day = 1:6

            tmp_day_blocks = day_blocks{day};

            %%% PRE AND POST SLEEP NREM COH

            pre_coh = [];
            post_coh = [];
            for pairs=1:256
                pre_coh = [pre_coh mean(pre_LFP_coh{pairs}(day,48:113))];
                post_coh = [post_coh mean(post_LFP_coh{pairs}(day,48:113))];
            end
            pre_coh = mean(reshape(pre_coh,16,16));
            post_coh = mean(reshape(post_coh,16,16));

            %%% LOAD BEH STATE, LFP, AND START/STOP

            load(['D:\Lemke_etal_2021\data\T107\sleep_activity\day_' num2str(day) '_block_1_sleep_activity.mat'])
            pre_beh_state = beh_state;
            data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T107\neural_data\Block-' num2str(tmp_day_blocks(pre_sleep_blocks{day}))],'Type',4);
            pre_lfp = data.streams.LFPs;
            pre_video_start_stop_lfp = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
            pre_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];
            pre_lfp = pre_lfp.data;

            pre_lfp_norm = [];
            for n = 1:16
                tmp_lfp = pre_lfp(n,:);
                tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                b = find(abs(tmp_lfp)>6);
                tmp_lfp = pre_lfp(n,:);
                tmp_lfp(b)=NaN;
                tmp_mean=nanmean(tmp_lfp);
                tmp_sd=nanstd(tmp_lfp);
                tmp_lfp(b)=tmp_mean;
                pre_lfp_norm =  [pre_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
            end
            pre_lfp_norm(1:16,:) = pre_lfp_norm(1:16,:)-repmat(mean(pre_lfp_norm(1:16,:),1),[size(pre_lfp_norm(1:16,:),1) 1]);

            load(['D:\Lemke_etal_2021\data\T107\sleep_activity\day_' num2str(day) '_block_2_sleep_activity.mat'])
            post_beh_state = beh_state;
            data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T107\neural_data\Block-' num2str(tmp_day_blocks(post_sleep_blocks{day}))],'Type',4);
            post_lfp = data.streams.LFPs;
            post_video_start_stop_lfp = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
            post_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];
            post_lfp = post_lfp.data;

            post_lfp_norm = [];
            for n = 1:16
                tmp_lfp = post_lfp(n,:);
                tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                b = find(abs(tmp_lfp)>6);
                tmp_lfp = post_lfp(n,:);
                tmp_lfp(b)=NaN;
                tmp_mean=nanmean(tmp_lfp);
                tmp_sd=nanstd(tmp_lfp);
                tmp_lfp(b)=tmp_mean;
                post_lfp_norm =  [post_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
            end
            post_lfp_norm(1:16,:) = post_lfp_norm(1:16,:)-repmat(mean(post_lfp_norm(1:16,:),1),[size(post_lfp_norm(1:16,:),1) 1]);

            %%% LOAD SPIKING DATA

            plx_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\T107\neural_data\' plx_files{day} '.txt']);
            plx_spiking = table2array(plx_spiking);

            tmp_ts = cell(32,length(tmp_day_blocks));
            block_count = 1;
            for block = tmp_day_blocks
                tmp_num = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T107\neural_data\Block-' num2str(block)],'Type',3);
                for chan = 1:32
                    tmp_ts{chan,block_count} = tmp_num.snips.eNeu.ts(tmp_num.snips.eNeu.chan==chan);
                end
                block_count = block_count + 1;
            end

            tmp_length = zeros(32,length(tmp_day_blocks));
            for chan = 1:32
                for block = 1:length(tmp_day_blocks)
                    tmp_length(chan,block) = length(tmp_ts{chan,block});
                end
            end
            tmp_length = [zeros(32,1) tmp_length];

            tmp_sort = cell(32,length(tmp_day_blocks));
            tmp_width = cell(32,length(tmp_day_blocks));
            tmp_peak2valley = cell(32,length(tmp_day_blocks));
            for chan = 1:32
                chan_sorts = plx_spiking(plx_spiking(:,1)==chan,:);
                for block = 1:length(tmp_day_blocks)
                    tmp_sort{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),2);
                    tmp_width{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),6);
                    tmp_peak2valley{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),4);
                end
            end

            pre_dls_spiking = cell(16,max(plx_spiking(:,2)));
            post_dls_spiking = cell(16,max(plx_spiking(:,2)));
            count = 1;
            for dls_chan = 17:32
                for dls_unit = 1:max(plx_spiking(:,2))
                    pre_dls_spiking{count,dls_unit} = tmp_ts{dls_chan,pre_sleep_blocks{day}}(tmp_sort{dls_chan,pre_sleep_blocks{day}}==dls_unit);
                    post_dls_spiking{count,dls_unit} = tmp_ts{dls_chan,post_sleep_blocks{day}}(tmp_sort{dls_chan,post_sleep_blocks{day}}==dls_unit);
                end
                count = count + 1;
            end

            %%% CALCULATE SPIKE COH

            coh_coh = [];

            for m1_lfp = 1:16

                tmp_pre_m1_lfp = pre_lfp_norm(m1_lfp,:);
                tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2));
                tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
                tmp_pre_m1_lfp = eegfilt(tmp_pre_m1_lfp,data.streams.LFPs.fs,4,8);
                tmp_pre_m1_lfp = angle(hilbert(tmp_pre_m1_lfp));

                tmp_post_m1_lfp = post_lfp_norm(m1_lfp,:);
                tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2));
                tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
                tmp_post_m1_lfp = eegfilt(tmp_post_m1_lfp,data.streams.LFPs.fs,4,8);
                tmp_post_m1_lfp = angle(hilbert(tmp_post_m1_lfp));

                tmp_chan_plv = [];

                for dls_chan = 1:size(pre_dls_spiking,1)
                    for dls_unit = 1:size(pre_dls_spiking,2)

                        if isempty(pre_dls_spiking{dls_chan,dls_unit}) || isempty(post_dls_spiking{dls_chan,dls_unit})
                            continue
                        end

                        pre_dls_unit = pre_dls_spiking{dls_chan,dls_unit};%(pre_dls_spiking{dls_chan,dls_unit}>pre_video_start_stop(1) & pre_dls_spiking{dls_chan,dls_unit}<pre_video_start_stop(2))-pre_video_start_stop(1);
                        pre_dls_unit = histc(pre_dls_unit,linspace(pre_video_start_stop(1),pre_video_start_stop(2),length(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2))));
                        pre_dls_unit = pre_dls_unit(pre_beh_state.nrem_interp==1);

                        post_dls_unit = post_dls_spiking{dls_chan,dls_unit};%(post_dls_spiking{dls_chan,dls_unit}>post_video_start_stop(1) & post_dls_spiking{dls_chan,dls_unit}<post_video_start_stop(2))-post_video_start_stop(1);
                        post_dls_unit = histc(post_dls_unit,linspace(post_video_start_stop(1),post_video_start_stop(2),length(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2))));
                        post_dls_unit = post_dls_unit(post_beh_state.nrem_interp==1);

                        [pval, z] = circ_rtest([tmp_pre_m1_lfp(pre_dls_unit==1) tmp_post_m1_lfp(post_dls_unit==1)]);
                        sd = circ_std([tmp_pre_m1_lfp(pre_dls_unit==1) tmp_post_m1_lfp(post_dls_unit==1)]');                        
                        tmp_chan_plv = [tmp_chan_plv; pval z sd pre_coh(m1_lfp)];

                    end
                end

                coh_coh = [coh_coh; mean(tmp_chan_plv)];

            end

            animal_unit_lfp_coh = [animal_unit_lfp_coh; coh_coh];
            animal_unit_lfp_coh_norm = [animal_unit_lfp_coh_norm; coh_coh(:,1) zscore(coh_coh(:,2)) zscore(coh_coh(:,3)) zscore(coh_coh(:,4))]; 

        end

        all_unit_coh_lfp_coh = [all_unit_coh_lfp_coh; animal_unit_lfp_coh];
        all_unit_coh_lfp_coh_day_norm = [all_unit_coh_lfp_coh_day_norm; animal_unit_lfp_coh_norm];
        
        figure;
            subplot(2,2,1); hold on;
                scatter(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,2),30,[0 0 0],'filled')
                [R,P] = corrcoef(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,2));
                [p] = polyfit(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,2),1);
                x=linspace(min(animal_unit_lfp_coh(:,4)),max(animal_unit_lfp_coh(:,4)),10);
                xlim([min(animal_unit_lfp_coh(:,4)) max(animal_unit_lfp_coh(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T107 - lfp coh vs. dls unit m1 lfp z val - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,2); hold on;
                scatter(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,2),30,[0 0 0],'filled');
                [R,P] = corrcoef(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,2));
                [p] = polyfit(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,2),1);
                x=linspace(min(animal_unit_lfp_coh_norm(:,4)),max(animal_unit_lfp_coh_norm(:,4)),10);
                xlim([min(animal_unit_lfp_coh_norm(:,4)) max(animal_unit_lfp_coh_norm(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T107 - norm lfp coh vs. norm dls unit m1 lfp z val - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,3); hold on;
                scatter(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,3),30,[0 0 0],'filled')
                [R,P] = corrcoef(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,3));
                [p] = polyfit(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,3),1);
                x=linspace(min(animal_unit_lfp_coh(:,4)),max(animal_unit_lfp_coh(:,4)),10);
                xlim([min(animal_unit_lfp_coh(:,4)) max(animal_unit_lfp_coh(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T107 - lfp coh vs. dls unit m1 lfp SD - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,4); hold on;
                scatter(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,3),30,[0 0 0],'filled');
                [R,P] = corrcoef(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,3));
                [p] = polyfit(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,3),1);
                x=linspace(min(animal_unit_lfp_coh_norm(:,4)),max(animal_unit_lfp_coh_norm(:,4)),10);
                xlim([min(animal_unit_lfp_coh_norm(:,4)) max(animal_unit_lfp_coh_norm(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T107 - norm lfp coh vs. norm dls unit m1 lfp SD - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
                
%%% T200

        animal_unit_lfp_coh = [];
        animal_unit_lfp_coh_norm = [];

        day_blocks = {[1],[2 7 8],[9 11 12 13 14],[15 16 17 18 19],[20 21 22 23 24],[25 26 27],[28 29 30],[32 33 34 35 36],[37 38 39 40 41],[42 43 44 45 46 47]};
        pre_sleep_blocks = {[1],[1],[2],[1],[1],[1],[1],[1],[1],[1]};
        post_sleep_blocks = {[],[3],[5],[5],[4],[3],[3],[5],[5],[6]};
        plx_files = {'T200_blocks_1-01','T200_blocks_2_7_8-01','T200_blocks_9_11_12_13_14-01','T200_blocks_15_16_17_18_19-01','T200_blocks_20_21_22_23_24-01','T200_blocks_25_26_27-01','T200_blocks_28_29_30-01','T200_blocks_32_33_34_35_36-01','T200_blocks_37_38_39_40_41-01','T200_blocks_42_43_44_45_46_47-01'};

        load(['D:\Lemke_etal_2021\data\T200\sleep_activity\LFP_coherence.mat'])

        for day = 2:10

            tmp_day_blocks = day_blocks{day};

            %%% PRE AND POST SLEEP NREM COH

            pre_coh = [];
            post_coh = [];
            for pairs=1:256
                pre_coh = [pre_coh mean(pre_LFP_coh{pairs}(day-1,48:113))];
                post_coh = [post_coh mean(post_LFP_coh{pairs}(day-1,48:113))];
            end
            pre_coh = mean(reshape(pre_coh,16,16));
            post_coh = mean(reshape(post_coh,16,16));

            %%% LOAD BEH STATE, LFP, AND START/STOP

            load(['D:\Lemke_etal_2021\data\T200\sleep_activity\day_' num2str(day) '_block_1_sleep_activity.mat'])
            pre_beh_state = beh_state;
            data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\Block-' num2str(tmp_day_blocks(pre_sleep_blocks{day}))],'Type',4);
            pre_lfp = data.streams.LFPs;
            pre_video_start_stop_lfp = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
            pre_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];
            pre_lfp = pre_lfp.data;

            pre_lfp_norm = [];
            for n = 17:32
                tmp_lfp = pre_lfp(n,:);
                tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                b = find(abs(tmp_lfp)>6);
                tmp_lfp = pre_lfp(n,:);
                tmp_lfp(b)=NaN;
                tmp_mean=nanmean(tmp_lfp);
                tmp_sd=nanstd(tmp_lfp);
                tmp_lfp(b)=tmp_mean;
                pre_lfp_norm =  [pre_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
            end
            pre_lfp_norm(1:16,:) = pre_lfp_norm(1:16,:)-repmat(mean(pre_lfp_norm(1:16,:),1),[size(pre_lfp_norm(1:16,:),1) 1]);

            load(['D:\Lemke_etal_2021\data\T200\sleep_activity\day_' num2str(day) '_block_2_sleep_activity.mat'])
            post_beh_state = beh_state;
            data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\Block-' num2str(tmp_day_blocks(post_sleep_blocks{day}))],'Type',4);
            post_lfp = data.streams.LFPs;
            post_video_start_stop_lfp = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
            post_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];
            post_lfp = post_lfp.data;

            post_lfp_norm = [];
            for n = 17:32
                tmp_lfp = post_lfp(n,:);
                tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                b = find(abs(tmp_lfp)>6);
                tmp_lfp = post_lfp(n,:);
                tmp_lfp(b)=NaN;
                tmp_mean=nanmean(tmp_lfp);
                tmp_sd=nanstd(tmp_lfp);
                tmp_lfp(b)=tmp_mean;
                post_lfp_norm =  [post_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
            end
            post_lfp_norm(1:16,:) = post_lfp_norm(1:16,:)-repmat(mean(post_lfp_norm(1:16,:),1),[size(post_lfp_norm(1:16,:),1) 1]);

            %%% LOAD SPIKING DATA

            plx_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\' plx_files{day} '.txt']);
            plx_spiking = table2array(plx_spiking);

            tmp_ts = cell(32,length(tmp_day_blocks));
            block_count = 1;
            for block = tmp_day_blocks
                tmp_num = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\Block-' num2str(block)],'Type',3);
                for chan = 1:32
                    tmp_ts{chan,block_count} = tmp_num.snips.eNeu.ts(tmp_num.snips.eNeu.chan==chan);
                end
                block_count = block_count + 1;
            end

            tmp_length = zeros(32,length(tmp_day_blocks));
            for chan = 1:32
                for block = 1:length(tmp_day_blocks)
                    tmp_length(chan,block) = length(tmp_ts{chan,block});
                end
            end
            tmp_length = [zeros(32,1) tmp_length];

            tmp_sort = cell(32,length(tmp_day_blocks));
            tmp_width = cell(32,length(tmp_day_blocks));
            tmp_peak2valley = cell(32,length(tmp_day_blocks));
            for chan = 1:32
                chan_sorts = plx_spiking(plx_spiking(:,1)==chan,:);
                for block = 1:length(tmp_day_blocks)
                    tmp_sort{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),2);
                    tmp_width{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),6);
                    tmp_peak2valley{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),4);
                end
            end

            pre_dls_spiking = cell(16,max(plx_spiking(:,2)));
            post_dls_spiking = cell(16,max(plx_spiking(:,2)));
            count = 1;
            for dls_chan = 1:16
                for dls_unit = 1:max(plx_spiking(:,2))
                    pre_dls_spiking{count,dls_unit} = tmp_ts{dls_chan,pre_sleep_blocks{day}}(tmp_sort{dls_chan,pre_sleep_blocks{day}}==dls_unit);
                    post_dls_spiking{count,dls_unit} = tmp_ts{dls_chan,post_sleep_blocks{day}}(tmp_sort{dls_chan,post_sleep_blocks{day}}==dls_unit);
                end
                count = count + 1;
            end

            %%% CALCULATE SPIKE COH

            coh_coh = [];

            for m1_lfp = 1:16

                tmp_pre_m1_lfp = pre_lfp_norm(m1_lfp,:);
                if day ~= 6
                    tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2));
                end
                tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
                tmp_pre_m1_lfp = eegfilt(tmp_pre_m1_lfp,data.streams.LFPs.fs,4,8);
                tmp_pre_m1_lfp = angle(hilbert(tmp_pre_m1_lfp));

                tmp_post_m1_lfp = post_lfp_norm(m1_lfp,:);
                tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2));
                tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
                tmp_post_m1_lfp = eegfilt(tmp_post_m1_lfp,data.streams.LFPs.fs,4,8);
                tmp_post_m1_lfp = angle(hilbert(tmp_post_m1_lfp));

                tmp_chan_plv = [];

                for dls_chan = 1:size(pre_dls_spiking,1)
                    for dls_unit = 1:size(pre_dls_spiking,2)

                        if isempty(pre_dls_spiking{dls_chan,dls_unit}) || isempty(post_dls_spiking{dls_chan,dls_unit})
                            continue
                        end

                        pre_dls_unit = pre_dls_spiking{dls_chan,dls_unit};%(pre_dls_spiking{dls_chan,dls_unit}>pre_video_start_stop(1) & pre_dls_spiking{dls_chan,dls_unit}<pre_video_start_stop(2))-pre_video_start_stop(1);
                        if day ~= 6
                            pre_dls_unit = histc(pre_dls_unit,linspace(pre_video_start_stop(1),pre_video_start_stop(2),length(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2))));
                        else
                            pre_dls_unit = histc(pre_dls_unit,linspace(1,length(pre_beh_state.nrem_interp),length(pre_beh_state.nrem_interp)));                        
                        end
                        pre_dls_unit = pre_dls_unit(pre_beh_state.nrem_interp==1);

                        post_dls_unit = post_dls_spiking{dls_chan,dls_unit};%(post_dls_spiking{dls_chan,dls_unit}>post_video_start_stop(1) & post_dls_spiking{dls_chan,dls_unit}<post_video_start_stop(2))-post_video_start_stop(1);
                        post_dls_unit = histc(post_dls_unit,linspace(post_video_start_stop(1),post_video_start_stop(2),length(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2))));
                        post_dls_unit = post_dls_unit(post_beh_state.nrem_interp==1);

                        [pval, z] = circ_rtest([tmp_pre_m1_lfp(pre_dls_unit==1) tmp_post_m1_lfp(post_dls_unit==1)]);
                        sd = circ_std([tmp_pre_m1_lfp(pre_dls_unit==1) tmp_post_m1_lfp(post_dls_unit==1)]');                        
                        tmp_chan_plv = [tmp_chan_plv; pval z sd pre_coh(m1_lfp)];

                    end
                end

                coh_coh = [coh_coh; mean(tmp_chan_plv)];

            end

            animal_unit_lfp_coh = [animal_unit_lfp_coh; coh_coh];
            animal_unit_lfp_coh_norm = [animal_unit_lfp_coh_norm; coh_coh(:,1) zscore(coh_coh(:,2)) zscore(coh_coh(:,3)) zscore(coh_coh(:,4))]; 

        end

        all_unit_coh_lfp_coh = [all_unit_coh_lfp_coh; animal_unit_lfp_coh];
        all_unit_coh_lfp_coh_day_norm = [all_unit_coh_lfp_coh_day_norm; animal_unit_lfp_coh_norm];
             
        figure;
            subplot(2,2,1); hold on;
                scatter(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,2),30,[0 0 0],'filled')
                [R,P] = corrcoef(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,2));
                [p] = polyfit(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,2),1);
                x=linspace(min(animal_unit_lfp_coh(:,4)),max(animal_unit_lfp_coh(:,4)),10);
                xlim([min(animal_unit_lfp_coh(:,4)) max(animal_unit_lfp_coh(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T200 - lfp coh vs. dls unit m1 lfp z val - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,2); hold on;
                scatter(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,2),30,[0 0 0],'filled');
                [R,P] = corrcoef(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,2));
                [p] = polyfit(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,2),1);
                x=linspace(min(animal_unit_lfp_coh_norm(:,4)),max(animal_unit_lfp_coh_norm(:,4)),10);
                xlim([min(animal_unit_lfp_coh_norm(:,4)) max(animal_unit_lfp_coh_norm(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T200 - norm lfp coh vs. norm dls unit m1 lfp z val - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,3); hold on;
                scatter(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,3),30,[0 0 0],'filled')
                [R,P] = corrcoef(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,3));
                [p] = polyfit(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,3),1);
                x=linspace(min(animal_unit_lfp_coh(:,4)),max(animal_unit_lfp_coh(:,4)),10);
                xlim([min(animal_unit_lfp_coh(:,4)) max(animal_unit_lfp_coh(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T200 - lfp coh vs. dls unit m1 lfp SD - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,4); hold on;
                scatter(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,3),30,[0 0 0],'filled');
                [R,P] = corrcoef(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,3));
                [p] = polyfit(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,3),1);
                x=linspace(min(animal_unit_lfp_coh_norm(:,4)),max(animal_unit_lfp_coh_norm(:,4)),10);
                xlim([min(animal_unit_lfp_coh_norm(:,4)) max(animal_unit_lfp_coh_norm(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T200 - norm lfp coh vs. norm dls unit m1 lfp SD - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
                
%%% T201

        animal_unit_lfp_coh = [];
        animal_unit_lfp_coh_norm = [];

        day_blocks = {[1],[2 3],[4 5 6 7],[8 9 10],[11 12 13],[14 15 16 17],[18 19 20],[21 22 23 24 25 26],[27 28 29], ...
            [30 31 32],[33 34 35 36 37],[38 39 40],[41 42 43 44],[45 46 47 48],[49 50 51 52],[53 54 55 56],[57 58 59 60 61]};
        pre_sleep_blocks = {[1],[1],[2],[1],[1],[2],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1]};
        post_sleep_blocks = {[],[2],[4],[3],[3],[4],[3],[6],[3],[3],[5],[3],[4],[4],[4],[4],[5]};
        plx_files = {'T201_blocks_1-01','T201_blocks_2_3-01','T201_blocks_4_5_6_7-01','T201_blocks_8_9_10-01','T201_blocks_11_12_13-01', ...
            'T201_blocks_14_15_16_17-01','T201_blocks_18_19_20-01','T201_blocks_21_22_23_24_25_26-01','T201_blocks_27_28_29-01', ...
            'T201_blocks_30_31_32-01','T201_blocks_33_34_35_36_37-01','T201_blocks_38_39_40-01','T201_blocks_41_42_43_44-01','T201_blocks_45_46_47_48-01', ...
            'T201_blocks_49_50_51_52-01','T201_blocks_53_54_55_56-01','T201_blocks_57_58_59_60_61-01'};

        load(['D:\Lemke_etal_2021\data\T201\sleep_activity\LFP_coherence.mat'])

        for day = 2:16

            tmp_day_blocks = day_blocks{day};

            %%% PRE AND POST SLEEP NREM COH

            pre_coh = [];
            post_coh = [];
            for pairs=1:512
                pre_coh = [pre_coh mean(pre_LFP_coh{pairs}(day-1,48:113))];
                post_coh = [post_coh mean(post_LFP_coh{pairs}(day-1,48:113))];
            end
            pre_coh = mean(reshape(pre_coh,16,32));
            post_coh = mean(reshape(post_coh,16,32));

            %%% LOAD BEH STATE, LFP, AND START/STOP

            load(['D:\Lemke_etal_2021\data\T201\sleep_activity\day_' num2str(day) '_block_1_sleep_activity.mat'])
            pre_beh_state = beh_state;
            if day == 13
                data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(tmp_day_blocks(pre_sleep_blocks{day}))],'Type',4,'T2',7200);
            else
                data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(tmp_day_blocks(pre_sleep_blocks{day}))],'Type',4);
            end
            pre_lfp = data.streams.LFPs;
            pre_video_start_stop_lfp = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
            pre_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];
            pre_lfp = pre_lfp.data;

            pre_lfp_norm = [];
            for n = 1:2:63
                tmp_lfp = pre_lfp(n,:);
                tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                b = find(abs(tmp_lfp)>6);
                tmp_lfp = pre_lfp(n,:);
                tmp_lfp(b)=NaN;
                tmp_mean=nanmean(tmp_lfp);
                tmp_sd=nanstd(tmp_lfp);
                tmp_lfp(b)=tmp_mean;
                pre_lfp_norm =  [pre_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
            end
            pre_lfp_norm(1:32,:) = pre_lfp_norm(1:32,:)-repmat(mean(pre_lfp_norm(1:32,:),1),[size(pre_lfp_norm(1:32,:),1) 1]);

            load(['D:\Lemke_etal_2021\data\T201\sleep_activity\day_' num2str(day) '_block_2_sleep_activity.mat'])
            post_beh_state = beh_state;
            data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(tmp_day_blocks(post_sleep_blocks{day}))],'Type',4);
            post_lfp = data.streams.LFPs;
            post_video_start_stop_lfp = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
            post_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];
            post_lfp = post_lfp.data;

            post_lfp_norm = [];
            for n = 1:2:63
                tmp_lfp = post_lfp(n,:);
                tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                b = find(abs(tmp_lfp)>6);
                tmp_lfp = post_lfp(n,:);
                tmp_lfp(b)=NaN;
                tmp_mean=nanmean(tmp_lfp);
                tmp_sd=nanstd(tmp_lfp);
                tmp_lfp(b)=tmp_mean;
                post_lfp_norm =  [post_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
            end
            post_lfp_norm(1:32,:) = post_lfp_norm(1:32,:)-repmat(mean(post_lfp_norm(1:32,:),1),[size(post_lfp_norm(1:32,:),1) 1]);

            %%% LOAD SPIKING DATA

            plx_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\' plx_files{day} '.txt']);
            plx_spiking = table2array(plx_spiking);

            tmp_ts = cell(64,length(tmp_day_blocks));
            block_count = 1;
            for block = tmp_day_blocks
                if block == 41
                    tmp_num = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(block)],'Type',3,'T2',7200);
                else
                    tmp_num = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(block)],'Type',3);
                end
                for chan = 1:64
                    tmp_ts{chan,block_count} = tmp_num.snips.eNeu.ts(tmp_num.snips.eNeu.chan==chan);
                end
                block_count = block_count + 1;
            end

            tmp_length = zeros(64,length(tmp_day_blocks));
            for chan = 1:64
                for block = 1:length(tmp_day_blocks)
                    tmp_length(chan,block) = length(tmp_ts{chan,block});
                end
            end
            tmp_length = [zeros(64,1) tmp_length];

            tmp_sort = cell(64,length(tmp_day_blocks));
            for chan = 1:64
                chan_sorts = plx_spiking(plx_spiking(:,1)==chan,:);
                for block = 1:length(tmp_day_blocks)
                    tmp_sort{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),2);
                end
            end

            pre_dls_spiking = cell(16,max(plx_spiking(:,2)));
            post_dls_spiking = cell(16,max(plx_spiking(:,2)));
            count = 1;
            for dls_chan = 34:2:64
                for dls_unit = 1:max(plx_spiking(:,2))
                    pre_dls_spiking{count,dls_unit} = tmp_ts{dls_chan,pre_sleep_blocks{day}}(tmp_sort{dls_chan,pre_sleep_blocks{day}}==dls_unit);
                    post_dls_spiking{count,dls_unit} = tmp_ts{dls_chan,post_sleep_blocks{day}}(tmp_sort{dls_chan,post_sleep_blocks{day}}==dls_unit);
                end
                count = count + 1;
            end

            dls_exist = [];
            count = 1;
            for dls_chan = 34:2:64
                for dls_unit = 1:max(plx_spiking(:,2))
                    dls_exist = [dls_exist isempty(pre_dls_spiking{count,dls_unit})];
                end
                count = count + 1;
            end
            if min(dls_exist)==1
                continue
            end

            %%% CALCULATE SPIKE COH

            coh_coh = [];

            for m1_lfp = 1:32

                tmp_pre_m1_lfp = pre_lfp_norm(m1_lfp,:);
                tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2));
                tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
                tmp_pre_m1_lfp = eegfilt(tmp_pre_m1_lfp,data.streams.LFPs.fs,4,8);
                tmp_pre_m1_lfp = angle(hilbert(tmp_pre_m1_lfp));

                tmp_post_m1_lfp = post_lfp_norm(m1_lfp,:);
                tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2));
                tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
                tmp_post_m1_lfp = eegfilt(tmp_post_m1_lfp,data.streams.LFPs.fs,4,8);
                tmp_post_m1_lfp = angle(hilbert(tmp_post_m1_lfp));

                tmp_chan_plv = [];

                for dls_chan = 1:size(pre_dls_spiking,1)
                    for dls_unit = 1:size(pre_dls_spiking,2)

                        if isempty(pre_dls_spiking{dls_chan,dls_unit}) || isempty(post_dls_spiking{dls_chan,dls_unit})
                            continue
                        end

                        pre_dls_unit = pre_dls_spiking{dls_chan,dls_unit};%(pre_dls_spiking{dls_chan,dls_unit}>pre_video_start_stop(1) & pre_dls_spiking{dls_chan,dls_unit}<pre_video_start_stop(2))-pre_video_start_stop(1);
                        pre_dls_unit = histc(pre_dls_unit,linspace(pre_video_start_stop(1),pre_video_start_stop(2),length(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2))));
                        pre_dls_unit = pre_dls_unit(pre_beh_state.nrem_interp==1);

                        post_dls_unit = post_dls_spiking{dls_chan,dls_unit};%(post_dls_spiking{dls_chan,dls_unit}>post_video_start_stop(1) & post_dls_spiking{dls_chan,dls_unit}<post_video_start_stop(2))-post_video_start_stop(1);
                        post_dls_unit = histc(post_dls_unit,linspace(post_video_start_stop(1),post_video_start_stop(2),length(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2))));
                        post_dls_unit = post_dls_unit(post_beh_state.nrem_interp==1);

                        [pval, z] = circ_rtest([tmp_pre_m1_lfp(pre_dls_unit==1) tmp_post_m1_lfp(post_dls_unit==1)]);
                        sd = circ_std([tmp_pre_m1_lfp(pre_dls_unit==1) tmp_post_m1_lfp(post_dls_unit==1)]');                        
                        tmp_chan_plv = [tmp_chan_plv; pval z sd pre_coh(m1_lfp)];

                    end
                end

                coh_coh = [coh_coh; mean(tmp_chan_plv)];

            end

            animal_unit_lfp_coh = [animal_unit_lfp_coh; coh_coh];
            animal_unit_lfp_coh_norm = [animal_unit_lfp_coh_norm; coh_coh(:,1) zscore(coh_coh(:,2)) zscore(coh_coh(:,3)) zscore(coh_coh(:,4))]; 

        end

        all_unit_coh_lfp_coh = [all_unit_coh_lfp_coh; animal_unit_lfp_coh];
        all_unit_coh_lfp_coh_day_norm = [all_unit_coh_lfp_coh_day_norm; animal_unit_lfp_coh_norm];
        
        figure;
            subplot(2,2,1); hold on;
                scatter(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,2),30,[0 0 0],'filled')
                [R,P] = corrcoef(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,2));
                [p] = polyfit(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,2),1);
                x=linspace(min(animal_unit_lfp_coh(:,4)),max(animal_unit_lfp_coh(:,4)),10);
                xlim([min(animal_unit_lfp_coh(:,4)) max(animal_unit_lfp_coh(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T201 - lfp coh vs. dls unit m1 lfp z val - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,2); hold on;
                scatter(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,2),30,[0 0 0],'filled');
                [R,P] = corrcoef(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,2));
                [p] = polyfit(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,2),1);
                x=linspace(min(animal_unit_lfp_coh_norm(:,4)),max(animal_unit_lfp_coh_norm(:,4)),10);
                xlim([min(animal_unit_lfp_coh_norm(:,4)) max(animal_unit_lfp_coh_norm(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T201 - norm lfp coh vs. norm dls unit m1 lfp z val - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,3); hold on;
                scatter(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,3),30,[0 0 0],'filled')
                [R,P] = corrcoef(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,3));
                [p] = polyfit(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,3),1);
                x=linspace(min(animal_unit_lfp_coh(:,4)),max(animal_unit_lfp_coh(:,4)),10);
                xlim([min(animal_unit_lfp_coh(:,4)) max(animal_unit_lfp_coh(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T201 - lfp coh vs. dls unit m1 lfp SD - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,4); hold on;
                scatter(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,3),30,[0 0 0],'filled');
                [R,P] = corrcoef(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,3));
                [p] = polyfit(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,3),1);
                x=linspace(min(animal_unit_lfp_coh_norm(:,4)),max(animal_unit_lfp_coh_norm(:,4)),10);
                xlim([min(animal_unit_lfp_coh_norm(:,4)) max(animal_unit_lfp_coh_norm(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T201 - norm lfp coh vs. norm dls unit m1 lfp SD - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
                           
%%% SLDD2

        animal_unit_lfp_coh = [];
        animal_unit_lfp_coh_norm = [];

        day_blocks = {{'SLDD2-191113-143329'}, ...
            {'SLDD2-191114-121834'}, ...
            {'SLDD2-191115-154323'}, ...
            {'SLDD2-191116-145119','SLDD2-191116-170045','SLDD2-191116-193647'}, ...
            {'SLDD2-191117-150115','SLDD2-191117-170253','SLDD2-191117-181931'}, ...
            {'SLDD2-191118-154100','SLDD2-191118-180452','SLDD2-191118-192527','SLDD2-191118-205741'}, ...
            {'SLDD2-191119-160043','SLDD2-191119-180437','SLDD2-191119-183709','SLDD2-191119-190631','SLDD2-191119-193444'}, ...
            {'SLDD2-191120-140404','SLDD2-191120-160657','SLDD2-191120-170507'}, ...
            {'SLDD2-191121-170233','SLDD2-191121-180559','SLDD2-191121-182818','SLDD2-191121-192200'}, ...
            {'SLDD2-191122-152109','SLDD2-191122-172516','SLDD2-191122-184340','SLDD2-191122-204428'}, ...
            {'SLDD2-191123-173250','SLDD2-191123-193545','SLDD2-191123-205347'}, ...
            {'SLDD2-191124-152841','SLDD2-191124-173146','SLDD2-191124-183235'}, ...
            {'SLDD2-191125-141910','SLDD2-191125-162007','SLDD2-191125-171539'}, ...
            {'SLDD2-191126-173517','SLDD2-191126-191048','SLDD2-191126-195858'}};

        pre_sleep_blocks = {[],[],[],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1]};
        post_sleep_blocks = {[],[],[],[3],[3],[4],[5],[3],[4],[3],[3],[3],[3],[3]};

        pre_sleep_plx_blocks = {[],[],[],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1]};
        post_sleep_plx_blocks = {[],[],[],[3],[3],[4],[5],[3],[4],[3],[3],[3],[3],[3]};

        plx_files = {'sldd2-191113-143329-01', ...
            'sldd2-191114-121834-01', ...
            'sldd2-191115-154323-01', ...
            'sldd2-191116-145119_mrg-01', ...
            'sldd2-191117-150115_mrg-01', ...
            'sldd2-191118-154100_mrg-01', ...
            'sldd2-191119-160043_mrg-01', ...
            'sldd2-191120-140404_mrg-01', ...
            'sldd2-191121-170233_mrg-01', ...
            'sldd2-191122-152109_mrg-01', ...
            'sldd2-191123-173250_mrg-01', ...
            'sldd2-191124-152841_mrg-01', ...
            'sldd2-191125-141910_mrg-01', ...
            'sldd2-191126-173517_mrg-01'};

        load(['D:\Lemke_etal_2021\data\SLDD2\sleep_activity\LFP_coherence.mat'])

        for day = 4:13

            tmp_sleep_blocks = day_blocks{day};

            %%% PRE AND POST SLEEP NREM COH

            pre_coh = [];
            post_coh = [];
            for pairs=1:1024
                pre_coh = [pre_coh mean(pre_LFP_coh{pairs}(day-3,40:95))];
                post_coh = [post_coh mean(post_LFP_coh{pairs}(day-3,40:95))];
            end
            pre_coh = mean(reshape(pre_coh,32,32));
            post_coh = mean(reshape(post_coh,32,32));

            %%% LOAD BEH STATE AND START/STOP

            load(['D:\Lemke_etal_2021\data\SLDD2\sleep_activity\day_' num2str(day) '_block_1_sleep_activity.mat'])
            pre_beh_state = beh_state;
            wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\TDT_tanks\' day_blocks{day}{pre_sleep_blocks{day}}],'TYPE',[2]);
            pre_lfp = [];
            for chan = 1:64
                data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\RS4_data\' day_blocks{day}{pre_sleep_blocks{day}}],'CHANNEL',chan);
                pre_lfp = [pre_lfp; downsample(data.RSn1.data,20)];
            end
            pre_video_start_stop_lfp = [min(wave.epocs.PtC0.onset)*(data.RSn1.fs/20) max(wave.epocs.PtC0.onset)*(data.RSn1.fs/20)];
            pre_video_start_stop = [min(wave.epocs.PtC0.onset) max(wave.epocs.PtC0.onset)];

            pre_lfp_norm = [];
            for n = 1:size(pre_lfp,1)
                tmp_lfp = double(pre_lfp(n,:));
                tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                b = find(abs(tmp_lfp)>6);
                tmp_lfp = double(pre_lfp(n,:));
                tmp_lfp(b)=NaN;
                tmp_mean=nanmean(tmp_lfp);
                tmp_sd=nanstd(tmp_lfp);
                tmp_lfp(b)=tmp_mean;
                pre_lfp_norm =  [pre_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
            end
            pre_lfp_norm([8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33],:) = pre_lfp_norm([8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33],:)-repmat(mean(pre_lfp_norm([8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33],:),1),[size(pre_lfp_norm([8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33],:),1) 1]);
            pre_lfp_norm([2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63],:) = pre_lfp_norm([2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63],:)-repmat(mean(pre_lfp_norm([2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63],:),1),[size(pre_lfp_norm([2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63],:),1) 1]);

            load(['D:\Lemke_etal_2021\data\SLDD2\sleep_activity\day_' num2str(day) '_block_2_sleep_activity.mat'])
            post_beh_state = beh_state;
            wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\TDT_tanks\' day_blocks{day}{post_sleep_blocks{day}}],'TYPE',[2]);
            post_lfp = [];
            for chan = 1:64
                data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\RS4_data\' day_blocks{day}{post_sleep_blocks{day}}],'CHANNEL',chan);
                post_lfp = [post_lfp; downsample(data.RSn1.data,20)];
            end
            post_video_start_stop_lfp = [min(wave.epocs.PtC0.onset)*(data.RSn1.fs/20) max(wave.epocs.PtC0.onset)*(data.RSn1.fs/20)];
            post_video_start_stop = [min(wave.epocs.PtC0.onset) max(wave.epocs.PtC0.onset)];

            post_lfp_norm = [];
            for n = 1:size(post_lfp,1)
                tmp_lfp = double(post_lfp(n,:));
                tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                b = find(abs(tmp_lfp)>6);
                tmp_lfp = double(post_lfp(n,:));
                tmp_lfp(b)=NaN;
                tmp_mean=nanmean(tmp_lfp);
                tmp_sd=nanstd(tmp_lfp);
                tmp_lfp(b)=tmp_mean;
                post_lfp_norm =  [post_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
            end
            post_lfp_norm([8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33],:) = post_lfp_norm([8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33],:)-repmat(mean(post_lfp_norm([8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33],:),1),[size(post_lfp_norm([8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33],:),1) 1]);
            post_lfp_norm([2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63],:) = post_lfp_norm([2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63],:)-repmat(mean(post_lfp_norm([2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63],:),1),[size(post_lfp_norm([2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63],:),1) 1]);

            %%% LOAD SPIKING DATA

            plx_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\PLX_sorts\day_' num2str(day) '\' plx_files{day} '.txt']);
            plx_spiking = table2array(plx_spiking);

            tmp_blocks = ls(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\PLX_sorts\day_' num2str(day) '\*.txt']);
            tmp_blocks_cell = cell(size(tmp_blocks,1),1);
            for block = 1:size(tmp_blocks,1)
                tmp_blocks_cell{block} = deblank(tmp_blocks(block,:));
            end
            tmp_ind = ~contains(tmp_blocks_cell,'mrg');
            tmp_blocks_cell = tmp_blocks_cell(tmp_ind);

            tmp_ts = cell(64,length(tmp_blocks_cell));
            for block = 1:length(tmp_blocks_cell)
                tmp_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\PLX_sorts\day_' num2str(day) '\' tmp_blocks_cell{block}]);
                tmp_spiking = table2array(tmp_spiking);
                for chan = 1:64
                    tmp_ts{chan,block} = tmp_spiking(tmp_spiking(:,1)==chan,2);
                end
            end

            tmp_length = zeros(64,length(tmp_blocks_cell));
            for chan = 1:64
                for block = 1:length(tmp_blocks_cell)
                    tmp_length(chan,block) = length(tmp_ts{chan,block});
                end
            end
            tmp_length = [zeros(64,1) tmp_length];

            tmp_sort = cell(64,length(tmp_blocks_cell));
            for chan = 1:64
                chan_sorts = plx_spiking(plx_spiking(:,1)==chan,:);
                for block = 1:length(tmp_blocks_cell)
                    tmp_sort{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),2);
                end
            end

            pre_dls_spiking = cell(32,max(plx_spiking(:,2)));
            post_dls_spiking = cell(32,max(plx_spiking(:,2)));
            count = 1;
            for dls_chan = [2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63]
                for dls_unit = 1:max(plx_spiking(:,2))
                    pre_dls_spiking{count,dls_unit} = tmp_ts{dls_chan,pre_sleep_plx_blocks{day}}(tmp_sort{dls_chan,pre_sleep_plx_blocks{day}}==dls_unit);
                    post_dls_spiking{count,dls_unit} = tmp_ts{dls_chan,post_sleep_plx_blocks{day}}(tmp_sort{dls_chan,post_sleep_plx_blocks{day}}==dls_unit);
                end
                count = count + 1;
            end

            %%% CALCULATE SPIKE COH

            coh_coh = [];

            m1_chan_count = 1;
            for m1_lfp = [8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33]

                tmp_pre_m1_lfp = pre_lfp_norm(m1_lfp,:);
                tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2));
                tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
                tmp_pre_m1_lfp = eegfilt(tmp_pre_m1_lfp,(data.RSn1.fs/20),4,8);
                tmp_pre_m1_lfp = angle(hilbert(tmp_pre_m1_lfp));

                tmp_post_m1_lfp = post_lfp_norm(m1_lfp,:);
                tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2));
                tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
                tmp_post_m1_lfp = eegfilt(tmp_post_m1_lfp,(data.RSn1.fs/20),4,8);
                tmp_post_m1_lfp = angle(hilbert(tmp_post_m1_lfp));

                tmp_chan_plv = [];

                for dls_chan = 1:size(pre_dls_spiking,1)
                    for dls_unit = 1:size(pre_dls_spiking,2)

                        if isempty(pre_dls_spiking{dls_chan,dls_unit}) || isempty(post_dls_spiking{dls_chan,dls_unit})
                            continue
                        end

                        pre_dls_unit = pre_dls_spiking{dls_chan,dls_unit};%(pre_dls_spiking{dls_chan,dls_unit}>pre_video_start_stop(1) & pre_dls_spiking{dls_chan,dls_unit}<pre_video_start_stop(2))-pre_video_start_stop(1);
                        pre_dls_unit = histc(pre_dls_unit,linspace(pre_video_start_stop(1),pre_video_start_stop(2),length(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2))));
                        pre_dls_unit = pre_dls_unit(pre_beh_state.nrem_interp==1);

                        post_dls_unit = post_dls_spiking{dls_chan,dls_unit};%(post_dls_spiking{dls_chan,dls_unit}>post_video_start_stop(1) & post_dls_spiking{dls_chan,dls_unit}<post_video_start_stop(2))-post_video_start_stop(1);
                        post_dls_unit = histc(post_dls_unit,linspace(post_video_start_stop(1),post_video_start_stop(2),length(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2))));
                        post_dls_unit = post_dls_unit(post_beh_state.nrem_interp==1);

                        [pval, z] = circ_rtest([tmp_pre_m1_lfp(pre_dls_unit==1) tmp_post_m1_lfp(post_dls_unit==1)]);
                        sd = circ_std([tmp_pre_m1_lfp(pre_dls_unit==1) tmp_post_m1_lfp(post_dls_unit==1)]');                        
                        tmp_chan_plv = [tmp_chan_plv; pval z sd pre_coh(m1_chan_count)];

                    end
                end

                coh_coh = [coh_coh; mean(tmp_chan_plv)];
                m1_chan_count = m1_chan_count + 1;

            end

            animal_unit_lfp_coh = [animal_unit_lfp_coh; coh_coh];
            animal_unit_lfp_coh_norm = [animal_unit_lfp_coh_norm; coh_coh(:,1) zscore(coh_coh(:,2)) zscore(coh_coh(:,3)) zscore(coh_coh(:,4))]; 

        end

        all_unit_coh_lfp_coh = [all_unit_coh_lfp_coh; animal_unit_lfp_coh];
        all_unit_coh_lfp_coh_day_norm = [all_unit_coh_lfp_coh_day_norm; animal_unit_lfp_coh_norm];
        
        figure;
            subplot(2,2,1); hold on;
                scatter(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,2),30,[0 0 0],'filled')
                [R,P] = corrcoef(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,2));
                [p] = polyfit(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,2),1);
                x=linspace(min(animal_unit_lfp_coh(:,4)),max(animal_unit_lfp_coh(:,4)),10);
                xlim([min(animal_unit_lfp_coh(:,4)) max(animal_unit_lfp_coh(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['SLDD2 - lfp coh vs. dls unit m1 lfp z val - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,2); hold on;
                scatter(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,2),30,[0 0 0],'filled');
                [R,P] = corrcoef(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,2));
                [p] = polyfit(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,2),1);
                x=linspace(min(animal_unit_lfp_coh_norm(:,4)),max(animal_unit_lfp_coh_norm(:,4)),10);
                xlim([min(animal_unit_lfp_coh_norm(:,4)) max(animal_unit_lfp_coh_norm(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['SLDD2 - norm lfp coh vs. norm dls unit m1 lfp z val - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,3); hold on;
                scatter(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,3),30,[0 0 0],'filled')
                [R,P] = corrcoef(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,3));
                [p] = polyfit(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,3),1);
                x=linspace(min(animal_unit_lfp_coh(:,4)),max(animal_unit_lfp_coh(:,4)),10);
                xlim([min(animal_unit_lfp_coh(:,4)) max(animal_unit_lfp_coh(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['SLDD2 - lfp coh vs. dls unit m1 lfp SD - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,4); hold on;
                scatter(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,3),30,[0 0 0],'filled');
                [R,P] = corrcoef(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,3));
                [p] = polyfit(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,3),1);
                x=linspace(min(animal_unit_lfp_coh_norm(:,4)),max(animal_unit_lfp_coh_norm(:,4)),10);
                xlim([min(animal_unit_lfp_coh_norm(:,4)) max(animal_unit_lfp_coh_norm(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['SLDD2 - norm lfp coh vs. norm dls unit m1 lfp SD - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
                                        
%%% T398

        animal_unit_lfp_coh = [];
        animal_unit_lfp_coh_norm = [];

        spike_path = {['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\kilosort\2.26.19'], ...
            ['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\kilosort\2.27.19'], ...
            ['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\kilosort\2.28.19'], ...
            ['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\kilosort\3.1.19'], ...
            ['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\kilosort\3.2.19'], ...
            ['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\kilosort\3.3.19'], ...
            ['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\kilosort\3.4.19'], ...
            ['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\kilosort\3.5.19'], ...
            ['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\kilosort\3.6.19']};

        pre_sleep = {'T398-190226-112203','T398-190227-113133','T398-190228-105336','T398-190301-112846','T398-190302-131102','T398-190303-140804','T398-190304-114807','T398-190305-130625','T398-190306-091252'};
        reach_blocks = {'T398-190226-151524','T398-190227-144600','T398-190228-125951','T398-190301-135817','T398-190302-151642','T398-190303-161420','T398-190304-135101','T398-190305-151059','T398-190306-111957'};
        post_sleep = {'T398-190226-160611','T398-190227-154618','T398-190228-141352','T398-190301-151942','T398-190302-162611','T398-190303-171207','T398-190304-145843','T398-190305-161031','T398-190306-124648'};

        % BANKS
        bank_1_2 = {[5 6 7 8 9 21 39 40 41 42 43 46 51 52 53 54 56 57 58 59 60 61 100 102 105 107 108 110 111 113 114 115 116 119 120 122 123], ...
            [129 130 133 140 141 153 159 160 166 168 170 173 176 179 180 181 187 188 194 195 196 197 198 199 200 201 202 203 205 0206 207 211 212 214 223 224 231]};
        bank_2_1 = {[1 2 3 5 13 19 25 31 32 38 40 42 44 45 48 51 52 53 60 67 68 69 70 71 72 73 74 75 77 78 79 83 84 86 95 96 103], ...
            [135 136 137 139 140 149 167 168 169 170 171 179 180 181 182 184 185 186 187 188 189 228 233 235 236 238 239 241 242 243 244 247 248 249 250 251 255]};

        % BANK ORDER (M1 then DLS)
        bank_order = {[1 2], [1 2], [1 2], [2 1], [1 2], [2 1], [1 2], [1 2], [2 1]};

        load(['D:\Lemke_etal_2021\data\T398\sleep_activity\LFP_coherence.mat'])

        for day = [3:5 7:9]

            %%% LOAD BEH STATE AND START/STOP

            load(['D:\Lemke_etal_2021\data\T398\sleep_activity\day_' num2str(day) '_block_1_sleep_activity.mat'])
            pre_beh_state = beh_state;
            data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\TDT_tanks\M1-DLS_256Probe-190222-144047\' pre_sleep{day}]);
            data_RS4 = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\RS4_data\M1-DLS_256Probe-190222-144047\' pre_sleep{day}],'CHANNEL',1);
            pre_video_start_stop_lfp = [min(data.epocs.PC1_.onset)*(data_RS4.RSn1.fs/20) max(data.epocs.PC1_.onset)*(data_RS4.RSn1.fs/20)];
            pre_video_start_stop = [min(data.epocs.PC1_.onset) max(data.epocs.PC1_.onset)];
            if day == 2
                pre_video_start_stop_lfp(2) = pre_video_start_stop_lfp(2)+(7200*(data_RS4.RSn1.fs/20));
                pre_video_start_stop(2) = pre_video_start_stop(2)+7200;
            end
            pre_beh_state.nrem_interp = pre_beh_state.nrem_interp(1:length(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2)));

            load(['D:\Lemke_etal_2021\data\T398\sleep_activity\day_' num2str(day) '_block_2_sleep_activity.mat'])
            post_beh_state = beh_state;
            data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\TDT_tanks\\M1-DLS_256Probe-190222-144047\' post_sleep{day}]);
            if day == 7
                post_video_start_stop = [10 7210];
                post_video_start_stop_lfp = [10*(data_RS4.RSn1.fs/20) 7210*(data_RS4.RSn1.fs/20)];
            else
                post_video_start_stop = [min(data.epocs.PC1_.onset) max(data.epocs.PC1_.onset)];
                post_video_start_stop_lfp = [min(data.epocs.PC1_.onset)*(data_RS4.RSn1.fs/20) max(data.epocs.PC1_.onset)*(data_RS4.RSn1.fs/20)];
            end
            clearvars data data_RS4 beh_state mean_lfp sleep_rhythms
            post_beh_state.nrem_interp = post_beh_state.nrem_interp(1:length(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2)));

            %%% LOAD LFP

            if bank_order{day}(1)==1
                M1_bank = bank_1_2{1};
                M1_offset = 1;
                DLS_bank = bank_1_2{2};
                DLS_offset = 129;
            elseif bank_order{day}(1)==2
                M1_bank = bank_2_1{2};
                M1_offset = 129;
                DLS_bank = bank_2_1{1};
                DLS_offset = 1;
            end

            M1_mapping = {(M1_offset+[15 16 14 17 13 18 12 19]), ...         %01
                (M1_offset+[11 20 10 21 9 22 8 23]), ...            %02
                [], ...                                             %03
                [], ...                                             %04
                (M1_offset+[63 32 62 33 61 34 60 35]), ...          %05
                (M1_offset+[59 36 58 37 57 38 56 39]), ...          %06
                [], ...                                             %07
                [], ...                                             %08
                (M1_offset+[79 80 78 81 77 82 76 83]), ...          %09
                (M1_offset+[75 84 74 85 73 86 72 87]), ...          %10
                [], ...                                             %11
                [], ...                                             %12
                (M1_offset+[127 96 126 97 125 98 124 99]), ...      %13
                (M1_offset+[123 100 122 101 121 102 120 103]), ...  %14
                [], ...                                             %15
                [], ...                                             %16
                (M1_offset+[111 112 110 113 109 114 108 115]), ...  %17
                (M1_offset+[107 116 106 117 105 118 104 119]), ...  %18
                [], ...                                             %19
                [], ...                                             %20
                (M1_offset+[95 64 94 65 93 66 92 67]), ...          %21
                (M1_offset+[91 68 90 69 89 70 88 71]), ...          %22
                [], ...                                             %23
                [], ...                                             %24
                (M1_offset+[47 48 46 49 45 50 44 51]), ...          %25
                (M1_offset+[43 52 42 53 41 54 40 55]), ...          %26
                [], ...                                             %27
                [], ...                                             %28
                (M1_offset+[31 0 30 1 29 2 28 3]), ...              %29
                (M1_offset+[27 4 26 5 25 6 24 7]), ...              %30
                [], ...                                             %31
                []};                                                %32

            DLS_mapping = {[], ...                                           %01
                [], ...                                             %02
                [], ...                                             %03
                [], ...                                             %04
                [], ...                                             %05
                [], ...                                             %06
                [], ...                                             %07
                [], ...                                             %08
                [], ...                                             %09
                [], ...                                             %10
                [], ...                                             %11
                [], ...                                             %12
                [], ...                                             %13
                [], ...                                             %14
                [], ...                                             %15
                [], ...                                             %16
                (DLS_offset+[15 16 14 17 13 18 12 19]), ...         %17
                (DLS_offset+[11 20 10 21 9 22 8 23]), ...           %18
                (DLS_offset+[7 24 6 25 5 26 4 27]), ...             %19
                (DLS_offset+[3 28 2 29 1 30 0 31]), ...             %20
                (DLS_offset+[63 32 62 33 61 34 60 35]), ...         %21
                (DLS_offset+[59 36 58 37 57 38 56 39]), ...         %22
                (DLS_offset+[55 40 54 41 53 42 52 43]), ...         %23
                (DLS_offset+[51 44 50 45 49 46 48 47]), ...         %24
                (DLS_offset+[79 80 78 81 77 82 76 83]), ...         %25
                (DLS_offset+[75 84 74 85 73 86 72 87]), ...         %26
                (DLS_offset+[71 88 70 89 69 90 68 91]), ...         %27
                (DLS_offset+[67 92 66 93 65 94 64 95]), ...         %28
                (DLS_offset+[127 96 126 97 125 96 124 95]), ...     %29
                (DLS_offset+[123 100 122 101 121 102 120 103]), ... %30
                (DLS_offset+[119 104 118 105 117 106 116 107]), ... %31
                (DLS_offset+[115 108 114 109 113 110 112 111])};    %32

            %%% PRE LFP

            pre_m1_lfp = [];
            for shank = [1 2 5 6 9 10 13 14 17 18 21 22 25 26 29 30]
                disp(shank);
                [tmp_chans, ~] = intersect(M1_mapping{shank},M1_bank);
                for electrode = 1:length(tmp_chans)
                    data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\RS4_data\M1-DLS_256Probe-190222-144047\' pre_sleep{day}],'CHANNEL',tmp_chans(electrode));
                    pre_m1_lfp = [pre_m1_lfp; downsample(data.RSn1.data,20)];
                end
            end

            pre_m1_lfp_norm = [];
            for n = 1:size(pre_m1_lfp,1)
                tmp_lfp = double(pre_m1_lfp(n,:));
                tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                b = find(abs(tmp_lfp)>6);
                tmp_lfp = double(pre_m1_lfp(n,:));
                tmp_lfp(b)=NaN;
                tmp_mean=nanmean(tmp_lfp);
                tmp_sd=nanstd(tmp_lfp);
                tmp_lfp(b)=tmp_mean;
                pre_m1_lfp_norm =  [pre_m1_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
            end
            pre_m1_lfp_norm = pre_m1_lfp_norm-repmat(mean(pre_m1_lfp_norm),[size(pre_m1_lfp_norm,1) 1]);

            %%% POST LFP

            post_m1_lfp = [];
            for shank = [1 2 5 6 9 10 13 14 17 18 21 22 25 26 29 30]
                disp(shank);
                [tmp_chans, ~] = intersect(M1_mapping{shank},M1_bank);
                for electrode = 1:length(tmp_chans)
                    data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\RS4_data\M1-DLS_256Probe-190222-144047\' post_sleep{day}],'CHANNEL',tmp_chans(electrode));
                    post_m1_lfp = [post_m1_lfp; downsample(data.RSn1.data,20)];
                end
            end

            post_m1_lfp_norm = [];
            for n = 1:size(post_m1_lfp,1)
                tmp_lfp = double(post_m1_lfp(n,:));
                tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                b = find(abs(tmp_lfp)>6);
                tmp_lfp = double(post_m1_lfp(n,:));
                tmp_lfp(b)=NaN;
                tmp_mean=nanmean(tmp_lfp);
                tmp_sd=nanstd(tmp_lfp);
                tmp_lfp(b)=tmp_mean;
                post_m1_lfp_norm =  [post_m1_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
            end
            post_m1_lfp_norm = post_m1_lfp_norm-repmat(mean(post_m1_lfp_norm),[size(post_m1_lfp_norm,1) 1]);

            lfp_samp_rate = data.RSn1.fs/20;

            %%% LOAD SPIKING DATA

            opts = delimitedTextImportOptions("NumVariables", 2);
            opts.DataLines = [2, Inf];
            opts.Delimiter = "\t";
            opts.VariableNames = ["cluster_id", "group"];
            opts.VariableTypes = ["double", "categorical"];
            opts = setvaropts(opts, 2, "EmptyFieldRule", "auto");
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";

            cd([spike_path{day} '\DLS']);
            dls_ts = readNPY('spike_times.npy');
            dls_id = readNPY('spike_clusters.npy');
            dls_unit_class = readtable('cluster_group.tsv', opts);
            dls_clust_id = zeros(1,max(dls_unit_class{:,1}));
            dls_clust_id_idx = dls_unit_class{:,1}'+1;
            idx_count = 1;
            for clust_id = dls_clust_id_idx
                tmp_qual = cellstr(dls_unit_class{idx_count,2});
                if strcmp(tmp_qual{1,1},'noise')
                    dls_clust_id(clust_id) = 0;
                else
                    dls_clust_id(clust_id) = 1;
                end
                idx_count = idx_count + 1;
            end
            dls_clust_id = find(dls_clust_id)-1;

            %%% CHECK UNIT %%%
            pre_sleep_length = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\RS4_data\M1-DLS_256Probe-190222-144047\' pre_sleep{day}],'CHANNEL',1);
            samp_freq = pre_sleep_length.RSn1.fs;
            pre_sleep_length = length(pre_sleep_length.RSn1.data);

            reach_length = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\RS4_data\M1-DLS_256Probe-190222-144047\' reach_blocks{day}],'CHANNEL',1);
            reach_length = length(reach_length.RSn1.data);

            post_sleep_length = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\RS4_data\M1-DLS_256Probe-190222-144047\' post_sleep{day}],'CHANNEL',1);
            post_sleep_length = length(post_sleep_length.RSn1.data);

            max(dls_ts)
            pre_sleep_length + reach_length + post_sleep_length
            %%% CHECK UNIT %%%

            pre_dls_spiking = cell(length(dls_clust_id),1);
            post_dls_spiking = cell(length(dls_clust_id),1);
            dls_count = 1;
            for dls_n = dls_clust_id
                tmp_spikes = dls_ts(dls_id==dls_n);
                pre_dls_spiking{dls_count,1} = double(tmp_spikes(tmp_spikes<pre_sleep_length))/samp_freq;
                post_dls_spiking{dls_count,1} = double(tmp_spikes(tmp_spikes>(pre_sleep_length+reach_length))-(pre_sleep_length+reach_length))/samp_freq;
                dls_count = dls_count+1;
            end

            %%% PRE AND POST SLEEP NREM COH

            pre_coh = [];
            post_coh = [];
            for pairs=1:1406
                pre_coh = [pre_coh mean(pre_LFP_coh{pairs}(day-2,40:95))];
                post_coh = [post_coh mean(post_LFP_coh{pairs}(day-2,40:95))];
            end
            pre_coh = mean(reshape(pre_coh,38,37));
            post_coh = mean(reshape(post_coh,38,37));

            %%% CALCULATE SPIKE COH

            coh_coh = [];

            for m1_lfp = 1:37

                tmp_pre_m1_lfp = pre_m1_lfp_norm(m1_lfp,:);
                tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2));
                tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
                tmp_pre_m1_lfp = eegfilt(tmp_pre_m1_lfp,lfp_samp_rate,4,8);
                tmp_pre_m1_lfp = angle(hilbert(tmp_pre_m1_lfp));

                tmp_post_m1_lfp = post_m1_lfp_norm(m1_lfp,:);
                tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2));
                tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
                tmp_post_m1_lfp = eegfilt(tmp_post_m1_lfp,lfp_samp_rate,4,8);
                tmp_post_m1_lfp = angle(hilbert(tmp_post_m1_lfp));

                tmp_chan_plv = [];

                for dls_unit = 1:length(pre_dls_spiking)

                    pre_dls_unit = pre_dls_spiking{dls_unit};%(pre_dls_spiking{dls_unit}>pre_video_start_stop(1) & pre_dls_spiking{dls_unit}<pre_video_start_stop(2))-pre_video_start_stop(1);
                    pre_dls_unit = histc(pre_dls_unit,linspace(pre_video_start_stop(1),pre_video_start_stop(2),length(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2))));
                    pre_dls_unit = pre_dls_unit(pre_beh_state.nrem_interp==1);

                    post_dls_unit = post_dls_spiking{dls_unit};%(post_dls_spiking{dls_unit}>post_video_start_stop(1) & post_dls_spiking{dls_unit}<post_video_start_stop(2))-post_video_start_stop(1);
                    post_dls_unit = histc(post_dls_unit,linspace(post_video_start_stop(1),post_video_start_stop(2),length(post_video_start_stop_lfp(1):post_video_start_stop_lfp(2))));
                    post_dls_unit = post_dls_unit(post_beh_state.nrem_interp==1);
                    
                    [pval, z] = circ_rtest([tmp_pre_m1_lfp(pre_dls_unit==1) tmp_post_m1_lfp(post_dls_unit==1)]);
                    sd = circ_std([tmp_pre_m1_lfp(pre_dls_unit==1) tmp_post_m1_lfp(post_dls_unit==1)]');
                    tmp_chan_plv = [tmp_chan_plv; pval z sd pre_coh(m1_lfp)];

                end

                coh_coh = [coh_coh; mean(tmp_chan_plv)];

            end

            animal_unit_lfp_coh = [animal_unit_lfp_coh; coh_coh];
            animal_unit_lfp_coh_norm = [animal_unit_lfp_coh_norm; coh_coh(:,1) zscore(coh_coh(:,2)) zscore(coh_coh(:,3)) zscore(coh_coh(:,4))]; 

        end

        all_unit_coh_lfp_coh = [all_unit_coh_lfp_coh; animal_unit_lfp_coh];
        all_unit_coh_lfp_coh_day_norm = [all_unit_coh_lfp_coh_day_norm; animal_unit_lfp_coh_norm];
        
        figure;
            subplot(2,2,1); hold on;
                scatter(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,2),30,[0 0 0],'filled')
                [R,P] = corrcoef(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,2));
                [p] = polyfit(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,2),1);
                x=linspace(min(animal_unit_lfp_coh(:,4)),max(animal_unit_lfp_coh(:,4)),10);
                xlim([min(animal_unit_lfp_coh(:,4)) max(animal_unit_lfp_coh(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T398 - lfp coh vs. dls unit m1 lfp z val - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,2); hold on;
                scatter(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,2),30,[0 0 0],'filled');
                [R,P] = corrcoef(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,2));
                [p] = polyfit(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,2),1);
                x=linspace(min(animal_unit_lfp_coh_norm(:,4)),max(animal_unit_lfp_coh_norm(:,4)),10);
                xlim([min(animal_unit_lfp_coh_norm(:,4)) max(animal_unit_lfp_coh_norm(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T398 - norm lfp coh vs. norm dls unit m1 lfp z val - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,3); hold on;
                scatter(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,3),30,[0 0 0],'filled')
                [R,P] = corrcoef(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,3));
                [p] = polyfit(animal_unit_lfp_coh(:,4),animal_unit_lfp_coh(:,3),1);
                x=linspace(min(animal_unit_lfp_coh(:,4)),max(animal_unit_lfp_coh(:,4)),10);
                xlim([min(animal_unit_lfp_coh(:,4)) max(animal_unit_lfp_coh(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T398 - lfp coh vs. dls unit m1 lfp SD - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,4); hold on;
                scatter(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,3),30,[0 0 0],'filled');
                [R,P] = corrcoef(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,3));
                [p] = polyfit(animal_unit_lfp_coh_norm(:,4),animal_unit_lfp_coh_norm(:,3),1);
                x=linspace(min(animal_unit_lfp_coh_norm(:,4)),max(animal_unit_lfp_coh_norm(:,4)),10);
                xlim([min(animal_unit_lfp_coh_norm(:,4)) max(animal_unit_lfp_coh_norm(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['T398 - norm lfp coh vs. norm dls unit m1 lfp SD - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
                           
%%% PLOT ALL

        figure;
            subplot(2,2,1); hold on;
                scatter(all_unit_coh_lfp_coh(:,4),all_unit_coh_lfp_coh(:,2),30,[0 0 0],'filled')
                [R,P] = corrcoef(all_unit_coh_lfp_coh(:,4),all_unit_coh_lfp_coh(:,2));
                [p] = polyfit(all_unit_coh_lfp_coh(:,4),all_unit_coh_lfp_coh(:,2),1);
                x=linspace(min(all_unit_coh_lfp_coh(:,4)),max(all_unit_coh_lfp_coh(:,4)),10);
                xlim([min(all_unit_coh_lfp_coh(:,4)) max(all_unit_coh_lfp_coh(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['ALL ANIMAL - lfp coh vs. dls unit m1 lfp z val - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,2); hold on;
                scatter(all_unit_coh_lfp_coh_day_norm(:,4),all_unit_coh_lfp_coh_day_norm(:,2),30,[0 0 0],'filled');
                [R,P] = corrcoef(all_unit_coh_lfp_coh_day_norm(:,4),all_unit_coh_lfp_coh_day_norm(:,2));
                [p] = polyfit(all_unit_coh_lfp_coh_day_norm(:,4),all_unit_coh_lfp_coh_day_norm(:,2),1);
                x=linspace(min(all_unit_coh_lfp_coh_day_norm(:,4)),max(all_unit_coh_lfp_coh_day_norm(:,4)),10);
                xlim([min(all_unit_coh_lfp_coh_day_norm(:,4)) max(all_unit_coh_lfp_coh_day_norm(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['ALL ANIMAL - norm lfp coh vs. norm dls unit m1 lfp z val - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,3); hold on;
                scatter(all_unit_coh_lfp_coh(:,4),all_unit_coh_lfp_coh(:,3),30,[0 0 0],'filled')
                [R,P] = corrcoef(all_unit_coh_lfp_coh(:,4),all_unit_coh_lfp_coh(:,3));
                [p] = polyfit(all_unit_coh_lfp_coh(:,4),all_unit_coh_lfp_coh(:,3),1);
                x=linspace(min(all_unit_coh_lfp_coh(:,4)),max(all_unit_coh_lfp_coh(:,4)),10);
                xlim([min(all_unit_coh_lfp_coh(:,4)) max(all_unit_coh_lfp_coh(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['ALL ANIMAL - lfp coh vs. dls unit m1 lfp SD - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
            subplot(2,2,4); hold on;
                scatter(all_unit_coh_lfp_coh_day_norm(:,4),all_unit_coh_lfp_coh_day_norm(:,3),30,[0 0 0],'filled');
                [R,P] = corrcoef(all_unit_coh_lfp_coh_day_norm(:,4),all_unit_coh_lfp_coh_day_norm(:,3));
                [p] = polyfit(all_unit_coh_lfp_coh_day_norm(:,4),all_unit_coh_lfp_coh_day_norm(:,3),1);
                x=linspace(min(all_unit_coh_lfp_coh_day_norm(:,4)),max(all_unit_coh_lfp_coh_day_norm(:,4)),10);
                xlim([min(all_unit_coh_lfp_coh_day_norm(:,4)) max(all_unit_coh_lfp_coh_day_norm(:,4))]);
                y = x*p(1)+p(2);
                plot(x,y,'color','r','LineStyle','--','LineWidth',3);
                title(['ALL ANIMAL - norm lfp coh vs. norm dls unit m1 lfp SD - R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
 
%% NEED RAW DATA NOT ON DRYAD | sFIG 1 | 4-8HZ LFP PHASE DIFFERENCE

    %%% INVIDUAL EXAMPLE

        day_blocks = {[2 3 7 8 9 11 12 13], ...
            [14 16 18 19 20 21], ...
            [22 24 25 26 27 29], ...
            [30 31 32 33 34 35 36 37 38], ...
            [39 43 44 45 47], ...
            [48 49 51 53 54 56], ...
            [57 58 59 60 61], ...
            [62 63 64 65 66]};
        pre_sleep_blocks = {[2],[1],[1],[1],[1],[1],[1],[1]};
        reach_blocks = {[3 4 6 7], [2 3 5],[2 3 6],[2 3 5 6 7 9],[2 3 5],[2 3 5],[2 3 5]};
        post_sleep_blocks = {[8],[4],[5],[8],[4],[5],[4],[4]};
        plx_files = {'T102_blocks_2_3_7_8_9_11_12_13-01','T102_blocks_14_16_18_19_20_21-01','T102_blocks_22_24_25_26_27_29-01','T102_blocks_30_31_32_33_34_35_36_37_38-01','T102_blocks_39_43_44_45_47-01','T102_blocks_48_49_51_53_54_56-01','T102_blocks_57_58_59_60_61-01','T102_blocks_62_63_64_65_66-01'};

        for day = 7
            tmp_day_blocks = day_blocks{day};

            %%% LOAD BEH STATE, LFP, AND START/STOP

                load(['D:\Lemke_etal_2021\data\T102\sleep_rhythms\day_' num2str(day) '\block_' num2str(tmp_day_blocks(pre_sleep_blocks{day})) '\sleep_classification_rhythms.mat'])
                pre_beh_state = beh_state;
                data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\Block-' num2str(tmp_day_blocks(pre_sleep_blocks{day}))],'Type',4);
                pre_lfp = data.streams.LFPs;
                pre_video_start_stop_lfp = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
                pre_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];
                pre_lfp = pre_lfp.data;

                pre_lfp_norm = [];
                for n = 1:32
                    tmp_lfp = pre_lfp(n,:);
                    tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                    b = find(abs(tmp_lfp)>6);
                    tmp_lfp = pre_lfp(n,:);
                    tmp_lfp(b)=NaN;
                    tmp_mean=nanmean(tmp_lfp);
                    tmp_sd=nanstd(tmp_lfp);
                    tmp_lfp(b)=tmp_mean;
                    pre_lfp_norm =  [pre_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
                end
                pre_lfp_norm(1:16,:) = pre_lfp_norm(1:16,:)-repmat(mean(pre_lfp_norm(1:16,:),1),[size(pre_lfp_norm(1:16,:),1) 1]);
                pre_lfp_norm(17:32,:) = pre_lfp_norm(17:32,:)-repmat(mean(pre_lfp_norm(17:32,:),1),[size(pre_lfp_norm(17:32,:),1) 1]);

            %%% LOAD SPIKING DATA

                plx_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\' plx_files{day} '.txt']);
                plx_spiking = table2array(plx_spiking);

                tmp_ts = cell(32,length(tmp_day_blocks));
                block_count = 1;
                for block = tmp_day_blocks
                    tmp_num = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\Block-' num2str(block)],'Type',3);
                    for chan = 1:32
                        tmp_ts{chan,block_count} = tmp_num.snips.eNeu.ts(tmp_num.snips.eNeu.chan==chan);
                    end
                    block_count = block_count + 1;
                end

                tmp_length = zeros(32,length(tmp_day_blocks));
                for chan = 1:32
                    for block = 1:length(tmp_day_blocks)
                        tmp_length(chan,block) = length(tmp_ts{chan,block});
                    end
                end
                tmp_length = [zeros(32,1) tmp_length];

                tmp_sort = cell(32,length(tmp_day_blocks));
                for chan = 1:32
                    chan_sorts = plx_spiking(plx_spiking(:,1)==chan,:);
                    for block = 1:length(tmp_day_blocks)
                        tmp_sort{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),2);
                    end
                end

                pre_dls_spiking = cell(16,max(plx_spiking(:,2)));
                post_dls_spiking = cell(16,max(plx_spiking(:,2)));
                count = 1;
                for dls_chan = 1:16
                    for dls_unit = 1:max(plx_spiking(:,2))
                        pre_dls_spiking{count,dls_unit} = tmp_ts{dls_chan,pre_sleep_blocks{day}}(tmp_sort{dls_chan,pre_sleep_blocks{day}}==dls_unit);
                        post_dls_spiking{count,dls_unit} = tmp_ts{dls_chan,post_sleep_blocks{day}}(tmp_sort{dls_chan,post_sleep_blocks{day}}==dls_unit);
                    end
                    count = count + 1;
                end

            %%% M1-DLS LFP COH

                phase_diff = [];
                for m1_lfp = 1:16
                    for dls_lfp = 17:32

                        tmp_pre_m1_lfp = pre_lfp_norm(m1_lfp,:);
                        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2));
                        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
                        tmp_pre_m1_lfp_filt = eegfilt(tmp_pre_m1_lfp,data.streams.LFPs.fs,4,8);
                        tmp_pre_m1_lfp_phase = angle(hilbert(tmp_pre_m1_lfp_filt));

                        tmp_pre_dls_lfp = pre_lfp_norm(dls_lfp,:);
                        tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2));
                        tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_beh_state.nrem_interp==1);
                        tmp_pre_dls_lfp_filt = eegfilt(tmp_pre_dls_lfp,data.streams.LFPs.fs,4,8);
                        tmp_pre_dls_lfp_phase = angle(hilbert(tmp_pre_dls_lfp_filt));

                        pair_phase_diff = [tmp_pre_dls_lfp_phase-tmp_pre_m1_lfp_phase];
                        pair_phase_diff(pair_phase_diff>pi) = pair_phase_diff(pair_phase_diff>pi)-(2*pi);
                        pair_phase_diff(pair_phase_diff<-pi) = pair_phase_diff(pair_phase_diff<-pi)+(2*pi);

                        [m, i] = max(histc(pair_phase_diff,[-pi:pi/24:pi]));
                        tmp_x = [-pi:pi/24:pi];
                        phase_diff = [phase_diff; m tmp_x(i)];

                    end
                end

                %     tmp = phase_diff(:,2)<0 & phase_diff(:,2)>-.5;
                %     tmp_phase_diff = phase_diff(:,1);
                %     tmp_phase_diff(tmp==0) = 0;
                %     [m,i] = sort(tmp_phase_diff)

                index = find(phase_diff(:,1) > prctile(phase_diff(:,1),75));
                index = index(23);
                count = 1;
                for m1_lfp = 1:16
                    for dls_lfp = 17:32
                        if ismember(count,index)
                            count = count + 1;
                        else
                            count = count + 1;
                            continue;
                        end

                        tmp_pre_m1_lfp = pre_lfp_norm(m1_lfp,:);
                        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2));
                        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
                        tmp_pre_m1_lfp_filt = eegfilt(tmp_pre_m1_lfp,data.streams.LFPs.fs,4,8);
                        tmp_pre_m1_lfp_phase = angle(hilbert(tmp_pre_m1_lfp_filt));

                        tmp_pre_dls_lfp = pre_lfp_norm(dls_lfp,:);
                        tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_video_start_stop_lfp(1):pre_video_start_stop_lfp(2));
                        tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_beh_state.nrem_interp==1);
                        tmp_pre_dls_lfp_filt = eegfilt(tmp_pre_dls_lfp,data.streams.LFPs.fs,4,8);
                        tmp_pre_dls_lfp_phase = angle(hilbert(tmp_pre_dls_lfp_filt));

                        figure;
                        ax1 = subplot(3,1,1); hold on;
                        plot(tmp_pre_m1_lfp(1:10000),'color',[1 0 0])
                        plot(tmp_pre_dls_lfp(1:10000),'color',[0 0 1])
                        ax2 = subplot(3,1,2); hold on;
                        plot(tmp_pre_m1_lfp_filt(1:10000),'color',[1 0 0])
                        plot(tmp_pre_dls_lfp_filt(1:10000),'color',[0 0 1])
                        ax3 = subplot(3,1,3); hold on;
                        plot(tmp_pre_m1_lfp_phase(1:10000),'color',[1 0 0])
                        plot(tmp_pre_dls_lfp_phase(1:10000),'color',[0 0 1])
                        linkaxes([ax1 ax2 ax3],'x')
                        xlim([2000 2000+(5*data.streams.LFPs.fs)])
                        %             print('-painters','-depsc','C:\Users\Stefan\Desktop\Lemke_etal_2021\paper\Figure_2\M1_DLS_LFP_phase_diff_example.esp')

                        pair_phase_diff = [tmp_pre_dls_lfp_phase-tmp_pre_m1_lfp_phase];
                        pair_phase_diff(pair_phase_diff>pi) = pair_phase_diff(pair_phase_diff>pi)-(2*pi);
                        pair_phase_diff(pair_phase_diff<-pi) = pair_phase_diff(pair_phase_diff<-pi)+(2*pi);

                        figure;
                        histogram(pair_phase_diff,[-pi:pi/24:pi],'DisplayStyle','Stairs','Normalization','Probability')
                        xlim([-pi pi])
                        %             print('-painters','-depsc','C:\Users\Stefan\Desktop\Lemke_etal_2021\paper\Figure_2\M1_DLS_LFP_phase_diff_histogram_example.esp')
                end
            end
        end

    %%% ALL ANIMAL MAX PHASE DIFF

        all_pair_max_phase_diff = [];

        %%% T102

            %%% FIND INCREASE CHANNELS

                load([data_path '\data\T102\sleep_activity\LFP_coherence.mat'])

                increase_channels = [];
                for pairs=1:256
                    if mean(mean([pre_LFP_coh{pairs}(8,48:113); post_LFP_coh{pairs}(8,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)])))
                        increase_channels = [increase_channels 1];
                    else
                        increase_channels = [increase_channels 0];
                    end
                end

            %%% FIND HIGH LFP COH PHASE DIFF

                day_blocks = {[2 3 7 8 9 11 12 13], ...
                    [14 16 18 19 20 21], ...
                    [22 24 25 26 27 29], ...
                    [30 31 32 33 34 35 36 37 38], ...
                    [39 43 44 45 47], ...
                    [48 49 51 53 54 56], ...
                    [57 58 59 60 61], ...
                    [62 63 64 65 66]};

                pre_sleep_blocks = {[2],[1],[1],[1],[1],[1],[1],[1]};
                reach_blocks = {[3 4 6 7], [2 3 5],[2 3 6],[2 3 5 6 7 9],[2 3 5],[2 3 5],[2 3 5]};
                post_sleep_blocks = {[8],[4],[5],[8],[4],[5],[4],[4]};

                for day = 1:8

                    tmp_day_blocks = day_blocks{day};

                    %%% LOAD BEH STATE, LFP, AND START/STOP

                    load(['D:\Lemke_etal_2021\data\T102\sleep_rhythms\day_' num2str(day) '\block_' num2str(tmp_day_blocks(pre_sleep_blocks{day})) '\sleep_classification_rhythms.mat'])
                    pre_beh_state = beh_state;
                    data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\Block-' num2str(tmp_day_blocks(pre_sleep_blocks{day}))],'Type',4);
                    pre_lfp = data.streams.LFPs;
                    pre_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
                    pre_lfp = pre_lfp.data;

                    pre_lfp_norm = [];
                    for n = 1:size(pre_lfp,1)
                        tmp_lfp = pre_lfp(n,:);
                        tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                        b = find(abs(tmp_lfp)>6);
                        tmp_lfp = pre_lfp(n,:);
                        tmp_lfp(b)=NaN;
                        tmp_mean=nanmean(tmp_lfp);
                        tmp_sd=nanstd(tmp_lfp);
                        tmp_lfp(b)=tmp_mean;
                        pre_lfp_norm =  [pre_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
                    end
                    pre_lfp_norm(1:16,:) = pre_lfp_norm(1:16,:)-repmat(mean(pre_lfp_norm(1:16,:),1),[size(pre_lfp_norm(1:16,:),1) 1]);
                    pre_lfp_norm(17:32,:) = pre_lfp_norm(17:32,:)-repmat(mean(pre_lfp_norm(17:32,:),1),[size(pre_lfp_norm(17:32,:),1) 1]);

                    load(['D:\Lemke_etal_2021\data\T102\sleep_rhythms\day_' num2str(day) '\block_' num2str(tmp_day_blocks(post_sleep_blocks{day})) '\sleep_classification_rhythms.mat'])
                    post_beh_state = beh_state;
                    data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\Block-' num2str(tmp_day_blocks(post_sleep_blocks{day}))],'Type',4);
                    post_lfp = data.streams.LFPs;
                    post_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
                    post_lfp = post_lfp.data;

                    post_lfp_norm = [];
                    for n = 1:size(post_lfp,1)
                        tmp_lfp = post_lfp(n,:);
                        tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                        b = find(abs(tmp_lfp)>6);
                        tmp_lfp = post_lfp(n,:);
                        tmp_lfp(b)=NaN;
                        tmp_mean=nanmean(tmp_lfp);
                        tmp_sd=nanstd(tmp_lfp);
                        tmp_lfp(b)=tmp_mean;
                        post_lfp_norm =  [post_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
                    end
                    post_lfp_norm(1:16,:) = post_lfp_norm(1:16,:)-repmat(mean(post_lfp_norm(1:16,:),1),[size(post_lfp_norm(1:16,:),1) 1]);
                    post_lfp_norm(17:32,:) = post_lfp_norm(17:32,:)-repmat(mean(post_lfp_norm(17:32,:),1),[size(post_lfp_norm(17:32,:),1) 1]);

                    pair_count = 1;
                    for m1_chan = 1:16

                        tmp_pre_m1_lfp = pre_lfp_norm(m1_chan,:);
                        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
                        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
                        tmp_pre_m1_lfp = eegfilt(tmp_pre_m1_lfp(1:data.streams.LFPs.fs*300),data.streams.LFPs.fs,4,8);
                        tmp_pre_m1_lfp = angle(hilbert(tmp_pre_m1_lfp));

                        tmp_post_m1_lfp = post_lfp_norm(m1_chan,:);
                        tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop(1):post_video_start_stop(2));
                        tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
                        tmp_post_m1_lfp = eegfilt(tmp_post_m1_lfp(1:data.streams.LFPs.fs*300),data.streams.LFPs.fs,4,8);
                        tmp_post_m1_lfp = angle(hilbert(tmp_post_m1_lfp));

                        for dls_chan = 17:32

                            if increase_channels(pair_count)==0
                                pair_count = pair_count + 1;
                                continue
                            end

                            pair_count = pair_count + 1;

                            tmp_pre_dls_lfp = pre_lfp_norm(dls_chan,:);
                            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
                            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_beh_state.nrem_interp==1);
                            tmp_pre_dls_lfp = eegfilt(tmp_pre_dls_lfp(1:data.streams.LFPs.fs*300),data.streams.LFPs.fs,4,8);
                            tmp_pre_dls_lfp = angle(hilbert(tmp_pre_dls_lfp));

                            tmp_post_dls_lfp = post_lfp_norm(dls_chan,:);
                            tmp_post_dls_lfp = tmp_post_dls_lfp(post_video_start_stop(1):post_video_start_stop(2));
                            tmp_post_dls_lfp = tmp_post_dls_lfp(post_beh_state.nrem_interp==1);
                            tmp_post_dls_lfp = eegfilt(tmp_post_dls_lfp(1:data.streams.LFPs.fs*300),data.streams.LFPs.fs,4,8);
                            tmp_post_dls_lfp = angle(hilbert(tmp_post_dls_lfp));

                            pair_phase_diff = [tmp_pre_dls_lfp-tmp_pre_m1_lfp tmp_post_dls_lfp-tmp_post_m1_lfp];
                            pair_phase_diff(pair_phase_diff>pi) = pair_phase_diff(pair_phase_diff>pi)-(2*pi);
                            pair_phase_diff(pair_phase_diff<-pi) = pair_phase_diff(pair_phase_diff<-pi)+(2*pi);

                            [~, i] = max(histc(pair_phase_diff,[-pi:pi/24:pi]));
                            tmp_x = [-pi:pi/24:pi];
                            all_pair_max_phase_diff = [all_pair_max_phase_diff tmp_x(i)];

                        end
                    end
                end

        %%% T107

            %%% FIND INCREASE CHANNELS

                load([data_path '\data\T107\sleep_activity\LFP_coherence.mat'])

                increase_channels = [];
                for pairs=1:256
                    if mean(mean([pre_LFP_coh{pairs}(6,48:113); post_LFP_coh{pairs}(6,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)])))
                        increase_channels = [increase_channels 1];
                    else
                        increase_channels = [increase_channels 0];
                    end
                end

            %%% FIND HIGH LFP COH PHASE DIFF

                day_blocks = {[1 2],[3 4 5 6 8],[9 10 11 12 13 14 15],[16 17 18 19 20 21 22 24 25],[26 27 30 31 36 39],[41 42 43 47 49 51 52]};

                pre_sleep_blocks = {[1],[1],[1],[3],[1],[1]};
                post_sleep_blocks = {[2],[4],[6],[8],[5],[4]};

                for day = 1:6

                    tmp_day_blocks = day_blocks{day};

                    %%% LOAD BEH STATE, LFP, AND START/STOP

                    load(['D:\Lemke_etal_2021\data\T107\sleep_rhythms\day_' num2str(day) '\block_' num2str(tmp_day_blocks(pre_sleep_blocks{day})) '\sleep_classification_rhythms.mat'])
                    pre_beh_state = beh_state;
                    data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T107\neural_data\Block-' num2str(tmp_day_blocks(pre_sleep_blocks{day}))],'Type',4);
                    pre_lfp = data.streams.LFPs;
                    pre_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
                    pre_lfp = pre_lfp.data;

                    pre_lfp_norm = [];
                    for n = 1:size(pre_lfp,1)
                        tmp_lfp = pre_lfp(n,:);
                        tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                        b = find(abs(tmp_lfp)>6);
                        tmp_lfp = pre_lfp(n,:);
                        tmp_lfp(b)=NaN;
                        tmp_mean=nanmean(tmp_lfp);
                        tmp_sd=nanstd(tmp_lfp);
                        tmp_lfp(b)=tmp_mean;
                        pre_lfp_norm =  [pre_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
                    end
                    pre_lfp_norm(1:16,:) = pre_lfp_norm(1:16,:)-repmat(mean(pre_lfp_norm(1:16,:),1),[size(pre_lfp_norm(1:16,:),1) 1]);
                    pre_lfp_norm(17:32,:) = pre_lfp_norm(17:32,:)-repmat(mean(pre_lfp_norm(17:32,:),1),[size(pre_lfp_norm(17:32,:),1) 1]);

                    load(['D:\Lemke_etal_2021\data\T107\sleep_rhythms\day_' num2str(day) '\block_' num2str(tmp_day_blocks(post_sleep_blocks{day})) '\sleep_classification_rhythms.mat'])
                    post_beh_state = beh_state;
                    data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T107\neural_data\Block-' num2str(tmp_day_blocks(post_sleep_blocks{day}))],'Type',4);
                    post_lfp = data.streams.LFPs;
                    post_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
                    post_lfp = post_lfp.data;

                    post_lfp_norm = [];
                    for n = 1:size(post_lfp,1)
                        tmp_lfp = post_lfp(n,:);
                        tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                        b = find(abs(tmp_lfp)>6);
                        tmp_lfp = post_lfp(n,:);
                        tmp_lfp(b)=NaN;
                        tmp_mean=nanmean(tmp_lfp);
                        tmp_sd=nanstd(tmp_lfp);
                        tmp_lfp(b)=tmp_mean;
                        post_lfp_norm =  [post_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
                    end
                    post_lfp_norm(1:16,:) = post_lfp_norm(1:16,:)-repmat(mean(post_lfp_norm(1:16,:),1),[size(post_lfp_norm(1:16,:),1) 1]);
                    post_lfp_norm(17:32,:) = post_lfp_norm(17:32,:)-repmat(mean(post_lfp_norm(17:32,:),1),[size(post_lfp_norm(17:32,:),1) 1]);

                    pair_count = 1;
                    for m1_chan = 1:16

                        tmp_pre_m1_lfp = pre_lfp_norm(m1_chan,:);
                        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
                        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
                        tmp_pre_m1_lfp = eegfilt(tmp_pre_m1_lfp(1:data.streams.LFPs.fs*300),data.streams.LFPs.fs,4,8);
                        tmp_pre_m1_lfp = angle(hilbert(tmp_pre_m1_lfp));

                        tmp_post_m1_lfp = post_lfp_norm(m1_chan,:);
                        tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop(1):post_video_start_stop(2));
                        tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
                        tmp_post_m1_lfp = eegfilt(tmp_post_m1_lfp(1:data.streams.LFPs.fs*300),data.streams.LFPs.fs,4,8);
                        tmp_post_m1_lfp = angle(hilbert(tmp_post_m1_lfp));

                        for dls_chan = 17:32

                            if increase_channels(pair_count)==0
                                pair_count = pair_count + 1;
                                continue
                            end

                            pair_count = pair_count + 1;

                            tmp_pre_dls_lfp = pre_lfp_norm(dls_chan,:);
                            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
                            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_beh_state.nrem_interp==1);
                            tmp_pre_dls_lfp = eegfilt(tmp_pre_dls_lfp(1:data.streams.LFPs.fs*300),data.streams.LFPs.fs,4,8);
                            tmp_pre_dls_lfp = angle(hilbert(tmp_pre_dls_lfp));

                            tmp_post_dls_lfp = post_lfp_norm(dls_chan,:);
                            tmp_post_dls_lfp = tmp_post_dls_lfp(post_video_start_stop(1):post_video_start_stop(2));
                            tmp_post_dls_lfp = tmp_post_dls_lfp(post_beh_state.nrem_interp==1);
                            tmp_post_dls_lfp = eegfilt(tmp_post_dls_lfp(1:data.streams.LFPs.fs*300),data.streams.LFPs.fs,4,8);
                            tmp_post_dls_lfp = angle(hilbert(tmp_post_dls_lfp));

                            pair_phase_diff = [tmp_pre_dls_lfp-tmp_pre_m1_lfp tmp_post_dls_lfp-tmp_post_m1_lfp];
                            pair_phase_diff(pair_phase_diff>pi) = pair_phase_diff(pair_phase_diff>pi)-(2*pi);
                            pair_phase_diff(pair_phase_diff<-pi) = pair_phase_diff(pair_phase_diff<-pi)+(2*pi);

                            [~, i] = max(histc(pair_phase_diff,[-pi:pi/24:pi]));
                            tmp_x = [-pi:pi/24:pi];
                            all_pair_max_phase_diff = [all_pair_max_phase_diff tmp_x(i)];

                        end
                    end
                end

        %%% T200

            %%% FIND INCREASE CHANNELS

                load([data_path '\data\T200\sleep_activity\LFP_coherence.mat'])

                increase_channels = [];
                for pairs=1:256
                    if mean(mean([pre_LFP_coh{pairs}(9,48:113); post_LFP_coh{pairs}(9,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)])))
                        increase_channels = [increase_channels 1];
                    else
                        increase_channels = [increase_channels 0];
                    end
                end

            %%% FIND HIGH LFP COH PHASE DIFF

                day_blocks = {[1],[2 7 8],[9 11 12 13 14],[15 16 17 18 19],[20 21 22 23 24],[25 26 27],[28 29 30],[32 33 34 35 36],[37 38 39 40 41],[42 43 44 45 46 47]};

                pre_sleep_blocks = {[1],[1],[2],[1],[1],[1],[1],[1],[1],[1]};
                post_sleep_blocks = {[],[3],[5],[5],[4],[3],[3],[5],[5],[6]};

                for day = 2:10

                    tmp_day_blocks = day_blocks{day};

                    %%% LOAD BEH STATE, LFP, AND START/STOP

                    load(['D:\Lemke_etal_2021\data\T200\sleep_rhythms\day_' num2str(day) '\block_' num2str(tmp_day_blocks(pre_sleep_blocks{day})) '\sleep_classification_rhythms.mat'])
                    pre_beh_state = beh_state;
                    data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\Block-' num2str(tmp_day_blocks(pre_sleep_blocks{day}))],'Type',4);
                    pre_lfp = data.streams.LFPs;
                    pre_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
                    pre_lfp = pre_lfp.data;

                    pre_lfp_norm = [];
                    for n = 1:size(pre_lfp,1)
                        tmp_lfp = pre_lfp(n,:);
                        tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                        b = find(abs(tmp_lfp)>6);
                        tmp_lfp = pre_lfp(n,:);
                        tmp_lfp(b)=NaN;
                        tmp_mean=nanmean(tmp_lfp);
                        tmp_sd=nanstd(tmp_lfp);
                        tmp_lfp(b)=tmp_mean;
                        pre_lfp_norm =  [pre_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
                    end
                    pre_lfp_norm(1:16,:) = pre_lfp_norm(1:16,:)-repmat(mean(pre_lfp_norm(1:16,:),1),[size(pre_lfp_norm(1:16,:),1) 1]);
                    pre_lfp_norm(17:32,:) = pre_lfp_norm(17:32,:)-repmat(mean(pre_lfp_norm(17:32,:),1),[size(pre_lfp_norm(17:32,:),1) 1]);

                    load(['D:\Lemke_etal_2021\data\T200\sleep_rhythms\day_' num2str(day) '\block_' num2str(tmp_day_blocks(post_sleep_blocks{day})) '\sleep_classification_rhythms.mat'])
                    post_beh_state = beh_state;
                    data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\Block-' num2str(tmp_day_blocks(post_sleep_blocks{day}))],'Type',4);
                    post_lfp = data.streams.LFPs;
                    post_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
                    post_lfp = post_lfp.data;

                    post_lfp_norm = [];
                    for n = 1:size(post_lfp,1)
                        tmp_lfp = post_lfp(n,:);
                        tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                        b = find(abs(tmp_lfp)>6);
                        tmp_lfp = post_lfp(n,:);
                        tmp_lfp(b)=NaN;
                        tmp_mean=nanmean(tmp_lfp);
                        tmp_sd=nanstd(tmp_lfp);
                        tmp_lfp(b)=tmp_mean;
                        post_lfp_norm =  [post_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
                    end
                    post_lfp_norm(1:16,:) = post_lfp_norm(1:16,:)-repmat(mean(post_lfp_norm(1:16,:),1),[size(post_lfp_norm(1:16,:),1) 1]);
                    post_lfp_norm(17:32,:) = post_lfp_norm(17:32,:)-repmat(mean(post_lfp_norm(17:32,:),1),[size(post_lfp_norm(17:32,:),1) 1]);

                    pair_count = 1;
                    for m1_chan = 17:32

                        tmp_pre_m1_lfp = pre_lfp_norm(m1_chan,:);
                        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
                        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
                        tmp_pre_m1_lfp = eegfilt(tmp_pre_m1_lfp(1:data.streams.LFPs.fs*300),data.streams.LFPs.fs,4,8);
                        tmp_pre_m1_lfp = angle(hilbert(tmp_pre_m1_lfp));

                        tmp_post_m1_lfp = post_lfp_norm(m1_chan,:);
                        tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop(1):post_video_start_stop(2));
                        tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
                        tmp_post_m1_lfp = eegfilt(tmp_post_m1_lfp(1:data.streams.LFPs.fs*300),data.streams.LFPs.fs,4,8);
                        tmp_post_m1_lfp = angle(hilbert(tmp_post_m1_lfp));

                        for dls_chan = 1:16

                            if increase_channels(pair_count)==0
                                pair_count = pair_count + 1;
                                continue
                            end

                            pair_count = pair_count + 1;

                            tmp_pre_dls_lfp = pre_lfp_norm(dls_chan,:);
                            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
                            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_beh_state.nrem_interp==1);
                            tmp_pre_dls_lfp = eegfilt(tmp_pre_dls_lfp(1:data.streams.LFPs.fs*300),data.streams.LFPs.fs,4,8);
                            tmp_pre_dls_lfp = angle(hilbert(tmp_pre_dls_lfp));

                            tmp_post_dls_lfp = post_lfp_norm(dls_chan,:);
                            tmp_post_dls_lfp = tmp_post_dls_lfp(post_video_start_stop(1):post_video_start_stop(2));
                            tmp_post_dls_lfp = tmp_post_dls_lfp(post_beh_state.nrem_interp==1);
                            tmp_post_dls_lfp = eegfilt(tmp_post_dls_lfp(1:data.streams.LFPs.fs*300),data.streams.LFPs.fs,4,8);
                            tmp_post_dls_lfp = angle(hilbert(tmp_post_dls_lfp));

                            pair_phase_diff = [tmp_pre_dls_lfp-tmp_pre_m1_lfp tmp_post_dls_lfp-tmp_post_m1_lfp];
                            pair_phase_diff(pair_phase_diff>pi) = pair_phase_diff(pair_phase_diff>pi)-(2*pi);
                            pair_phase_diff(pair_phase_diff<-pi) = pair_phase_diff(pair_phase_diff<-pi)+(2*pi);

                            [~, i] = max(histc(pair_phase_diff,[-pi:pi/24:pi]));
                            tmp_x = [-pi:pi/24:pi];
                            all_pair_max_phase_diff = [all_pair_max_phase_diff tmp_x(i)];

                        end
                    end
                end

        %%% T201

            %%% FIND INCREASE CHANNELS

                load([data_path '\data\T201\sleep_activity\LFP_coherence.mat'])

                increase_channels = [];
                for pairs=1:512
                    if mean(mean([pre_LFP_coh{pairs}(16,48:113); post_LFP_coh{pairs}(16,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)])))
                        increase_channels = [increase_channels 1];
                    else
                        increase_channels = [increase_channels 0];
                    end
                end

            %%% FIND HIGH LFP COH PHASE DIFF

            day_blocks = {[1],[2 3],[4 5 6 7],[8 9 10],[11 12 13],[14 15 16 17],[18 19 20],[21 22 23 24 25 26],[27 28 29], ...
                [30 31 32],[33 34 35 36 37],[38 39 40],[41 42 43 44],[45 46 47 48],[49 50 51 52],[53 54 55 56],[57 58 59 60 61]};

            pre_sleep_blocks = {[1],[1],[2],[1],[1],[2],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1]};
            post_sleep_blocks = {[],[2],[4],[3],[3],[4],[3],[6],[3],[3],[5],[3],[4],[4],[4],[4],[5]};

                for day = 13:16

                    tmp_day_blocks = day_blocks{day};

                    %%% LOAD BEH STATE, LFP, AND START/STOP

                    load(['D:\Lemke_etal_2021\data\T201\sleep_rhythms\day_' num2str(day) '\block_' num2str(tmp_day_blocks(pre_sleep_blocks{day})) '\sleep_classification_rhythms.mat'])
                    pre_beh_state = beh_state;
                    if day == 13
                        data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(tmp_day_blocks(pre_sleep_blocks{day}))],'Type',4,'T2',7200);
                    else
                        data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(tmp_day_blocks(pre_sleep_blocks{day}))],'Type',4);
                    end
                    pre_lfp = data.streams.LFPs;
                    pre_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
                    pre_lfp = pre_lfp.data;

                    pre_lfp_norm = [];
                    for n = 1:size(pre_lfp,1)
                        tmp_lfp = pre_lfp(n,:);
                        tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                        b = find(abs(tmp_lfp)>6);
                        tmp_lfp = pre_lfp(n,:);
                        tmp_lfp(b)=NaN;
                        tmp_mean=nanmean(tmp_lfp);
                        tmp_sd=nanstd(tmp_lfp);
                        tmp_lfp(b)=tmp_mean;
                        pre_lfp_norm =  [pre_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
                    end
                    pre_lfp_norm(1:2:63,:) = pre_lfp_norm(1:2:63,:)-repmat(mean(pre_lfp_norm(1:2:63,:),1),[size(pre_lfp_norm(1:2:63,:),1) 1]);
                    pre_lfp_norm(34:2:64,:) = pre_lfp_norm(34:2:64,:)-repmat(mean(pre_lfp_norm(34:2:64,:),1),[size(pre_lfp_norm(34:2:64,:),1) 1]);

                    load(['D:\Lemke_etal_2021\data\T201\sleep_rhythms\day_' num2str(day) '\block_' num2str(tmp_day_blocks(post_sleep_blocks{day})) '\sleep_classification_rhythms.mat'])
                    post_beh_state = beh_state;
                    data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(tmp_day_blocks(post_sleep_blocks{day}))],'Type',4);
                    post_lfp = data.streams.LFPs;
                    post_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first') find(data.streams.Wave.data(1,:)>1.4,1,'last')];
                    post_lfp = post_lfp.data;

                    post_lfp_norm = [];
                    for n = 1:size(post_lfp,1)
                        tmp_lfp = post_lfp(n,:);
                        tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                        b = find(abs(tmp_lfp)>6);
                        tmp_lfp = post_lfp(n,:);
                        tmp_lfp(b)=NaN;
                        tmp_mean=nanmean(tmp_lfp);
                        tmp_sd=nanstd(tmp_lfp);
                        tmp_lfp(b)=tmp_mean;
                        post_lfp_norm =  [post_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
                    end
                    post_lfp_norm(1:2:63,:) = post_lfp_norm(1:2:63,:)-repmat(mean(post_lfp_norm(1:2:63,:),1),[size(post_lfp_norm(1:2:63,:),1) 1]);
                    post_lfp_norm(34:2:64,:) = post_lfp_norm(34:2:64,:)-repmat(mean(post_lfp_norm(34:2:64,:),1),[size(post_lfp_norm(34:2:64,:),1) 1]);

                    pair_count = 1;
                    for m1_chan = 1:2:63

                        tmp_pre_m1_lfp = pre_lfp_norm(m1_chan,:);
                        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
                        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
                        tmp_pre_m1_lfp = eegfilt(tmp_pre_m1_lfp(1:data.streams.LFPs.fs*300),data.streams.LFPs.fs,4,8);
                        tmp_pre_m1_lfp = angle(hilbert(tmp_pre_m1_lfp));

                        tmp_post_m1_lfp = post_lfp_norm(m1_chan,:);
                        tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop(1):post_video_start_stop(2));
                        tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
                        tmp_post_m1_lfp = eegfilt(tmp_post_m1_lfp(1:data.streams.LFPs.fs*300),data.streams.LFPs.fs,4,8);
                        tmp_post_m1_lfp = angle(hilbert(tmp_post_m1_lfp));

                        for dls_chan = 34:2:64

                            if increase_channels(pair_count)==0
                                pair_count = pair_count + 1;
                                continue
                            end

                            pair_count = pair_count + 1;

                            tmp_pre_dls_lfp = pre_lfp_norm(dls_chan,:);
                            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
                            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_beh_state.nrem_interp==1);
                            tmp_pre_dls_lfp = eegfilt(tmp_pre_dls_lfp(1:data.streams.LFPs.fs*300),data.streams.LFPs.fs,4,8);
                            tmp_pre_dls_lfp = angle(hilbert(tmp_pre_dls_lfp));

                            tmp_post_dls_lfp = post_lfp_norm(dls_chan,:);
                            tmp_post_dls_lfp = tmp_post_dls_lfp(post_video_start_stop(1):post_video_start_stop(2));
                            tmp_post_dls_lfp = tmp_post_dls_lfp(post_beh_state.nrem_interp==1);
                            tmp_post_dls_lfp = eegfilt(tmp_post_dls_lfp(1:data.streams.LFPs.fs*300),data.streams.LFPs.fs,4,8);
                            tmp_post_dls_lfp = angle(hilbert(tmp_post_dls_lfp));

                            pair_phase_diff = [tmp_pre_dls_lfp-tmp_pre_m1_lfp tmp_post_dls_lfp-tmp_post_m1_lfp];
                            pair_phase_diff(pair_phase_diff>pi) = pair_phase_diff(pair_phase_diff>pi)-(2*pi);
                            pair_phase_diff(pair_phase_diff<-pi) = pair_phase_diff(pair_phase_diff<-pi)+(2*pi);

                            [~, i] = max(histc(pair_phase_diff,[-pi:pi/24:pi]));
                            tmp_x = [-pi:pi/24:pi];
                            all_pair_max_phase_diff = [all_pair_max_phase_diff tmp_x(i)];

                        end
                    end
                end

        %%% SLDD2

            %%% FIND INCREASE CHANNELS

                load([data_path '\data\SLDD2\sleep_activity\LFP_coherence.mat'])

                increase_channels = [];
                for pairs=1:1024
                    if mean(mean([pre_LFP_coh{pairs}(11,48:113); post_LFP_coh{pairs}(11,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)])))
                        increase_channels = [increase_channels 1];
                    else
                        increase_channels = [increase_channels 0];
                    end
                end   

            %%% FIND HIGH LFP COH PHASE DIFF

            day_blocks = {{'SLDD2-191113-143329'}, ...
                            {'SLDD2-191114-121834'}, ...
                            {'SLDD2-191115-154323'}, ...
                            {'SLDD2-191116-145119','SLDD2-191116-170045','SLDD2-191116-193647'}, ...
                            {'SLDD2-191117-150115','SLDD2-191117-170253','SLDD2-191117-181931'}, ...
                            {'SLDD2-191118-154100','SLDD2-191118-180452','SLDD2-191118-192527','SLDD2-191118-205741'}, ...
                            {'SLDD2-191119-160043','SLDD2-191119-180437','SLDD2-191119-183709','SLDD2-191119-190631','SLDD2-191119-193444'}, ...
                            {'SLDD2-191120-140404','SLDD2-191120-160657','SLDD2-191120-170507'}, ...
                            {'SLDD2-191121-170233','SLDD2-191121-180559','SLDD2-191121-182818','SLDD2-191121-192200'}, ...
                            {'SLDD2-191122-152109','SLDD2-191122-172516','SLDD2-191122-184340','SLDD2-191122-204428'}, ...
                            {'SLDD2-191123-173250','SLDD2-191123-193545','SLDD2-191123-205347'}, ...
                            {'SLDD2-191124-152841','SLDD2-191124-173146','SLDD2-191124-183235'}, ...
                            {'SLDD2-191125-141910','SLDD2-191125-162007','SLDD2-191125-171539'}, ...
                            {'SLDD2-191126-173517','SLDD2-191126-191048','SLDD2-191126-195858'}};

            pre_sleep_blocks = {[],[],[],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1]};
            post_sleep_blocks = {[],[],[],[3],[3],[4],[5],[3],[4],[3],[3],[3],[3],[3]};

            pre_sleep_plx_blocks = {[],[],[],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1]};
            post_sleep_plx_blocks = {[],[],[],[3],[3],[4],[5],[3],[4],[3],[3],[3],[3],[3]};

            plx_files = {'sldd2-191113-143329-01', ...
                        'sldd2-191114-121834-01', ...
                        'sldd2-191115-154323-01', ...
                        'sldd2-191116-145119_mrg-01', ...
                        'sldd2-191117-150115_mrg-01', ...
                        'sldd2-191118-154100_mrg-01', ...
                        'sldd2-191119-160043_mrg-01', ...
                        'sldd2-191120-140404_mrg-01', ...
                        'sldd2-191121-170233_mrg-01', ...
                        'sldd2-191122-152109_mrg-01', ...
                        'sldd2-191123-173250_mrg-01', ...
                        'sldd2-191124-152841_mrg-01', ...
                        'sldd2-191125-141910_mrg-01', ...
                        'sldd2-191126-173517_mrg-01'};

            for day = 4:length(pre_sleep_blocks)

                tmp_sleep_blocks = day_blocks{day};

                %%% LOAD BEH STATE AND START/STOP

                    load(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\Lemke_etal_sleep_2\data\SLDD2\sleep_rhythms\day_' num2str(day) '\block_' day_blocks{day}{pre_sleep_blocks{day}} '\sleep_classification_rhythms.mat'])
                    pre_beh_state = beh_state;
                    wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\TDT_tanks\' day_blocks{day}{pre_sleep_blocks{day}}],'TYPE',[2]);
                    pre_lfp = [];
                    for chan = 1:64
                        data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\RS4_data\' day_blocks{day}{pre_sleep_blocks{day}}],'CHANNEL',chan);
                        pre_lfp = [pre_lfp; downsample(data.RSn1.data,20)];
                    end
                    pre_video_start_stop = [min(wave.epocs.PtC0.onset)*(data.RSn1.fs/20) max(wave.epocs.PtC0.onset)*(data.RSn1.fs/20)];

                    pre_lfp_norm = [];
                    for n = 1:size(pre_lfp,1)
                        tmp_lfp = double(pre_lfp(n,:));
                        tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                        b = find(abs(tmp_lfp)>6);
                        tmp_lfp = double(pre_lfp(n,:));
                        tmp_lfp(b)=NaN;
                        tmp_mean=nanmean(tmp_lfp);
                        tmp_sd=nanstd(tmp_lfp);
                        tmp_lfp(b)=tmp_mean;
                        pre_lfp_norm =  [pre_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
                    end
                    pre_lfp_norm([8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33],:) = pre_lfp_norm([8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33],:)-repmat(mean(pre_lfp_norm([8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33],:),1),[size(pre_lfp_norm([8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33],:),1) 1]);
                    pre_lfp_norm([2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63],:) = pre_lfp_norm([2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63],:)-repmat(mean(pre_lfp_norm([2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63],:),1),[size(pre_lfp_norm([2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63],:),1) 1]);

                    load(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\Lemke_etal_sleep_2\data\SLDD2\sleep_rhythms\day_' num2str(day) '\block_' day_blocks{day}{post_sleep_blocks{day}} '\sleep_classification_rhythms.mat'])
                    post_beh_state = beh_state;
                    wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\TDT_tanks\' day_blocks{day}{post_sleep_blocks{day}}],'TYPE',[2]);
                    post_lfp = [];
                    for chan = 1:64
                        data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\RS4_data\' day_blocks{day}{post_sleep_blocks{day}}],'CHANNEL',chan);
                        post_lfp = [post_lfp; downsample(data.RSn1.data,20)];
                    end
                    post_video_start_stop = [min(wave.epocs.PtC0.onset)*(data.RSn1.fs/20) max(wave.epocs.PtC0.onset)*(data.RSn1.fs/20)];

                    post_lfp_norm = [];
                    for n = 1:size(post_lfp,1)
                        tmp_lfp = double(post_lfp(n,:));
                        tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                        b = find(abs(tmp_lfp)>6);
                        tmp_lfp = double(post_lfp(n,:));
                        tmp_lfp(b)=NaN;
                        tmp_mean=nanmean(tmp_lfp);
                        tmp_sd=nanstd(tmp_lfp);
                        tmp_lfp(b)=tmp_mean;
                        post_lfp_norm =  [post_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
                    end
                    post_lfp_norm([8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33],:) = post_lfp_norm([8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33],:)-repmat(mean(post_lfp_norm([8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33],:),1),[size(post_lfp_norm([8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33],:),1) 1]);
                    post_lfp_norm([2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63],:) = post_lfp_norm([2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63],:)-repmat(mean(post_lfp_norm([2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63],:),1),[size(post_lfp_norm([2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63],:),1) 1]);

                    pair_count = 1;
                    for m1_chan = [8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33]

                        tmp_pre_m1_lfp = pre_lfp_norm(m1_chan,:);
                        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
                        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
                        tmp_pre_m1_lfp = eegfilt(tmp_pre_m1_lfp(1:data.RSn1.fs/20*300),data.RSn1.fs/20,4,8);
                        tmp_pre_m1_lfp = angle(hilbert(tmp_pre_m1_lfp));

                        tmp_post_m1_lfp = post_lfp_norm(m1_chan,:);
                        tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop(1):post_video_start_stop(2));
                        tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
                        tmp_post_m1_lfp = eegfilt(tmp_post_m1_lfp(1:data.RSn1.fs/20*300),data.RSn1.fs/20,4,8);
                        tmp_post_m1_lfp = angle(hilbert(tmp_post_m1_lfp));

                        for dls_chan = [2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63]

                            if increase_channels(pair_count)==0
                                pair_count = pair_count + 1;
                                continue
                            end

                            pair_count = pair_count + 1;

                            tmp_pre_dls_lfp = pre_lfp_norm(dls_chan,:);
                            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
                            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_beh_state.nrem_interp==1);
                            tmp_pre_dls_lfp = eegfilt(tmp_pre_dls_lfp(1:data.RSn1.fs/20*300),data.RSn1.fs/20,4,8);
                            tmp_pre_dls_lfp = angle(hilbert(tmp_pre_dls_lfp));

                            tmp_post_dls_lfp = post_lfp_norm(dls_chan,:);
                            tmp_post_dls_lfp = tmp_post_dls_lfp(post_video_start_stop(1):post_video_start_stop(2));
                            tmp_post_dls_lfp = tmp_post_dls_lfp(post_beh_state.nrem_interp==1);
                            tmp_post_dls_lfp = eegfilt(tmp_post_dls_lfp(1:data.RSn1.fs/20*300),data.RSn1.fs/20,4,8);
                            tmp_post_dls_lfp = angle(hilbert(tmp_post_dls_lfp));

                            pair_phase_diff = [tmp_pre_dls_lfp-tmp_pre_m1_lfp tmp_post_dls_lfp-tmp_post_m1_lfp];
                            pair_phase_diff(pair_phase_diff>pi) = pair_phase_diff(pair_phase_diff>pi)-(2*pi);
                            pair_phase_diff(pair_phase_diff<-pi) = pair_phase_diff(pair_phase_diff<-pi)+(2*pi);

                            [~, i] = max(histc(pair_phase_diff,[-pi:pi/24:pi]));
                            tmp_x = [-pi:pi/24:pi];
                            all_pair_max_phase_diff = [all_pair_max_phase_diff tmp_x(i)];

                        end
                    end
                end

        %%% T398

            %%% FIND INCREASE CHANNELS

                load([data_path '\data\T398\sleep_activity\LFP_coherence.mat'])

                increase_channels = [];
                for pairs=1:1406
                    if mean(mean([pre_LFP_coh{pairs}(7,48:113); post_LFP_coh{pairs}(7,48:113)])) > (mean(mean([pre_LFP_coh{pairs}(1,48:113); post_LFP_coh{pairs}(1,48:113)])))
                        increase_channels = [increase_channels 1];
                    else
                        increase_channels = [increase_channels 0];
                    end
                end   

            %%% FIND HIGH LFP COH PHASE DIFF

            sleep_blocks = {{'T398-190226-112203','T398-190226-160611'}, ...
                {'T398-190227-113133','T398-190227-154618'}, ...
                {'T398-190228-105336','T398-190228-141352'}, ...
                {'T398-190301-112846','T398-190301-151942'}, ...
                {'T398-190302-131102','T398-190302-162611'}, ...
                {'T398-190303-140804','T398-190303-171207'}, ...
                {'T398-190304-114807','T398-190304-145843'}, ...
                {'T398-190305-130625','T398-190305-161031'}, ...
                {'T398-190306-091252','T398-190306-124648'}};

            sleep_vids = {{'T398-(2019_2_26)-(11h22m)-SPONT.avi','T398-(2019_2_26)-(16h7m)-SPONT.avi'}, ...
                {'T398-(2019_2_27)-(11h32m)-SPONT.avi','T398-(2019_2_27)-(15h46m)-SPONT.avi'}, ...
                {'T398-(2019_2_28)-(10h54m)-SPONT.avi','T398-(2019_2_28)-(14h14m)-SPONT.avi'}, ...
                {'T398-(2019_3_1)-(11h29m)-SPONT.avi','T398-(2019_3_1)-(15h28m)-SPONT.avi'}, ...
                {'T398-(2019_3_2)-(13h11m)-SPONT.avi','T398-(2019_3_2)-(16h26m)-SPONT.avi'}, ...
                {'T398-(2019_3_3)-(14h8m)-SPONT','T398-(2019_3_3)-(17h12m)-SPONT'}, ...
                {'T398-(2019_3_4)-(11h48m)-SPONT.avi','T398-(2019_3_4)-(14h59m)-SPONT.avi'}, ...
                {'T398-(2019_3_5)-(13h7m)-SPONT.avi','T398-(2019_3_5)-(16h10m)-SPONT.avi'}, ...
                {'T398-(2019_3_6)-(9h13m)-SPONT.avi','T398-(2019_3_6)-(12h47m)-SPONT.avi'}};

            % BANKS
            bank_1_2 = {[5 6 7 8 9 21 39 40 41 42 43 46 51 52 53 54 56 57 58 59 60 61 100 102 105 107 108 110 111 113 114 115 116 119 120 122 123], ...
                [129 130 133 140 141 153 159 160 166 168 170 173 176 179 180 181 187 188 194 195 196 197 198 199 200 201 202 203 205 0206 207 211 212 214 223 224 231]};
            bank_2_1 = {[1 2 3 5 13 19 25 31 32 38 40 42 44 45 48 51 52 53 60 67 68 69 70 71 72 73 74 75 77 78 79 83 84 86 95 96 103], ...
                [135 136 137 139 140 149 167 168 169 170 171 179 180 181 182 184 185 186 187 188 189 228 233 235 236 238 239 241 242 243 244 247 248 249 250 251 255]};

            % BANK ORDER (M1 then DLS)
            bank_order = {[1 2], [1 2], [1 2], [2 1], [1 2], [2 1], [1 2], [1 2], [2 1]};

            for day = [5:9]
                for block = 1:length(sleep_blocks{day})

                    if bank_order{day}(1)==1
                        M1_bank = bank_1_2{1};
                        M1_offset = 1;
                        DLS_bank = bank_1_2{2};
                        DLS_offset = 129;
                    elseif bank_order{day}(1)==2
                        M1_bank = bank_2_1{2};
                        M1_offset = 129;
                        DLS_bank = bank_2_1{1};
                        DLS_offset = 1;
                    end

                    M1_mapping = {(M1_offset+[15 16 14 17 13 18 12 19]), ...%01
                        (M1_offset+[11 20 10 21 9 22 8 23]), ...            %02
                        [], ...                                             %03
                        [], ...                                             %04
                        (M1_offset+[63 32 62 33 61 34 60 35]), ...          %05
                        (M1_offset+[59 36 58 37 57 38 56 39]), ...          %06
                        [], ...                                             %07
                        [], ...                                             %08
                        (M1_offset+[79 80 78 81 77 82 76 83]), ...          %09
                        (M1_offset+[75 84 74 85 73 86 72 87]), ...          %10
                        [], ...                                             %11
                        [], ...                                             %12
                        (M1_offset+[127 96 126 97 125 98 124 99]), ...      %13
                        (M1_offset+[123 100 122 101 121 102 120 103]), ...  %14
                        [], ...                                             %15
                        [], ...                                             %16
                        (M1_offset+[111 112 110 113 109 114 108 115]), ...  %17
                        (M1_offset+[107 116 106 117 105 118 104 119]), ...  %18
                        [], ...                                             %19
                        [], ...                                             %20
                        (M1_offset+[95 64 94 65 93 66 92 67]), ...          %21
                        (M1_offset+[91 68 90 69 89 70 88 71]), ...          %22
                        [], ...                                             %23
                        [], ...                                             %24
                        (M1_offset+[47 48 46 49 45 50 44 51]), ...          %25
                        (M1_offset+[43 52 42 53 41 54 40 55]), ...          %26
                        [], ...                                             %27
                        [], ...                                             %28
                        (M1_offset+[31 0 30 1 29 2 28 3]), ...              %29
                        (M1_offset+[27 4 26 5 25 6 24 7]), ...              %30
                        [], ...                                             %31
                        []};                                                %32

                    DLS_mapping = {[], ...                                  %01
                        [], ...                                             %02
                        [], ...                                             %03
                        [], ...                                             %04
                        [], ...                                             %05
                        [], ...                                             %06
                        [], ...                                             %07
                        [], ...                                             %08
                        [], ...                                             %09
                        [], ...                                             %10
                        [], ...                                             %11
                        [], ...                                             %12
                        [], ...                                             %13
                        [], ...                                             %14
                        [], ...                                             %15
                        [], ...                                             %16
                        (DLS_offset+[15 16 14 17 13 18 12 19]), ...         %17
                        (DLS_offset+[11 20 10 21 9 22 8 23]), ...           %18
                        (DLS_offset+[7 24 6 25 5 26 4 27]), ...             %19
                        (DLS_offset+[3 28 2 29 1 30 0 31]), ...             %20
                        (DLS_offset+[63 32 62 33 61 34 60 35]), ...         %21
                        (DLS_offset+[59 36 58 37 57 38 56 39]), ...         %22
                        (DLS_offset+[55 40 54 41 53 42 52 43]), ...         %23
                        (DLS_offset+[51 44 50 45 49 46 48 47]), ...         %24
                        (DLS_offset+[79 80 78 81 77 82 76 83]), ...         %25
                        (DLS_offset+[75 84 74 85 73 86 72 87]), ...         %26
                        (DLS_offset+[71 88 70 89 69 90 68 91]), ...         %27
                        (DLS_offset+[67 92 66 93 65 94 64 95]), ...         %28
                        (DLS_offset+[127 96 126 97 125 96 124 95]), ...     %29
                        (DLS_offset+[123 100 122 101 121 102 120 103]), ... %30
                        (DLS_offset+[119 104 118 105 117 106 116 107]), ... %31
                        (DLS_offset+[115 108 114 109 113 110 112 111])};    %32

                    m1_lfp = [];
                    for shank = [1 2 5 6 9 10 13 14 17 18 21 22 25 26 29 30]
                        disp(shank);
                        [tmp_chans, ~] = intersect(M1_mapping{shank},M1_bank);
                        for electrode = 1:length(tmp_chans)
                            data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\RS4_data\M1-DLS_256Probe-190222-144047\' sleep_blocks{day}{block}],'CHANNEL',tmp_chans(electrode));
                            m1_lfp = [m1_lfp; downsample(data.RSn1.data,20)];
                        end
                    end
                    dls_lfp = [];
                    for shank = 17:32
                        disp(shank);
                        [tmp_chans, ~] = intersect(DLS_mapping{shank},DLS_bank);
                        for electrode = 1:length(tmp_chans)
                            data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\RS4_data\M1-DLS_256Probe-190222-144047\' sleep_blocks{day}{block}],'CHANNEL',tmp_chans(electrode));
                            dls_lfp = [dls_lfp; downsample(data.RSn1.data,20)];
                        end
                    end
                    lfp_samp_rate = data.RSn1.fs/20;

                    load(['D:\Lemke_etal_2021\data\T398\sleep_rhythms\day_' num2str(day) '\block_' sleep_blocks{day}{block} '\sleep_classification_rhythms.mat'])

                    if sum(beh_state.nrem_interp==1)<lfp_samp_rate*300
                        continue
                    end

                    if day==7 && block==2
                        video_start_stop = [10*(data.RSn1.fs/20) 7210*(data.RSn1.fs/20)];
                    else
                        wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\TDT_tanks\M1-DLS_256Probe-190222-144047\' sleep_blocks{day}{block}],'TYPE',[2]);
                        video_start_stop = [min(wave.epocs.PC1_.onset)*(data.RSn1.fs/20) max(wave.epocs.PC1_.onset)*(data.RSn1.fs/20)];
                        if day == 2 && block == 1
                            video_start_stop(2) = video_start_stop(2)+7200*(data.RSn1.fs/20);
                        end
                    end
                    beh_state.nrem_interp = beh_state.nrem_interp(1:length(video_start_stop(1):video_start_stop(2)));

                    m1_lfp_norm = [];
                    for n = 1:size(m1_lfp,1)
                        tmp_lfp = double(m1_lfp(n,:));
                        tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                        b = find(abs(tmp_lfp)>6);
                        tmp_lfp = double(m1_lfp(n,:));
                        tmp_lfp(b)=NaN;
                        tmp_mean=nanmean(tmp_lfp);
                        tmp_sd=nanstd(tmp_lfp);
                        tmp_lfp(b)=tmp_mean;
                        m1_lfp_norm =  [m1_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
                    end
                    m1_lfp_norm = m1_lfp_norm-repmat(mean(m1_lfp_norm),[size(m1_lfp_norm,1) 1]);

                    dls_lfp_norm = [];
                    for n = 1:size(dls_lfp,1)
                        tmp_lfp = double(dls_lfp(n,:));
                        tmp_lfp = ((tmp_lfp-mean(tmp_lfp))/std(tmp_lfp));
                        b = find(abs(tmp_lfp)>6);
                        tmp_lfp = double(dls_lfp(n,:));
                        tmp_lfp(b)=NaN;
                        tmp_mean=nanmean(tmp_lfp);
                        tmp_sd=nanstd(tmp_lfp);
                        tmp_lfp(b)=tmp_mean;
                        dls_lfp_norm =  [dls_lfp_norm; ((tmp_lfp-tmp_mean)'/tmp_sd)'];
                    end
                    dls_lfp_norm = dls_lfp_norm-repmat(mean(dls_lfp_norm),[size(dls_lfp_norm,1) 1]);

                    pair_count = 1;

                    for m1_chan = 1:size(m1_lfp_norm,1)

                        disp(['M1 chan ' num2str(m1_chan)])

                        tmp_m1_lfp = m1_lfp_norm(m1_chan,:);
                        if day ~= 6
                            tmp_m1_lfp = tmp_m1_lfp(video_start_stop(1):video_start_stop(2));
                        end
                        tmp_m1_lfp = tmp_m1_lfp(beh_state.nrem_interp==1);
                        tmp_m1_lfp = eegfilt(tmp_m1_lfp(1:lfp_samp_rate*300),data.RSn1.fs/20,4,8);
                        tmp_m1_lfp = angle(hilbert(tmp_m1_lfp));

                        for dls_chan = 1:size(dls_lfp_norm,1)

                            tmp_dls_lfp = dls_lfp_norm(dls_chan,:);
                            if day ~= 6
                                tmp_dls_lfp = tmp_dls_lfp(video_start_stop(1):video_start_stop(2));
                            end
                            tmp_dls_lfp = tmp_dls_lfp(beh_state.nrem_interp==1);
                            tmp_dls_lfp = eegfilt(tmp_dls_lfp(1:lfp_samp_rate*300),data.RSn1.fs/20,4,8);
                            tmp_dls_lfp = angle(hilbert(tmp_dls_lfp));

                            pair_phase_diff = [tmp_dls_lfp-tmp_m1_lfp];
                            pair_phase_diff(pair_phase_diff>pi) = pair_phase_diff(pair_phase_diff>pi)-(2*pi);
                            pair_phase_diff(pair_phase_diff<-pi) = pair_phase_diff(pair_phase_diff<-pi)+(2*pi);

                            [~, i] = max(histc(pair_phase_diff,[-pi:pi/24:pi]));
                            tmp_x = [-pi:pi/24:pi];
                            all_pair_max_phase_diff = [all_pair_max_phase_diff tmp_x(i)];

                        end
                    end
                end
            end

        figure;
        histogram(all_pair_max_phase_diff,[-pi:pi/12:pi],'DisplayStyle','Stairs','Normalization','Probability') 
        ylim([0 0.1])
                
%% sFIG 2 | WHAT FREQ LFP COHERENCE CORRELATES TO BEHAVIORAL IMPROVEMENT? 

%%% BEHAVIOR
all_animal_corr = [];

    %%% LC3        
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

        x_corr = [];
        y_corr = [];
        for n = 1:length(xvel)
            x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
            y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
        end
        animal_corr = mean([x_corr; y_corr]);
        all_animal_corr = [all_animal_corr zscore(animal_corr)];
    
    %%% LC4
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

        x_corr = [];
        y_corr = [];
        for n = 1:length(xvel)
            x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
            y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
        end
        animal_corr = mean([x_corr; y_corr]);
        all_animal_corr = [all_animal_corr zscore(animal_corr)];    
        
    %%% LC6
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
        
        x_corr = [];
        y_corr = [];
        for n = 1:length(xvel)
            x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
            y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
        end
        animal_corr = mean([x_corr; y_corr]);
        all_animal_corr = [all_animal_corr zscore(animal_corr)];
        
    %%% LC5
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
      
        x_corr = [];
        y_corr = [];
        for n = 1:length(xvel)
            x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
            y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
        end
        animal_corr = mean([x_corr; y_corr]);
        all_animal_corr = [all_animal_corr zscore(animal_corr)];
        
    %%% LC2
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
        
        x_corr = [];
        y_corr = [];
        for n = 1:length(xvel)
            x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
            y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
        end
        animal_corr = mean([x_corr; y_corr]);
        all_animal_corr = [all_animal_corr zscore(animal_corr)];
        
	%%% LC1
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
        
        x_corr = [];
        y_corr = [];
        for n = 1:length(xvel)
            x_corr = [x_corr corr(nanmean(xvel{n}(:,:))',nanmean(xvel{end}(:,:))')];
            y_corr = [y_corr corr(nanmean(yvel{n}(:,:))',nanmean(yvel{end}(:,:))')];
        end
        animal_corr = mean([x_corr; y_corr]);
        all_animal_corr = [all_animal_corr zscore(animal_corr(1:7))];
        
%%% SLEEP COHERENCE
all_animal_coh = [];

    %%% LC3  
    load([data_path '\LC3\sleep_activity\LFP_coherence.mat'])
    freq_bin = [1 33; 33 65; 65 97; 97 129; 129 161; 161 193; 193 225];
    
    tmp_day_coh = [];
    for day=1:8
        tmp_bin_coh = [];
        for bin = 1:7
            tmp_pre = [];
            tmp_post = [];
            for pairs = 1:size(pre_LFP_coh,2)
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,freq_bin(bin,1):freq_bin(bin,2)))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,freq_bin(bin,1):freq_bin(bin,2)))];
            end
            tmp_bin_coh = [tmp_bin_coh mean([tmp_pre tmp_post])];
        end
        tmp_day_coh = [tmp_day_coh; tmp_bin_coh];
    end
    for bin = 1:7
       tmp_day_coh(:,bin) = zscore(tmp_day_coh(:,bin));
    end
    all_animal_coh = [all_animal_coh; tmp_day_coh];
    
    %%% LC4
    load([data_path '\LC4\sleep_activity\LFP_coherence.mat'])
    freq_bin = [1 33; 33 65; 65 97; 97 129; 129 161; 161 193; 193 225];
    
    tmp_day_coh = [];
    for day=2:6
        tmp_bin_coh = [];
        for bin = 1:7
            tmp_pre = [];
            tmp_post = [];
            for pairs = 1:size(pre_LFP_coh,2)
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,freq_bin(bin,1):freq_bin(bin,2)))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,freq_bin(bin,1):freq_bin(bin,2)))];
            end
            tmp_bin_coh = [tmp_bin_coh mean([tmp_pre tmp_post])];
        end
        tmp_day_coh = [tmp_day_coh; tmp_bin_coh];
    end
    for bin = 1:7
       tmp_day_coh(:,bin) = zscore(tmp_day_coh(:,bin));
    end
    all_animal_coh = [all_animal_coh; tmp_day_coh];
        
    %%% LC6  
    load([data_path '\LC6\sleep_activity\LFP_coherence.mat'])
    freq_bin = [1 33; 33 65; 65 97; 97 129; 129 161; 161 193; 193 225];
    
    tmp_day_coh = [];
    for day=1:9
        tmp_bin_coh = [];
        for bin = 1:7
            tmp_pre = [];
            tmp_post = [];
            for pairs = 1:size(pre_LFP_coh,2)
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,freq_bin(bin,1):freq_bin(bin,2)))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,freq_bin(bin,1):freq_bin(bin,2)))];
            end
            tmp_bin_coh = [tmp_bin_coh mean([tmp_pre tmp_post])];
        end
        tmp_day_coh = [tmp_day_coh; tmp_bin_coh];
    end
    for bin = 1:7
       tmp_day_coh(:,bin) = zscore(tmp_day_coh(:,bin));
    end
    all_animal_coh = [all_animal_coh; tmp_day_coh];
    
    %%% LC5    
    load([data_path '\LC5\sleep_activity\LFP_coherence.mat'])
    freq_bin = [1 33; 33 65; 65 97; 97 129; 129 161; 161 193; 193 225];
    
    tmp_day_coh = [];
    for day=2:15
        tmp_bin_coh = [];
        for bin = 1:7
            tmp_pre = [];
            tmp_post = [];
            for pairs = 1:size(pre_LFP_coh,2)
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,freq_bin(bin,1):freq_bin(bin,2)))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,freq_bin(bin,1):freq_bin(bin,2)))];
            end
            tmp_bin_coh = [tmp_bin_coh mean([tmp_pre tmp_post])];
        end
        tmp_day_coh = [tmp_day_coh; tmp_bin_coh];
    end
    for bin = 1:7
       tmp_day_coh(:,bin) = zscore(tmp_day_coh(:,bin));
    end
    all_animal_coh = [all_animal_coh; tmp_day_coh];

    %%% LC2      
    load([data_path '\LC2\sleep_activity\LFP_coherence.mat'])
    freq_bin = [1 28;28 55;55 81;81 108;108 135;135 162;162 188];
    
    tmp_day_coh = [];
    for day=1:10
        tmp_bin_coh = [];
        for bin = 1:7
            tmp_pre = [];
            tmp_post = [];
            for pairs = 1:size(pre_LFP_coh,2)
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,freq_bin(bin,1):freq_bin(bin,2)))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,freq_bin(bin,1):freq_bin(bin,2)))];
            end
            tmp_bin_coh = [tmp_bin_coh mean([tmp_pre tmp_post])];
        end
        tmp_day_coh = [tmp_day_coh; tmp_bin_coh];
    end
    for bin = 1:7
       tmp_day_coh(:,bin) = zscore(tmp_day_coh(:,bin));
    end
    all_animal_coh = [all_animal_coh; tmp_day_coh];
    
    %%% LC1       
    load([data_path '\LC1\sleep_activity\LFP_coherence.mat'])
    freq_bin = [1 28;28 55;55 81;81 108;108 135;135 162;162 188];
    
    tmp_day_coh = [];
    for day=1:7
        tmp_bin_coh = [];
        for bin = 1:7
            tmp_pre = [];
            tmp_post = [];
            for pairs = 1:size(pre_LFP_coh,2)
                tmp_pre = [tmp_pre mean(pre_LFP_coh{pairs}(day,freq_bin(bin,1):freq_bin(bin,2)))];
                tmp_post = [tmp_post mean(post_LFP_coh{pairs}(day,freq_bin(bin,1):freq_bin(bin,2)))];
            end
            tmp_bin_coh = [tmp_bin_coh mean([tmp_pre tmp_post])];
        end
        tmp_day_coh = [tmp_day_coh; tmp_bin_coh];
    end
    for bin = 1:7
       tmp_day_coh(:,bin) = zscore(tmp_day_coh(:,bin));
    end
    all_animal_coh = [all_animal_coh; tmp_day_coh];

fig = figure;
    
    subplot(2,4,1); hold on;
        scatter(all_animal_coh(:,1),all_animal_corr',30,[0 0 0],'filled')
        [R,P] = corrcoef(all_animal_coh(:,1),all_animal_corr')
        [p] = polyfit(all_animal_coh(:,1),all_animal_corr',1)
        x=[-2.5:.1:2.5]
        xlim([-2.5 2.5])
        y = x*p(1)+p(2)
        plot(x,y,'color','r','LineStyle','--','LineWidth',3)
        title(['BEHAVIOR VS 1-3Hz LFP COH- R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
    
    subplot(2,4,2); hold on;
        scatter(all_animal_coh(:,2),all_animal_corr',30,[0 0 0],'filled')
        [R,P] = corrcoef(all_animal_coh(:,2),all_animal_corr')
        [p] = polyfit(all_animal_coh(:,2),all_animal_corr',1)
        x=[-2.5:.1:2.5]
        xlim([-2.5 2.5])
        y = x*p(1)+p(2)
        plot(x,y,'color','r','LineStyle','--','LineWidth',3)
        title(['BEHAVIOR VS 3-5Hz LFP COH- R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
    
    subplot(2,4,3); hold on;
        scatter(all_animal_coh(:,3),all_animal_corr',30,[0 0 0],'filled')
        [R,P] = corrcoef(all_animal_coh(:,3),all_animal_corr')
        [p] = polyfit(all_animal_coh(:,3),all_animal_corr',1)
        x=[-2.5:.1:2.5]
        xlim([-2.5 2.5])
        y = x*p(1)+p(2)
        plot(x,y,'color','r','LineStyle','--','LineWidth',3)
        title(['BEHAVIOR VS 5-7Hz LFP COH- R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
    
    subplot(2,4,4); hold on;
        scatter(all_animal_coh(:,4),all_animal_corr',30,[0 0 0],'filled')
        [R,P] = corrcoef(all_animal_coh(:,4),all_animal_corr')
        [p] = polyfit(all_animal_coh(:,4),all_animal_corr',1)
        x=[-2.5:.1:2.5]
        xlim([-2.5 2.5])
        y = x*p(1)+p(2)
        plot(x,y,'color','r','LineStyle','--','LineWidth',3)
        title(['BEHAVIOR VS 7-9Hz LFP COH- R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
    
    subplot(2,4,5); hold on;
        scatter(all_animal_coh(:,5),all_animal_corr',30,[0 0 0],'filled')
        [R,P] = corrcoef(all_animal_coh(:,5),all_animal_corr')
        [p] = polyfit(all_animal_coh(:,5),all_animal_corr',1)
        x=[-2.5:.1:2.5]
        xlim([-2.5 2.5])
        y = x*p(1)+p(2)
        plot(x,y,'color','r','LineStyle','--','LineWidth',3)
        title(['BEHAVIOR VS 9-11Hz LFP COH- R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
    
    subplot(2,4,6); hold on;
        scatter(all_animal_coh(:,6),all_animal_corr',30,[0 0 0],'filled')
        [R,P] = corrcoef(all_animal_coh(:,6),all_animal_corr')
        [p] = polyfit(all_animal_coh(:,6),all_animal_corr',1)
        x=[-2.5:.1:2.5]
        xlim([-2.5 2.5])
        y = x*p(1)+p(2)
        plot(x,y,'color','r','LineStyle','--','LineWidth',3)
        title(['BEHAVIOR VS 11-13Hz LFP COH- R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
    
    subplot(2,4,7); hold on;
        scatter(all_animal_coh(:,7),all_animal_corr',30,[0 0 0],'filled')
        [R,P] = corrcoef(all_animal_coh(:,7),all_animal_corr')
        [p] = polyfit(all_animal_coh(:,7),all_animal_corr',1)
        x=[-2.5:.1:2.5]
        xlim([-2.5 2.5])
        y = x*p(1)+p(2)
        plot(x,y,'color','r','LineStyle','--','LineWidth',3)
        title(['BEHAVIOR VS 13-15Hz LFP COH- R:' num2str(R(1,2)) ' - P:' num2str(P(1,2))]);
