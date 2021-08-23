%% T102

clc; clear all; close all;

day_blocks = {[2 3 7 8 9 11 12 13], ...
              [14 16 18 19 20 21], ...
              [22 24 25 26 27 29], ...
              [30 31 32 33 34 35 36 37 38], ...
              [39 43 44 45 47], ...
              [48 49 51 53 54 56], ...
              [57 58 59 60 61], ...
              [62 63 64 65 66]};

pre_sleep_blocks = {[2],[1],[1],[1],[1],[1],[1],[1]};
post_sleep_blocks = {[5],[4],[5],[4],[4],[5],[4],[4]};

plx_files = {'T102_blocks_2_3_7_8_9_11_12_13-01','T102_blocks_14_16_18_19_20_21-01','T102_blocks_22_24_25_26_27_29-01','T102_blocks_30_31_32_33_34_35_36_37_38-01','T102_blocks_39_43_44_45_47-01','T102_blocks_48_49_51_53_54_56-01','T102_blocks_57_58_59_60_61-01','T102_blocks_62_63_64_65_66-01'};
save_path = 'C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T102\cross_correlations';

for day = 1:8
    
    tmp_day_blocks = day_blocks{day};
    
    %%% LOAD BEH STATE AND START/STOP
    
        load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T102\sleep_activity\day_' num2str(day) '_block_1_sleep_activity.mat'])
        pre_beh_state = beh_state;
        pre_sleep_rhythms = sleep_rhythms;
        pre_lfp = mean_lfp;
        data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\Block-' num2str(tmp_day_blocks(pre_sleep_blocks{day}))],'Type',4);
        pre_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];

        load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T102\sleep_activity\day_' num2str(day) '_block_2_sleep_activity.mat'])
        post_beh_state = beh_state;
        post_sleep_rhythms = sleep_rhythms;
        post_lfp = mean_lfp;
        data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\Block-' num2str(tmp_day_blocks(post_sleep_blocks{day}))],'Type',4);
        post_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];

        clearvars data beh_state mean_lfp sleep_rhythms
    
    %%% LOAD SPIKING DATA
   
        plx_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\' plx_files{day} '.txt']);
        plx_spiking = table2array(plx_spiking);
        
        tmp_ts = cell(32,length(tmp_day_blocks));
        block_count = 1;
        for block = tmp_day_blocks
            tmp_num = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\Block-' num2str(block)],'Type',3,'NODATA',1);
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

        pre_m1_spiking = cell(16,max(plx_spiking(:,2)));
        post_m1_spiking = cell(16,max(plx_spiking(:,2)));
        count = 1;
        for m1_chan = 1:16
            for m1_unit = 1:max(plx_spiking(:,2))
                pre_m1_spiking{count,m1_unit} = tmp_ts{m1_chan,pre_sleep_blocks{day}}(tmp_sort{m1_chan,pre_sleep_blocks{day}}==m1_unit);
                post_m1_spiking{count,m1_unit} = tmp_ts{m1_chan,post_sleep_blocks{day}}(tmp_sort{m1_chan,post_sleep_blocks{day}}==m1_unit);
            end
            count = count + 1;
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
        
        [all_pairs, all_m1_spindle, all_dls_spindle, all_m1_SO, all_dls_SO, all_m1_delta, all_dls_delta] = CC_by_sleep_rhythm(pre_m1_spiking, pre_dls_spiking, post_m1_spiking, post_dls_spiking, pre_lfp, post_lfp, pre_video_start_stop, post_video_start_stop,pre_sleep_rhythms,post_sleep_rhythms);

   %%% SAVE

       if ~exist([save_path '\day_' num2str(day) '\'])
           mkdir([save_path '\day_' num2str(day) '\'])
       end
       cd([save_path '\day_' num2str(day) '\']);
       save('CC_by_sleep_rhythm.mat','all_pairs','-v7.3');

    pause(0.1);
   
end

%% T107

clc; clear all; close all;

day_blocks = {[1 2],[3 4 5 6 8],[9 10 11 12 13 14 15],[16 17 18 19 20 21 22 24 25],[26 27 30 31 36 39],[41 42 43 47 49 51 52]};
pre_sleep_blocks = {[1],[1],[1],[3],[1],[1]};
post_sleep_blocks = {[2],[4],[6],[8],[5],[4]};

plx_files = {'T107_blocks_1_2-01','T107_blocks_3_4_5_6_8-01','T107_blocks_9_10_11_12_13_14_15-01','T107_blocks_16_17_18_19_20_21_22_24_25-01','T107_blocks_26_27_30_31_36_39-01','T107_blocks_41_42_43_47_49_51_52-01'};
save_path = 'C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T107\cross_correlations';

for day = 1:6
    
    tmp_day_blocks = day_blocks{day};
    
    %%% LOAD BEH STATE AND START/STOP
    
        load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T107\sleep_activity\day_' num2str(day) '_block_1_sleep_activity.mat'])
        pre_beh_state = beh_state;
        pre_sleep_rhythms = sleep_rhythms;
        pre_lfp = mean_lfp;
        data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T107\neural_data\Block-' num2str(tmp_day_blocks(pre_sleep_blocks{day}))],'Type',4);
        pre_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];

        load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T107\sleep_activity\day_' num2str(day) '_block_2_sleep_activity.mat'])
        post_beh_state = beh_state;
        post_sleep_rhythms = sleep_rhythms;
        post_lfp = mean_lfp;
        data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T107\neural_data\Block-' num2str(tmp_day_blocks(post_sleep_blocks{day}))],'Type',4);
        post_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];

        clearvars data beh_state mean_lfp sleep_rhythms
    
    %%% LOAD SPIKING DATA
   
        plx_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\T107\neural_data\' plx_files{day} '.txt']);
        plx_spiking = table2array(plx_spiking);
        
        tmp_ts = cell(32,length(tmp_day_blocks));
        block_count = 1;
        for block = tmp_day_blocks
            tmp_num = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T107\neural_data\Block-' num2str(block)],'Type',3,'NODATA',1);
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

        pre_m1_spiking = cell(16,max(plx_spiking(:,2)));
        post_m1_spiking = cell(16,max(plx_spiking(:,2)));
        count = 1;
        for m1_chan = 1:16
            for m1_unit = 1:max(plx_spiking(:,2))
                pre_m1_spiking{count,m1_unit} = tmp_ts{m1_chan,pre_sleep_blocks{day}}(tmp_sort{m1_chan,pre_sleep_blocks{day}}==m1_unit);
                post_m1_spiking{count,m1_unit} = tmp_ts{m1_chan,post_sleep_blocks{day}}(tmp_sort{m1_chan,post_sleep_blocks{day}}==m1_unit);
            end
            count = count + 1;
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
        
        [all_pairs, all_m1_spindle, all_dls_spindle, all_m1_SO, all_dls_SO, all_m1_delta, all_dls_delta] = CC_by_sleep_rhythm(pre_m1_spiking, pre_dls_spiking, post_m1_spiking, post_dls_spiking, pre_lfp, post_lfp, pre_video_start_stop, post_video_start_stop,pre_sleep_rhythms,post_sleep_rhythms);

   %%% SAVE

       if ~exist([save_path '\day_' num2str(day) '\'])
           mkdir([save_path '\day_' num2str(day) '\'])
       end
       cd([save_path '\day_' num2str(day) '\']);
       save('CC_by_sleep_rhythm.mat','all_pairs','all_m1_spindle', 'all_dls_spindle', 'all_m1_SO', 'all_dls_SO', 'all_m1_delta', 'all_dls_delta','-v7.3');

    pause(0.1)
end

%% T200

clc; clear all; close all;

day_blocks = {[1],[2 7 8],[9 11 12 13 14],[15 16 17 18 19],[20 21 22 23 24],[25 26 27],[28 29 30],[32 33 34 35 36],[37 38 39 40 41],[42 43 44 45 46 47]};
pre_sleep_blocks = {[1],[1],[2],[1],[1],[1],[1],[1],[1],[1]};
post_sleep_blocks = {[],[3],[5],[5],[4],[3],[3],[5],[5],[6]};

plx_files = {'T200_blocks_1-01','T200_blocks_2_7_8-01','T200_blocks_9_11_12_13_14-01','T200_blocks_15_16_17_18_19-01','T200_blocks_20_21_22_23_24-01','T200_blocks_25_26_27-01','T200_blocks_28_29_30-01','T200_blocks_32_33_34_35_36-01','T200_blocks_37_38_39_40_41-01','T200_blocks_42_43_44_45_46_47-01'};
save_path = 'C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T200\cross_correlations';

for day = 2:10
    
    tmp_day_blocks = day_blocks{day};
    
    %%% LOAD BEH STATE AND START/STOP
    
        load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T200\sleep_activity\day_' num2str(day) '_block_1_sleep_activity.mat'])
        pre_beh_state = beh_state;
        pre_sleep_rhythms = sleep_rhythms;
        pre_lfp = mean_lfp;
        data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\Block-' num2str(tmp_day_blocks(pre_sleep_blocks{day}))],'Type',4);
        pre_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];
        
        load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T200\sleep_activity\day_' num2str(day) '_block_2_sleep_activity.mat'])
        post_beh_state = beh_state;
        post_sleep_rhythms = sleep_rhythms;
        post_lfp = mean_lfp;
        data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\Block-' num2str(tmp_day_blocks(post_sleep_blocks{day}))],'Type',4);
        post_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];

        clearvars data beh_state mean_lfp sleep_rhythms
    
    %%% LOAD SPIKING DATA
   
        plx_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\' plx_files{day} '.txt']);
        plx_spiking = table2array(plx_spiking);
        
        tmp_ts = cell(32,length(tmp_day_blocks));
        block_count = 1;
        for block = tmp_day_blocks
            tmp_num = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\Block-' num2str(block)],'Type',3,'NODATA',1);
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

        pre_m1_spiking = cell(16,max(plx_spiking(:,2)));
        post_m1_spiking = cell(16,max(plx_spiking(:,2)));
        count = 1;
        for m1_chan = 17:32
            for m1_unit = 1:max(plx_spiking(:,2))
                pre_m1_spiking{count,m1_unit} = tmp_ts{m1_chan,pre_sleep_blocks{day}}(tmp_sort{m1_chan,pre_sleep_blocks{day}}==m1_unit);
                post_m1_spiking{count,m1_unit} = tmp_ts{m1_chan,post_sleep_blocks{day}}(tmp_sort{m1_chan,post_sleep_blocks{day}}==m1_unit);
            end
            count = count + 1;
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
        
    [all_pairs, all_m1_spindle, all_dls_spindle, all_m1_SO, all_dls_SO, all_m1_delta, all_dls_delta] = CC_by_sleep_rhythm(pre_m1_spiking, pre_dls_spiking, post_m1_spiking, post_dls_spiking, pre_lfp, post_lfp, pre_video_start_stop, post_video_start_stop,pre_sleep_rhythms,post_sleep_rhythms);

    %%% SAVE

       if ~exist([save_path '\day_' num2str(day) '\'])
           mkdir([save_path '\day_' num2str(day) '\'])
       end
       cd([save_path '\day_' num2str(day) '\']);
       save('CC_by_sleep_rhythm.mat','all_pairs','all_m1_spindle', 'all_dls_spindle', 'all_m1_SO', 'all_dls_SO', 'all_m1_delta', 'all_dls_delta','-v7.3');

    pause(0.1)
   
end

%% T201

clc; clear all; close all;

day_blocks = {[1],[2 3],[4 5 6 7],[8 9 10],[11 12 13],[14 15 16 17],[18 19 20],[21 22 23 24 25 26],[27 28 29], ...
    [30 31 32],[33 34 35 36 37],[38 39 40],[41 42 43 44],[45 46 47 48],[49 50 51 52],[53 54 55 56],[57 58 59 60 61]};

pre_sleep_blocks = {[1],[1],[2],[1],[1],[2],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1]};
post_sleep_blocks = {[],[2],[4],[3],[3],[4],[3],[6],[3],[3],[5],[3],[4],[4],[4],[4],[5]};

plx_files = {'T201_blocks_1-01','T201_blocks_2_3-01','T201_blocks_4_5_6_7-01','T201_blocks_8_9_10-01','T201_blocks_11_12_13-01', ...
    'T201_blocks_14_15_16_17-01','T201_blocks_18_19_20-01','T201_blocks_21_22_23_24_25_26-01','T201_blocks_27_28_29-01', ...
    'T201_blocks_30_31_32-01','T201_blocks_33_34_35_36_37-01','T201_blocks_38_39_40-01','T201_blocks_41_42_43_44-01','T201_blocks_45_46_47_48-01', ...
    'T201_blocks_49_50_51_52-01','T201_blocks_53_54_55_56-01','T201_blocks_57_58_59_60_61-01'};

save_path = 'C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T201\cross_correlations';

for day = 2
    
    tmp_day_blocks = day_blocks{day};
    
    %%% LOAD BEH STATE AND START/STOP
    
        load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T201\sleep_activity\day_' num2str(day) '_block_1_sleep_activity.mat'])
        pre_beh_state = beh_state;
        pre_sleep_rhythms = sleep_rhythms;
        pre_lfp = mean_lfp;
        if day == 13
            data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(tmp_day_blocks(pre_sleep_blocks{day}))],'Type',4,'T2',7200);
            pre_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];
        else
            data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(tmp_day_blocks(pre_sleep_blocks{day}))],'Type',4);
            pre_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];
        end
        
        load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T201\sleep_activity\day_' num2str(day) '_block_2_sleep_activity.mat'])
        post_beh_state = beh_state;
        post_sleep_rhythms = sleep_rhythms;
        post_lfp = mean_lfp;
        data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(tmp_day_blocks(post_sleep_blocks{day}))],'Type',4);
        post_video_start_stop = [find(data.streams.Wave.data(1,:)>1.4,1,'first')/data.streams.Wave.fs find(data.streams.Wave.data(1,:)>1.4,1,'last')/data.streams.Wave.fs];

        clearvars data beh_state mean_lfp sleep_rhythms
    
    %%% LOAD SPIKING DATA
   
        plx_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\' plx_files{day} '.txt']);
        plx_spiking = table2array(plx_spiking);
        
        tmp_ts = cell(64,length(tmp_day_blocks));
        block_count = 1;
        for block = tmp_day_blocks
            if block == 41
                tmp_num = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(block)],'Type',3,'T2',7200);
            else
                tmp_num = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(block)],'Type',3,'NODATA',1);
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

        pre_m1_spiking = cell(32,max(plx_spiking(:,2)));
        post_m1_spiking = cell(32,max(plx_spiking(:,2)));
        count = 1;
        for m1_chan = 1:2:63
            for m1_unit = 1:max(plx_spiking(:,2))
                pre_m1_spiking{count,m1_unit} = tmp_ts{m1_chan,pre_sleep_blocks{day}}(tmp_sort{m1_chan,pre_sleep_blocks{day}}==m1_unit);
                post_m1_spiking{count,m1_unit} = tmp_ts{m1_chan,post_sleep_blocks{day}}(tmp_sort{m1_chan,post_sleep_blocks{day}}==m1_unit);
            end
            count = count + 1;
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
        
        [all_pairs, all_m1_spindle, all_dls_spindle, all_m1_SO, all_dls_SO, all_m1_delta, all_dls_delta] = CC_by_sleep_rhythm(pre_m1_spiking, pre_dls_spiking, post_m1_spiking, post_dls_spiking, pre_lfp, post_lfp, pre_video_start_stop, post_video_start_stop,pre_sleep_rhythms,post_sleep_rhythms);

    %%% SAVE

       if ~exist([save_path '\day_' num2str(day) '\'])
           mkdir([save_path '\day_' num2str(day) '\'])
       end
       cd([save_path '\day_' num2str(day) '\']);
       save('CC_by_sleep_rhythm.mat','all_pairs','all_m1_spindle', 'all_dls_spindle', 'all_m1_SO', 'all_dls_SO', 'all_m1_delta', 'all_dls_delta','-v7.3');

end

%% SLDD2

clc; clear all; close all;

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

save_path = 'C:\Users\Stefan\Desktop\Lemke_etal_2021\data\SLDD2\cross_correlations';

for day = [4:13]
    
    tmp_sleep_blocks = day_blocks{day};
    
    %%% LOAD BEH STATE AND START/STOP
    
        load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\SLDD2\sleep_activity\day_' num2str(day) '_block_1_sleep_activity.mat'])
        pre_beh_state = beh_state;
        pre_sleep_rhythms = sleep_rhythms;
        pre_lfp = mean_lfp;
        wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\TDT_tanks\' day_blocks{day}{pre_sleep_blocks{day}}],'TYPE',[2]);
        pre_video_start_stop = [min(wave.epocs.PtC0.onset) max(wave.epocs.PtC0.onset)];

        load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\SLDD2\sleep_activity\day_' num2str(day) '_block_2_sleep_activity.mat'])
        post_beh_state = beh_state;
        post_sleep_rhythms = sleep_rhythms;
        post_lfp = mean_lfp;
        wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\TDT_tanks\' day_blocks{day}{post_sleep_blocks{day}}],'TYPE',[2]);
        post_video_start_stop = [min(wave.epocs.PtC0.onset) max(wave.epocs.PtC0.onset)];

        clearvars wave beh_state mean_lfp sleep_rhythms
    
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
        tmp_width = cell(64,length(tmp_blocks_cell));
        tmp_peak2valley = cell(64,length(tmp_blocks_cell));
        for chan = 1:64
            chan_sorts = plx_spiking(plx_spiking(:,1)==chan,:);
            for block = 1:length(tmp_blocks_cell)
                tmp_sort{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),2);
                tmp_width{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),6);
                tmp_peak2valley{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),4);
            end
        end
        
        pre_m1_spiking = cell(32,max(plx_spiking(:,2)));
        post_m1_spiking = cell(32,max(plx_spiking(:,2)));
        count = 1;
        for m1_chan = [8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33]
            for m1_unit = 1:max(plx_spiking(:,2))
                pre_m1_spiking{count,m1_unit} = tmp_ts{m1_chan,pre_sleep_plx_blocks{day}}(tmp_sort{m1_chan,pre_sleep_plx_blocks{day}}==m1_unit);
                post_m1_spiking{count,m1_unit} = tmp_ts{m1_chan,post_sleep_plx_blocks{day}}(tmp_sort{m1_chan,post_sleep_plx_blocks{day}}==m1_unit);
            end
            count = count + 1;
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
        
    [all_pairs, all_m1_spindle, all_dls_spindle, all_m1_SO, all_dls_SO, all_m1_delta, all_dls_delta] = CC_by_sleep_rhythm(pre_m1_spiking, pre_dls_spiking, post_m1_spiking, post_dls_spiking, pre_lfp, post_lfp, pre_video_start_stop, post_video_start_stop,pre_sleep_rhythms,post_sleep_rhythms);

    %%% SAVE

       if ~exist([save_path '\day_' num2str(day) '\'])
           mkdir([save_path '\day_' num2str(day) '\'])
       end
       cd([save_path '\day_' num2str(day) '\']);
       save('CC_by_sleep_rhythm.mat','all_pairs','all_m1_spindle', 'all_dls_spindle', 'all_m1_SO', 'all_dls_SO', 'all_m1_delta', 'all_dls_delta','-v7.3');

    pause(0.1)
   
end

%% T398

clc; clear all; close all;

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

save_path = 'C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T398\cross_correlations';

for day = [1:9]

    %%% LOAD BEH STATE AND START/STOP

        load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T398\sleep_activity\day_' num2str(day) '_block_1_sleep_activity.mat'])
        pre_beh_state = beh_state;
        pre_sleep_rhythms = sleep_rhythms;
        pre_lfp = mean_lfp;
        data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\TDT_tanks\M1-DLS_256Probe-190222-144047\' pre_sleep{day}]);
        pre_video_start_stop = [min(data.epocs.PC1_.onset) max(data.epocs.PC1_.onset)];
        if day == 2
            pre_video_start_stop(2) = pre_video_start_stop(2)+7200;
        end
        
        load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T398\sleep_activity\day_' num2str(day) '_block_2_sleep_activity.mat'])
        post_beh_state = beh_state;
        post_sleep_rhythms = sleep_rhythms;
        post_lfp = mean_lfp;
        data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\TDT_tanks\M1-DLS_256Probe-190222-144047\' post_sleep{day}]);
        if day == 7
            post_video_start_stop = [10 7210];
        else
            post_video_start_stop = [min(data.epocs.PC1_.onset) max(data.epocs.PC1_.onset)];
        end
        clearvars data beh_state mean_lfp sleep_rhythms
    
    %%% LOAD SPIKING DATA

        opts = delimitedTextImportOptions("NumVariables", 2);
        opts.DataLines = [2, Inf];
        opts.Delimiter = "\t";
        opts.VariableNames = ["cluster_id", "group"];
        opts.VariableTypes = ["double", "categorical"];
        opts = setvaropts(opts, 2, "EmptyFieldRule", "auto");
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        
        cd([spike_path{day} '\M1']);
        m1_ts = readNPY('spike_times.npy');
        m1_id = readNPY('spike_clusters.npy');
        m1_unit_class = readtable('cluster_group.tsv', opts);
        m1_clust_id = zeros(1,max(m1_unit_class{:,1}));
        m1_clust_id_idx = m1_unit_class{:,1}'+1;
        idx_count = 1;
        for clust_id = m1_clust_id_idx
            tmp_qual = cellstr(m1_unit_class{idx_count,2});
            if strcmp(tmp_qual{1,1},'noise')
                m1_clust_id(clust_id) = 0;
            else
                m1_clust_id(clust_id) = 1;
            end
            idx_count = idx_count + 1;
        end
        m1_clust_id = find(m1_clust_id)-1;

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

        max(m1_ts)
        max(dls_ts)
        pre_sleep_length + reach_length + post_sleep_length
        %%% CHECK UNIT %%%
        
        pre_m1_spiking = cell(length(m1_clust_id),1);
        post_m1_spiking = cell(length(m1_clust_id),1);
        m1_count = 1;
        for m1_n = m1_clust_id
            tmp_spikes = m1_ts(m1_id==m1_n);
            pre_m1_spiking{m1_count,1} = double(tmp_spikes(tmp_spikes<pre_sleep_length))/samp_freq;
            post_m1_spiking{m1_count,1} = double(tmp_spikes(tmp_spikes>(pre_sleep_length+reach_length))-(pre_sleep_length+reach_length))/samp_freq;
            m1_count = m1_count+1;
        end
        
        pre_dls_spiking = cell(length(dls_clust_id),1);
        post_dls_spiking = cell(length(dls_clust_id),1);
        dls_count = 1;
        for dls_n = dls_clust_id
            tmp_spikes = dls_ts(dls_id==dls_n);
            pre_dls_spiking{dls_count,1} = double(tmp_spikes(tmp_spikes<pre_sleep_length))/samp_freq;
            post_dls_spiking{dls_count,1} = double(tmp_spikes(tmp_spikes>(pre_sleep_length+reach_length))-(pre_sleep_length+reach_length))/samp_freq;
            dls_count = dls_count+1;
        end
        
    [all_pairs, all_m1_spindle, all_dls_spindle, all_m1_SO, all_dls_SO, all_m1_delta, all_dls_delta] = CC_by_sleep_rhythm(pre_m1_spiking, pre_dls_spiking, post_m1_spiking, post_dls_spiking, pre_lfp, post_lfp, pre_video_start_stop, post_video_start_stop,pre_sleep_rhythms,post_sleep_rhythms);

    %%% SAVE

       if ~exist([save_path '\day_' num2str(day) '\'])
           mkdir([save_path '\day_' num2str(day) '\'])
       end
       cd([save_path '\day_' num2str(day) '\']);
       save('CC_by_sleep_rhythm.mat','all_pairs','all_m1_spindle', 'all_dls_spindle', 'all_m1_SO', 'all_dls_SO', 'all_m1_delta', 'all_dls_delta','-v7.3');

    pause(0.1)
   
end
            