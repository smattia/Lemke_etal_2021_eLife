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
reach_blocks = {[3 4 6 7], [2 3 5],[2 3 6],[2 3 5 6 7 9],[2 3 5],[3 4 6],[2 3],[2 3 5]};

plx_files = {'T102_blocks_2_3_7_8_9_11_12_13-01', ...
             'T102_blocks_14_16_18_19_20_21-01', ...
             'T102_blocks_22_24_25_26_27_29-01', ...
             'T102_blocks_30_31_32_33_34_35_36_37_38-01', ...
             'T102_blocks_39_43_44_45_47-01', ...
             'T102_blocks_48_49_51_53_54_56-01', ...
             'T102_blocks_57_58_59_60_61-01', ...
             'T102_blocks_62_63_64_65_66-01'};
         
day_dirs = {['L:\videos\DLC\T102_kin3\day_1\'], ...
            ['L:\videos\DLC\T102_kin3\day_2\'], ...
            ['L:\videos\DLC\T102_kin3\day_3\'], ...
            ['L:\videos\DLC\T102_kin3\day_4\'], ...
            ['L:\videos\DLC\T102_kin3\day_5\'], ...
            ['L:\videos\DLC\T102_kin3\day_6\'], ...
            ['L:\videos\DLC\T102_kin3\day_7\'], ...
            ['L:\videos\DLC\T102_kin3\day_8\']};

time_each_block = {{'13h26m','13h42m','16h19m','16h32m'}, ...
                   {'13h44m','14h21m','16h38m'}, ...
                   {'12h50m','13h4m','16h4m'}, ...
                   {'12h2m','12h15m','15h21m','15h34m','15h45m','18h6m'}, ...
                   {'12h43m','13h4m','16h25m'}, ...
                   {'13h56m','14h16m','16h58m'}, ...
                   {'12h14m','12h29m','12h42m'}, ...
                   {'11h52m','12h6m','12h19m'}};
        
save_path = '\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T102\reach_mod';

for day = 1:8
    
    tmp_day_blocks = day_blocks{day};
    
    %%% LOAD SPIKING DATA
   
        plx_spiking = readtable(['\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T102\T102\' plx_files{day} '.txt']);
        plx_spiking = table2array(plx_spiking);
        
        tmp_ts = cell(32,length(tmp_day_blocks));
        block_count = 1;
        for block = tmp_day_blocks
            tmp_num = TDTbin2mat(['\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T102\T102\Block-' num2str(block)],'Type',3,'NODATA',1);
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

        m1_spiking = cell(length(reach_blocks{day}),16,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for m1_chan = 1:16
                for m1_unit = 1:max(plx_spiking(:,2))
                    m1_spiking{blocks,count,m1_unit} = tmp_ts{m1_chan,reach_blocks{day}(blocks)}(tmp_sort{m1_chan,reach_blocks{day}(blocks)}==m1_unit);
                end
                count = count + 1;
            end
        end
        
        dls_spiking = cell(length(reach_blocks{day}),16,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for dls_chan = 17:32
                for dls_unit = 1:max(plx_spiking(:,2))
                    dls_spiking{blocks,dls_chan,dls_unit} = tmp_ts{dls_chan,reach_blocks{day}(blocks)}(tmp_sort{dls_chan,reach_blocks{day}(blocks)}==dls_unit);
                end
                count = count + 1;
            end
        end
        
	%%% LOAD BEHAVIORAL MARKERS
    
    beh_markers = cell(length(reach_blocks{day}),1);
    
    for blocks = 1:length(reach_blocks{day})
                
        %%% GET TDT TRIAL START TIMES
        clearvars final_pulses pellet_drops start_times
        wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T102\T102\Block-' num2str(day_blocks{day}(reach_blocks{day}(blocks)))],'TYPE',4);
        fs = wave.streams.Wave.fs;
        wave = wave.streams.Wave.data(1,:);
        thr = max(wave)*0.4;
        temp=wave>thr;
        pulses = find(temp==1);
        pulse_diff = diff(pulses);
        start_pulses = find(pulse_diff>2);
        start_pulses = start_pulses+1;
        start_pulses = [1 start_pulses];
        for n=1:length(start_pulses)
            final_pulses(n) = pulses(start_pulses(n));
        end
        diff_final_pulses = diff(final_pulses);
        num_pulses = length(diff_final_pulses);
        diff_final_pulses = [301 diff_final_pulses 301 301];
        p_count= 1;
        s_count = 1;
        for n=1:num_pulses+1
            % pellet drops
            if (diff_final_pulses(n) > 200 && diff_final_pulses(n+1) >300)
                pellet_drops(p_count) = final_pulses(n);
                p_count = p_count + 1;
            end
            % start times
            if (diff_final_pulses(n) > 300 && diff_final_pulses(n+1) < 20 && diff_final_pulses(n+2)>200)
                start_times(s_count) = final_pulses(n+1);
                s_count = s_count + 1;
            end
        end
        
        %%% GET KINEMATICS
        dlc_files = dir(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T102\reach_videos\day_' num2str(day) '\trajectories\*' time_each_block{day}{blocks} '*corrected*']);
        [~, i] = sort_nat({dlc_files(:).name});
        dlc_files = dlc_files(i);
        pellet_position = [];
        for file = 1:length(dlc_files)
            trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T102\reach_videos\day_' num2str(day) '\trajectories\' dlc_files(file).name]);
            tmp_pellet_x = trial_kin(1:25,5);
            tmp_pellet_y = trial_kin(1:25,6);
            tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
            tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
            pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
            pellet_position = [pellet_position; pellet_location];
        end
        pellet_position = nanmean(pellet_position);
    
        mat_files = dir(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T102\reach_videos\day_' num2str(day) '\mat_videos\*' time_each_block{day}{blocks} '*']);
        [~, i] = sort_nat({mat_files(:).name});
        mat_files = mat_files(i);
        
        tmp_beh_markers = [];
        tmp_dlc_frame = [];
        tmp_wave_time = []; 
        for file = 1:length(mat_files)
            trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T102\reach_videos\day_' num2str(day) '\trajectories\' dlc_files(file).name]);
            frame_times = load(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T102\reach_videos\day_' num2str(day) '\mat_videos\' mat_files(file).name]);
            frame_times = frame_times.t_lat;
            frame_times = frame_times*fs + start_times(file);
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
            if (i<31 || i>size(trial_kin,1)-31)
                tmp_beh_markers = [tmp_beh_markers nan];
            else
                tmp_beh_markers = [tmp_beh_markers frame_times(i)];
            end
            tmp_dlc_frame = [tmp_dlc_frame size(trial_kin,1)];
            tmp_wave_time = [tmp_wave_time pellet_drops(file)-start_times(file)];
        end
        beh_markers{blocks} = tmp_beh_markers/fs;
        
        figure 
            subplot(1,3,[1 2])
                hold on
                plot(wave);
                scatter(start_times,ones(1,length(start_times)))
                scatter(pellet_drops,ones(1,length(pellet_drops)))
                scatter(tmp_beh_markers,2*ones(1,length(tmp_beh_markers)))
                title(['Day ' num2str(day) ' Block ' num2str(blocks) ' - vids: ' num2str(length(mat_files)) ' - dlc files: ' num2str(length(dlc_files)) ' - wave starts: ' num2str(length(start_times))])
            subplot(1,3,3)
                hold on
                plot(zscore(tmp_dlc_frame));
                plot(zscore(tmp_wave_time));
                title('DLC frames vs. wave length');
    end

    [m1_reach, dls_reach] = reach_modulation(m1_spiking,dls_spiking,beh_markers);
    
    %%% SAVE
    if ~exist([save_path '\day_' num2str(day) '\'])
        mkdir([save_path '\day_' num2str(day) '\'])
    end
    cd([save_path '\day_' num2str(day) '\']);
    save('reach_mod.mat','m1_reach','dls_reach','-v7.3');
   
end

%% T107

clc; clear all; close all;

day_blocks = {[1 2], ...
              [3 4 5 6 8], ...
              [9 10 11 12 13 14 15], ...
              [16 17 18 19 20 21 22 24 25], ...
              [26 27 30 31 36 39], ...
              [41 42 43 47 49 51 52]};
reach_blocks = {[], [2 5],[2 4 5 7],[4 5 6 7 9],[3 4 6],[2 5 6 7]};

plx_files = {'T107_blocks_1_2-01', ...
             'T107_blocks_3_4_5_6_8-01', ...
             'T107_blocks_9_10_11_12_13_14_15-01', ...
             'T107_blocks_16_17_18_19_20_21_22_24_25-01', ...
             'T107_blocks_26_27_30_31_36_39-01', ...
             'T107_blocks_41_42_43_47_49_51_52-01'};

day_dirs = {[], ...
    ['L:\videos\DLC\T107_kin\day_2\'], ...
    ['L:\videos\DLC\T107_kin\day_3\'], ...
    ['L:\videos\DLC\T107_kin\day_4\'], ...
    ['L:\videos\DLC\T107_kin\day_5\'], ...
    ['L:\videos\DLC\T107_kin\day_6\']};

time_each_block = {{}, ...
                   {'10h43m','15h21m'}, ...
                   {'11h3m','14h55m','15h15m','15h30m'}, ...
                   {'12h10m','12h25m','12h42m','12h59m','15h11m'}, ...
                   {'14h12m','14h27m','16h56m'}, ...
                   {'13h33m','16h11m','16h23m','16h33m'}};
               
save_path = '\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T107\reach_mod';

for day = 2:6
    
    tmp_day_blocks = day_blocks{day};
    
    %%% LOAD SPIKING DATA
   
        plx_spiking = readtable(['\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T107\T107\' plx_files{day} '.txt']);
        plx_spiking = table2array(plx_spiking);
        
        tmp_ts = cell(32,length(tmp_day_blocks));
        block_count = 1;
        for block = tmp_day_blocks
            tmp_num = TDTbin2mat(['\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T107\T107\Block-' num2str(block)],'Type',3,'NODATA',1);
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

        m1_spiking = cell(length(reach_blocks{day}),16,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for m1_chan = 1:16
                for m1_unit = 1:max(plx_spiking(:,2))
                    m1_spiking{blocks,count,m1_unit} = tmp_ts{m1_chan,reach_blocks{day}(blocks)}(tmp_sort{m1_chan,reach_blocks{day}(blocks)}==m1_unit);
                end
                count = count + 1;
            end
        end
        
        dls_spiking = cell(length(reach_blocks{day}),16,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for dls_chan = 17:32
                for dls_unit = 1:max(plx_spiking(:,2))
                    dls_spiking{blocks,dls_chan,dls_unit} = tmp_ts{dls_chan,reach_blocks{day}(blocks)}(tmp_sort{dls_chan,reach_blocks{day}(blocks)}==dls_unit);
                end
                count = count + 1;
            end
        end
        
	%%% LOAD BEHAVIORAL MARKERS
    
    beh_markers = cell(length(reach_blocks{day}),1);
    
    for blocks = 1:length(reach_blocks{day})
                
        %%% GET TDT TRIAL START TIMES
        clearvars final_pulses pellet_drops start_times
        wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T107\T107\Block-' num2str(day_blocks{day}(reach_blocks{day}(blocks)))],'TYPE',4);
        fs = wave.streams.Wave.fs;
        wave = wave.streams.Wave.data(1,:);
        thr = max(wave)*0.4;
        temp=wave>thr;
        pulses = find(temp==1);
        pulse_diff = diff(pulses);
        start_pulses = find(pulse_diff>2);
        start_pulses = start_pulses+1;
        start_pulses = [1 start_pulses];
        for n=1:length(start_pulses)
            final_pulses(n) = pulses(start_pulses(n));
        end
        diff_final_pulses = diff(final_pulses);
        num_pulses = length(diff_final_pulses);
        diff_final_pulses = [301 diff_final_pulses 301 301];
        p_count= 1;
        s_count = 1;
        for n=1:num_pulses+1
            % pellet drops
            if (diff_final_pulses(n) > 200 && diff_final_pulses(n+1) >300)
                pellet_drops(p_count) = final_pulses(n);
                p_count = p_count + 1;
            end
            % start times
            if (diff_final_pulses(n) > 300 && diff_final_pulses(n+1) < 20 && diff_final_pulses(n+2)>200)
                start_times(s_count) = final_pulses(n+1);
                s_count = s_count + 1;
            end
        end
        
        %%% GET KINEMATICS
        dlc_files = dir(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T107\reach_videos\day_' num2str(day) '\trajectories\*' time_each_block{day}{blocks} '*corrected*']);
        [~, i] = sort_nat({dlc_files(:).name});
        dlc_files = dlc_files(i);
        pellet_position = [];
        for file = 1:length(dlc_files)
            trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T107\reach_videos\day_' num2str(day) '\trajectories\' dlc_files(file).name]);
            tmp_pellet_x = trial_kin(1:25,5);
            tmp_pellet_y = trial_kin(1:25,6);
            tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
            tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
            pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
            pellet_position = [pellet_position; pellet_location];
        end
        pellet_position = nanmean(pellet_position);
    
        mat_files = dir(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T107\reach_videos\day_' num2str(day) '\mat_videos\*' time_each_block{day}{blocks} '*']);
        [~, i] = sort_nat({mat_files(:).name});
        mat_files = mat_files(i);
        
        tmp_beh_markers = [];
        tmp_dlc_frame = [];
        tmp_wave_time = [];
        for file = 1:length(mat_files)
            trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T107\reach_videos\day_' num2str(day) '\trajectories\' dlc_files(file).name]);
            frame_times = load(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T107\reach_videos\day_' num2str(day) '\mat_videos\' mat_files(file).name]);
            frame_times = frame_times.t_lat;
            frame_times = frame_times*fs + start_times(file);
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
            if (i<31 || i>size(trial_kin,1)-31)
                tmp_beh_markers = [tmp_beh_markers nan];
            else
                tmp_beh_markers = [tmp_beh_markers frame_times(i)];
            end
            tmp_dlc_frame = [tmp_dlc_frame size(trial_kin,1)];
            tmp_wave_time = [tmp_wave_time pellet_drops(file)-start_times(file)];
        end
        beh_markers{blocks} = tmp_beh_markers/fs;
        
        figure 
            subplot(1,3,[1 2])
                hold on
                plot(wave);
                scatter(start_times,ones(1,length(start_times)))
                scatter(pellet_drops,ones(1,length(pellet_drops)))
                scatter(tmp_beh_markers,2*ones(1,length(tmp_beh_markers)))
                title(['Day ' num2str(day) ' Block ' num2str(blocks) ' - vids: ' num2str(length(mat_files)) ' - dlc files: ' num2str(length(dlc_files)) ' - wave starts: ' num2str(length(start_times))])
            subplot(1,3,3)
                hold on
                plot(zscore(tmp_dlc_frame));
                plot(zscore(tmp_wave_time));
                title('DLC frames vs. wave length');

    end

    [m1_reach, dls_reach] = reach_modulation(m1_spiking,dls_spiking,beh_markers);
    
    %%% SAVE
    if ~exist([save_path '\day_' num2str(day) '\'])
        mkdir([save_path '\day_' num2str(day) '\'])
    end
    cd([save_path '\day_' num2str(day) '\']);
    save('reach_mod.mat','m1_reach','dls_reach','-v7.3');

end

%% T200

clc; clear all; close all;

day_blocks = {[1], ...
              [2 7 8], ...
              [9 11 12 13 14], ...
              [15 16 17 18 19], ...
              [20 21 22 23 24], ...
              [25 26 27], ...
              [28 29 30], ...
              [32 33 34 35 36], ...
              [37 38 39 40 41], ...
              [42 43 44 45 46 47]};
reach_blocks = {[],[2],[3 4],[3 4],[2 3],[2],[2],[2 3],[2 3],[2]};

plx_files = {'T200_blocks_1-01', ...
             'T200_blocks_2_7_8-01', ...
             'T200_blocks_9_11_12_13_14-01', ...
             'T200_blocks_15_16_17_18_19-01', ...
             'T200_blocks_20_21_22_23_24-01', ...
             'T200_blocks_25_26_27-01', ...
             'T200_blocks_28_29_30-01', ...
             'T200_blocks_32_33_34_35_36-01', ...
             'T200_blocks_37_38_39_40_41-01', ...
             'T200_blocks_42_43_44_45_46_47-01'};

save_path = '\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T200\reach_mod';

for day = 2:10
    
    tmp_day_blocks = day_blocks{day};
    
    %%% LOAD SPIKING DATA
   
        plx_spiking = readtable(['\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T200\T200\' plx_files{day} '.txt']);
        plx_spiking = table2array(plx_spiking);
        
        tmp_ts = cell(32,length(tmp_day_blocks));
        block_count = 1;
        for block = tmp_day_blocks
            tmp_num = TDTbin2mat(['\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T200\T200\Block-' num2str(block)],'Type',3,'NODATA',1);
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
        
        m1_spiking = cell(length(reach_blocks{day}),16,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for m1_chan = 17:32
                for m1_unit = 1:max(plx_spiking(:,2))
                    m1_spiking{blocks,count,m1_unit} = tmp_ts{m1_chan,reach_blocks{day}(blocks)}(tmp_sort{m1_chan,reach_blocks{day}(blocks)}==m1_unit);
                end
                count = count + 1;
            end
        end
        
        dls_spiking = cell(length(reach_blocks{day}),16,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for dls_chan = 1:16
                for dls_unit = 1:max(plx_spiking(:,2))
                    dls_spiking{blocks,dls_chan,dls_unit} = tmp_ts{dls_chan,reach_blocks{day}(blocks)}(tmp_sort{dls_chan,reach_blocks{day}(blocks)}==dls_unit);
                end
                count = count + 1;
            end
        end

	%%% LOAD BEHAVIORAL MARKERS
    
        beh_markers = cell(length(reach_blocks{day}),1);
        dlc_files = dir(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T200\reach_videos\day_' num2str(day) '\trajectories\*corrected*']);
        [~, i] = sort_nat({dlc_files.name});
        dlc_files = dlc_files(i);     

        for blocks = 1:length(reach_blocks{day})

            tmp_block = tmp_day_blocks(reach_blocks{day}(blocks));

            %%% LOAD FRAME TIMES
            Vid = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T200\T200\Block-' num2str(tmp_block)],'STORE','Vid0');
            [~,i] = find(Vid.scalars.Vid0.data>1e8);
            for n_i = 1:length(i)
                Vid.scalars.Vid0.data(i(n_i)) = round(mean(Vid.scalars.Vid0.data(i(n_i)-1):Vid.scalars.Vid0.data(i(n_i)+1)));
            end  
            Vid.scalars.Vid0.data = Vid.scalars.Vid0.data+1;
            first_frames = find(diff(Vid.scalars.Vid0.ts)>8);
            trial_starts_idx = [Vid.scalars.Vid0.data(1) Vid.scalars.Vid0.data(first_frames+1)];
            trial_ends_idx = [Vid.scalars.Vid0.data(first_frames) Vid.scalars.Vid0.data(end)];      

            %%% GET KINEMATICS
            kin_all = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T200\reach_videos\day_' num2str(day) '\trajectories\' dlc_files(blocks).name]);

            pellet_position = [];
            for n_trial = 1:length(trial_starts_idx)
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

            tmp_beh_markers = [];
            for n_trial = 1:length(trial_starts_idx)
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
                if (i<31 || i>size(trial_kin,1)-31)
                    tmp_beh_markers = [tmp_beh_markers nan];
                else
                    tmp_beh_markers = [tmp_beh_markers interp_trial_ts(i)];
                end
           end

            beh_markers{blocks} = tmp_beh_markers;

        end

    [m1_reach, dls_reach] = reach_modulation(m1_spiking,dls_spiking,beh_markers);
    
    %%% SAVE
    if ~exist([save_path '\day_' num2str(day) '\'])
        mkdir([save_path '\day_' num2str(day) '\'])
    end
    cd([save_path '\day_' num2str(day) '\']);
    save('reach_mod.mat','m1_reach','dls_reach','-v7.3');
   
end

%% T201

clc; clear all; close all;

day_blocks = {[1], ...
              [2 3], ...
              [4 5 6 7], ...
              [8 9 10], ...
              [11 12 13], ...
              [14 15 16 17], ...
              [18 19 20], ...
              [21 22 23 24 25 26], ...
              [27 28 29], ...
              [30 31 32], ...
              [33 34 35 36 37], ...
              [38 39 40], ...
              [41 42 43 44], ...
              [45 46 47 48], ...
              [49 50 51 52], ...
              [53 54 55 56], ...
              [57 58 59 60 61]};

reach_blocks = {[], ...
                [], ...
                [3], ...
                [2], ...
                [2], ...
                [3], ...
                [2], ...
                [2 3 4], ...
                [2], ...
                [2], ...
                [3], ...
                [2], ...
                [2], ...
                [2], ...
                [2], ...
                [2], ...
                [3]};

plx_files = {'T201_blocks_1-01', ...
             'T201_blocks_2_3-01', ...
             'T201_blocks_4_5_6_7-01', ...
             'T201_blocks_8_9_10-01', ...
             'T201_blocks_11_12_13-01', ...
             'T201_blocks_14_15_16_17-01', ...
             'T201_blocks_18_19_20-01', ...
             'T201_blocks_21_22_23_24_25_26-01', ...
             'T201_blocks_27_28_29-01', ... 
             'T201_blocks_30_31_32-01', ...
             'T201_blocks_33_34_35_36_37-01', ...
             'T201_blocks_38_39_40-01', ...
             'T201_blocks_41_42_43_44-01', ...
             'T201_blocks_45_46_47_48-01', ...
             'T201_blocks_49_50_51_52-01', ...
             'T201_blocks_53_54_55_56-01', ...
             'T201_blocks_57_58_59_60_61-01'};

save_path = '\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T201\reach_mod';

for day = 4:17
    
    tmp_day_blocks = day_blocks{day};
    
    %%% LOAD SPIKING DATA
   
        plx_spiking = readtable(['\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T201\T201\' plx_files{day} '.txt']);
        plx_spiking = table2array(plx_spiking);
        
        tmp_ts = cell(64,length(tmp_day_blocks));
        block_count = 1;
        for block = tmp_day_blocks
            if block == 41
                tmp_num = TDTbin2mat(['\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T201\T201\Block-' num2str(block)],'Type',3,'T2',7200);
            else
                tmp_num = TDTbin2mat(['\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T201\T201\Block-' num2str(block)],'Type',3,'NODATA',1);
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
        
        m1_spiking = cell(length(reach_blocks{day}),32,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for m1_chan = 1:2:63
                for m1_unit = 1:max(plx_spiking(:,2))
                    m1_spiking{blocks,count,m1_unit} = tmp_ts{m1_chan,reach_blocks{day}(blocks)}(tmp_sort{m1_chan,reach_blocks{day}(blocks)}==m1_unit);
                end
                count = count + 1;
            end
        end
        
        dls_spiking = cell(length(reach_blocks{day}),16,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for dls_chan = 34:2:64
                for dls_unit = 1:max(plx_spiking(:,2))
                    dls_spiking{blocks,dls_chan,dls_unit} = tmp_ts{dls_chan,reach_blocks{day}(blocks)}(tmp_sort{dls_chan,reach_blocks{day}(blocks)}==dls_unit);
                end
                count = count + 1;
            end
        end

	%%% LOAD BEHAVIORAL MARKERS
    
        beh_markers = cell(length(reach_blocks{day}),1);
        dlc_files = dir(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T201\reach_videos\day_' num2str(day) '\trajectories\*corrected*']);
        [~, i] = sort_nat({dlc_files.name});
        dlc_files = dlc_files(i);     

        for blocks = 1:length(reach_blocks{day})

            tmp_block = tmp_day_blocks(reach_blocks{day}(blocks));

            %%% LOAD FRAME TIMES
            Vid = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T201\T201\Block-' num2str(tmp_block)],'STORE','Vid0');
            [~,i] = find(Vid.scalars.Vid0.data>1e8);
            for n_i = 1:length(i)
                Vid.scalars.Vid0.data(i(n_i)) = round(mean(Vid.scalars.Vid0.data(i(n_i)-1):Vid.scalars.Vid0.data(i(n_i)+1)));
            end  
            Vid.scalars.Vid0.data = Vid.scalars.Vid0.data+1;
            first_frames = find(diff(Vid.scalars.Vid0.ts)>8);
            trial_starts_idx = [Vid.scalars.Vid0.data(1) Vid.scalars.Vid0.data(first_frames+1)];
            trial_ends_idx = [Vid.scalars.Vid0.data(first_frames) Vid.scalars.Vid0.data(end)];      

            %%% GET KINEMATICS
            kin_all = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T201\reach_videos\day_' num2str(day) '\trajectories\' dlc_files(blocks).name]);

            if day==11 || day ==12 || day ==13 || day ==15 || day ==16 || day ==17
                trial_starts_idx = trial_starts_idx(1:end-1);
                trial_ends_idx = trial_ends_idx(1:end-1);
            end
                
            pellet_position = [];
            for n_trial = 1:length(trial_starts_idx)
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

            tmp_beh_markers = [];
            for n_trial = 1:length(trial_starts_idx)
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
                if (i<31 || i>size(trial_kin,1)-31)
                    tmp_beh_markers = [tmp_beh_markers nan];
                else
                    tmp_beh_markers = [tmp_beh_markers interp_trial_ts(i)];
                end
           end

            beh_markers{blocks} = tmp_beh_markers;

        end

    [m1_reach, dls_reach] = reach_modulation(m1_spiking,dls_spiking,beh_markers);
    
%     pause;
    
    %%% SAVE
    if ~exist([save_path '\day_' num2str(day) '\'])
        mkdir([save_path '\day_' num2str(day) '\'])
    end
    cd([save_path '\day_' num2str(day) '\']);
    save('reach_mod.mat','m1_reach','dls_reach','-v7.3');
   
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
            
reach_blocks = {[], ...
                [], ...
                [], ...
                [2], ...
                [2], ...
                [2], ...
                [2 3 4], ...
                [2], ...
                [3], ...
                [2], ...
                [2], ...
                [2], ...
                [2], ...
                [2]};

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

save_path = '\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\SLDD2\reach_mod';

for day = [4:14]
        
    %%% LOAD SPIKING DATA
   
        plx_spiking = readtable(['\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\SLDD2\sorted_spikes\day_' num2str(day) '\' plx_files{day} '.txt']);
        plx_spiking = table2array(plx_spiking);
        
        tmp_blocks = ls(['\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\SLDD2\sorted_spikes\day_' num2str(day) '\*.txt']);
        tmp_blocks_cell = cell(size(tmp_blocks,1),1);
        for block = 1:size(tmp_blocks,1)
            tmp_blocks_cell{block} = deblank(tmp_blocks(block,:));
        end
        tmp_ind = ~contains(tmp_blocks_cell,'mrg');
        tmp_blocks_cell = tmp_blocks_cell(tmp_ind);

        tmp_ts = cell(64,length(tmp_blocks_cell));
        for block = 1:length(tmp_blocks_cell)
            disp(['\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\SLDD2\sorted_spikes\day_' num2str(day) '\' tmp_blocks_cell{block}]);
            tmp_spiking = readtable(['\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\SLDD2\sorted_spikes\day_' num2str(day) '\' tmp_blocks_cell{block}]);
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
        
        m1_spiking = cell(length(reach_blocks{day}),32,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for m1_chan = [8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33]
                for m1_unit = 1:max(plx_spiking(:,2))
                    m1_spiking{blocks,count,m1_unit} = tmp_ts{m1_chan,reach_blocks{day}(blocks)}(tmp_sort{m1_chan,reach_blocks{day}(blocks)}==m1_unit);
                end
                count = count + 1;
            end
        end
        
        dls_spiking = cell(length(reach_blocks{day}),32,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for dls_chan = [2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63]
                for dls_unit = 1:max(plx_spiking(:,2))
                    dls_spiking{blocks,dls_chan,dls_unit} = tmp_ts{dls_chan,reach_blocks{day}(blocks)}(tmp_sort{dls_chan,reach_blocks{day}(blocks)}==dls_unit);
                end
                count = count + 1;
            end
        end

	%%% LOAD BEHAVIORAL MARKERS
    
        beh_markers = cell(length(reach_blocks{day}),1);

        for blocks = 1:length(reach_blocks{day})

            tmp_block = day_blocks{day}{reach_blocks{day}(blocks)};

            frames = TDTbin2mat(['\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\SLDD2\TDT\' tmp_block]);
            frames = frames.epocs.PtC1.onset;
            trial_starts = [1; find(diff(frames)>1)+1];
            trial_ends = [find(diff(frames)>1); length(frames)];

            if day == 7
                if blocks == 1
                    files = dir(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\SLDD2\reach_videos\day_' num2str(day-3) '\trajectories\*18h4m*corrected*']);
                    [~, i] = sort_nat({files(:).name});
                    files = files(i);
                    trial_starts = trial_starts([1:20 22:end]);
                    trial_ends = trial_ends([1:20 22:end]);
                elseif blocks == 2
                    files = dir(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\SLDD2\reach_videos\day_' num2str(day-3) '\trajectories\*18h37m*corrected*']);
                    [~, i] = sort_nat({files(:).name});
                    files = files(i);
                else
                    files = dir(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\SLDD2\reach_videos\day_' num2str(day-3) '\trajectories\*19h6m*corrected*']);
                    [~, i] = sort_nat({files(:).name});
                    files = files(i);
                end
            else
                files = dir(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\SLDD2\reach_videos\day_' num2str(day-3) '\trajectories\*corrected*']);
                [~, i] = sort_nat({files(:).name});
                files = files(i);
            end
            
            if day==7 && blocks ==2
                files = files(1:end-1);
            end
            
            %%% FIND PELLET POSITION
            pellet_position = [];
            for file = 1:length(files)
                trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\SLDD2\reach_videos\day_' num2str(day-3) '\trajectories\' files(file).name]);
                tmp_pellet_x = trial_kin(1:25,8);
                tmp_pellet_y = trial_kin(1:25,9);
                tmp_pellet_x(trial_kin(1:25,10)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,10)<.99) = NaN;
                pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                pellet_position = [pellet_position; pellet_location];
            end
            pellet_position = nanmean(pellet_position);
   
            tmp_tdt_length = [];
            tmp_dlc_length = [];
            tmp_beh_markers = [];
            for file = 1:length(files)
                frame_times = frames(trial_starts(file):trial_ends(file));
                trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\SLDD2\reach_videos\day_' num2str(day-3) '\trajectories\' files(file).name]);
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
                if (i<31 || i>size(trial_kin,1)-31)
                    tmp_beh_markers = [tmp_beh_markers nan];
                else
                    tmp_beh_markers = [tmp_beh_markers frame_times(i)];
                end
                tmp_dlc_length = [tmp_dlc_length size(trial_kin,1)];
                tmp_tdt_length = [tmp_tdt_length frames(trial_ends(file))-frames(trial_starts(file))];
            end
            beh_markers{blocks} = tmp_beh_markers;
            
            figure
                hold on
                plot(zscore(tmp_dlc_length));
                plot(zscore(tmp_tdt_length));
                title(['Day ' num2str(day) ' Block ' num2str(blocks) ' - dlc files: ' num2str(length(files)) ' - wave starts: ' num2str(length(trial_starts))])
        end

    [m1_reach, dls_reach] = reach_modulation(m1_spiking,dls_spiking,beh_markers);
    
    %%% SAVE
    if ~exist([save_path '\day_' num2str(day) '\'])
        mkdir([save_path '\day_' num2str(day) '\'])
    end
    cd([save_path '\day_' num2str(day) '\']);
    save('reach_mod.mat','m1_reach','dls_reach','-v7.3');
   
    pause(0.1);
    
end

%% T398
clc; clear all; close all;

spike_path = {['E:\Kilosort\T398_kilosort\2.26.19'], ...
    ['E:\Kilosort\T398_kilosort\2.27.19'], ...
    ['E:\Kilosort\T398_kilosort\2.28.19'], ...
    ['E:\Kilosort\T398_kilosort\3.1.19'], ...
    ['E:\Kilosort\T398_kilosort\3.2.19'], ...
    ['E:\Kilosort\T398_kilosort\3.3.19'], ...
    ['E:\Kilosort\T398_kilosort\3.4.19'], ...
    ['E:\Kilosort\T398_kilosort\3.5.19'], ...
    ['E:\Kilosort\T398_kilosort\3.6.19']};

pre_sleep = {'T398-190226-112203','T398-190227-113133','T398-190228-105336','T398-190301-112846','T398-190302-131102','T398-190303-140804','T398-190304-114807','T398-190305-130625','T398-190306-091252'};
reach_blocks = {'T398-190226-151524','T398-190227-144600','T398-190228-125951','T398-190301-135817','T398-190302-151642','T398-190303-161420','T398-190304-135101','T398-190305-151059','T398-190306-111957'};
post_sleep = {'T398-190226-160611','T398-190227-154618','T398-190228-141352','T398-190301-151942','T398-190302-162611','T398-190303-171207','T398-190304-145843','T398-190305-161031','T398-190306-124648'};

save_path = '\\mycloudpr4100\public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T398\reach_mod';

for day = [3:5 7:9]
      
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
        pre_sleep_length = TDTbin2mat(['F:\T398\RS4\M1-DLS_256Probe-190222-144047\' pre_sleep{day}],'CHANNEL',1);
        samp_freq = pre_sleep_length.RSn1.fs;
        pre_sleep_length = length(pre_sleep_length.RSn1.data);

        reach_length = TDTbin2mat(['F:\T398\RS4\M1-DLS_256Probe-190222-144047\' reach_blocks{day}],'CHANNEL',1);
        reach_length = length(reach_length.RSn1.data);

        post_sleep_length = TDTbin2mat(['F:\T398\RS4\M1-DLS_256Probe-190222-144047\' post_sleep{day}],'CHANNEL',1);
        post_sleep_length = length(post_sleep_length.RSn1.data);

        max(m1_ts)
        max(dls_ts)
        pre_sleep_length + reach_length + post_sleep_length
        %%% CHECK UNIT %%%
        
        m1_spiking = cell(length(reach_blocks{day}),length(m1_clust_id),1);
        for blocks = 1:length(reach_blocks{day})
            m1_count = 1;
            for m1_n = m1_clust_id
                tmp_spikes = m1_ts(m1_id==m1_n);
                m1_spiking{blocks,m1_count,1} = double(tmp_spikes(tmp_spikes>pre_sleep_length & tmp_spikes<pre_sleep_length+reach_length)-pre_sleep_length)/samp_freq;
                m1_count = m1_count+1;
            end
        end
        
        dls_spiking = cell(length(reach_blocks{day}),length(dls_clust_id),1);
        for blocks = 1:length(reach_blocks{day})
            dls_count = 1;
            for dls_n = dls_clust_id
                tmp_spikes = dls_ts(dls_id==dls_n);
                dls_spiking{blocks,dls_count,1} = double(tmp_spikes(tmp_spikes>pre_sleep_length & tmp_spikes<pre_sleep_length+reach_length)-pre_sleep_length)/samp_freq;
                dls_count = dls_count+1;
            end
        end
        
	%%% LOAD BEHAVIORAL MARKERS
    
        beh_markers = cell(1,1);
        
        files = dir(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T398\reach_videos\day_' num2str(day) '\trajectories\*T398*']);
        [~, i] = sort_nat({files(:).name});
        files = files(i);
                
        %%% FIND PELLET POSITION
        pellet_position = [];
        for file = 1:length(files)
            trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T398\reach_videos\day_' num2str(day) '\trajectories\' files(file).name]);
            tmp_pellet_x = trial_kin(1:25,8);
            tmp_pellet_y = trial_kin(1:25,9);
            tmp_pellet_x(trial_kin(1:25,10)<.99) = NaN;
            tmp_pellet_y(trial_kin(1:25,10)<.99) = NaN;
            pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
            pellet_position = [pellet_position; pellet_location];
        end
        pellet_position = nanmean(pellet_position);

        frames = TDTbin2mat(['\\mycloudpr4100\public\SL_Backups\T398\TDT\M1-DLS_256Probe-190222-144047\' reach_blocks{day}]);
        frames = frames.epocs.PC2_.onset;
        trial_starts = [1; find(diff(frames)>1)+1];
        trial_ends = [find(diff(frames)>1); length(frames)];
        
        if day == 3
            trial_starts = trial_starts(2:end-1);
            trial_ends = trial_ends(2:end-1);
        end
        
        tmp_tdt_length = [];
        tmp_dlc_length = [];
        tmp_beh_markers = [];
        for file = 1:length(files)
            frame_times = frames(trial_starts(file):trial_ends(file));
            trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\M1_DLS_SLEEP\data\T398\reach_videos\day_' num2str(day) '\trajectories\' files(file).name]);
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
            if (i<31 || i>size(trial_kin,1)-31)
                tmp_beh_markers = [tmp_beh_markers nan];
            else
                tmp_beh_markers = [tmp_beh_markers frame_times(i)];
            end
            tmp_dlc_length = [tmp_dlc_length size(trial_kin,1)];
            tmp_tdt_length = [tmp_tdt_length frames(trial_ends(file))-frames(trial_starts(file))];
        end
        beh_markers{1} = tmp_beh_markers;
            
        %%% CHECK %%%
        figure
        hold on
        plot(zscore(tmp_dlc_length));
        plot(zscore(tmp_tdt_length));
        title(['Day ' num2str(day) ' Block ' num2str(blocks) ' - dlc files: ' num2str(length(files)) ' - wave starts: ' num2str(length(trial_starts))])
        %%% CHECK %%%

    [m1_reach, dls_reach] = reach_modulation(m1_spiking,dls_spiking,beh_markers);
    
    %%% SAVE
    if ~exist([save_path '\day_' num2str(day) '\'])
        mkdir([save_path '\day_' num2str(day) '\'])
    end
    cd([save_path '\day_' num2str(day) '\']);
    save('reach_mod.mat','m1_reach','dls_reach','-v7.3');
   
    pause(0.1);
end
