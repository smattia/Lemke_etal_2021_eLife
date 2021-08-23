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

pre_LFP_coh = cell(1,256);
post_LFP_coh = cell(1,256);

figure
hold on

for day = 1:8
    
    tmp_day_blocks = day_blocks{day};
    
    %%% LOAD BEH STATE, LFP, AND START/STOP

        load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T102\sleep_activity\day_' num2str(day) '_block_1_sleep_activity.mat'])
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

        if day == 6
            load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T102\sleep_activity\day_' num2str(day) '_block_3_sleep_activity.mat'])
        else
            load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T102\sleep_activity\day_' num2str(day) '_block_2_sleep_activity.mat'])
        end
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

    %%% FUNCTION START
    
    params.Fs = data.streams.Wave.fs;
    params.tapers =[3 5];
    params.fpass = [1 15];
    params.err = [2 0.05];
    params.trialave = 1;
    params.pad = 0;
    movingwin = [10 10];
    
    pre_coh = [];
    post_coh = [];
    
    pair_count = 1;
    
    for m1_chan = 1:16
        
        disp(['M1 chan ' num2str(m1_chan)])
        
        tmp_pre_m1_lfp = pre_lfp_norm(m1_chan,:);
        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
        
        tmp_post_m1_lfp = post_lfp_norm(m1_chan,:);
        tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop(1):post_video_start_stop(2));
        tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
        
        for dls_chan = 17:32
            
            tmp_pre_dls_lfp = pre_lfp_norm(dls_chan,:);
            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_beh_state.nrem_interp==1);
            
            tmp_post_dls_lfp = post_lfp_norm(dls_chan,:);
            tmp_post_dls_lfp = tmp_post_dls_lfp(post_video_start_stop(1):post_video_start_stop(2));
            tmp_post_dls_lfp = tmp_post_dls_lfp(post_beh_state.nrem_interp==1);
            
            [tmp_coh,~,~,~,~,times,freqs,~,~,~]=cohgramc(tmp_pre_m1_lfp',tmp_pre_dls_lfp',movingwin,params);
            pre_LFP_coh{pair_count} = [pre_LFP_coh{pair_count}; mean(tmp_coh)];
            pre_coh = [pre_coh; mean(tmp_coh)];
            
            [tmp_coh,~,~,~,~,times,freqs,~,~,~]=cohgramc(tmp_post_m1_lfp',tmp_post_dls_lfp',movingwin,params);
            post_LFP_coh{pair_count} = [post_LFP_coh{pair_count}; mean(tmp_coh)];
            post_coh = [post_coh; mean(tmp_coh)];
            
            pair_count = pair_count + 1;
            
        end
    end
    
    plot(mean(pre_coh),'color',[day/8 0 1-(day/8)])
    plot(mean(post_coh),'color',[day/8 0 1-(day/8)])
    
end

%%% SAVE

    save_path = 'C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T102\sleep_activity';
    if ~exist([save_path '\'])
        mkdir([save_path '\'])
    end
    cd([save_path '\']);
    save('LFP_coherence.mat','pre_LFP_coh','post_LFP_coh','times','freqs','-v7.3');

%% T107

clc; clear all; close all;

day_blocks = {[1 2],[3 4 5 6 8],[9 10 11 12 13 14 15],[16 17 18 19 20 21 22 24 25],[26 27 30 31 36 39],[41 42 43 47 49 51 52]};
pre_sleep_blocks = {[1],[1],[1],[3],[1],[1]};
post_sleep_blocks = {[2],[4],[6],[8],[5],[4]};

pre_LFP_coh = cell(1,256);
post_LFP_coh = cell(1,256);

figure
hold on

for day = 1:6
    
    tmp_day_blocks = day_blocks{day};
    
    %%% LOAD BEH STATE, LFP, AND START/STOP
    
    load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T107\sleep_activity\day_' num2str(day) '_block_1_sleep_activity.mat'])
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
    
    load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T107\sleep_activity\day_' num2str(day) '_block_2_sleep_activity.mat'])
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
    
    %%% FUNCTION START
    
    params.Fs= data.streams.Wave.fs;
    params.tapers=[3 5];
    params.fpass = [1 15];
    params.err = [2 0.05];
    params.trialave = 1;
    params.pad = 0;
    movingwin = [10 10];
    
    pre_coh = [];
    post_coh = [];
    
    pair_count = 1;
    
    for m1_chan = 1:16
        
        disp(['M1 chan ' num2str(m1_chan)])
        
        tmp_pre_m1_lfp = pre_lfp_norm(m1_chan,:);
        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
        
        tmp_post_m1_lfp = post_lfp_norm(m1_chan,:);
        tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop(1):post_video_start_stop(2));
        tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
        
        for dls_chan = 17:32
            
            tmp_pre_dls_lfp = pre_lfp_norm(dls_chan,:);
            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_beh_state.nrem_interp==1);
            
            tmp_post_dls_lfp = post_lfp_norm(dls_chan,:);
            tmp_post_dls_lfp = tmp_post_dls_lfp(post_video_start_stop(1):post_video_start_stop(2));
            tmp_post_dls_lfp = tmp_post_dls_lfp(post_beh_state.nrem_interp==1);
            
            [tmp_coh,~,~,~,~,times,freqs,~,~,~]=cohgramc(tmp_pre_m1_lfp',tmp_pre_dls_lfp',movingwin,params);
            pre_LFP_coh{pair_count} = [pre_LFP_coh{pair_count}; mean(tmp_coh)];
            pre_coh = [pre_coh; mean(tmp_coh)];
            
            [tmp_coh,~,~,~,~,times,freqs,~,~,~]=cohgramc(tmp_post_m1_lfp',tmp_post_dls_lfp',movingwin,params);
            post_LFP_coh{pair_count} = [post_LFP_coh{pair_count}; mean(tmp_coh)];
            post_coh = [post_coh; mean(tmp_coh)];
            
            pair_count = pair_count + 1;
            
        end
    end
    
    plot(mean(pre_coh),'color',[day/6 0 1-(day/6)])
    plot(mean(post_coh),'color',[day/6 0 1-(day/6)])
    
end

%%% SAVE

    save_path = 'C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T107\sleep_activity';
    if ~exist([save_path '\'])
        mkdir([save_path '\'])
    end
    cd([save_path '\']);
    save('LFP_coherence.mat','pre_LFP_coh','post_LFP_coh','times','freqs','-v7.3');

%% T200

clc; clear all; close all;

day_blocks = {[1],[2 7 8],[9 11 12 13 14],[15 16 17 18 19],[20 21 22 23 24],[25 26 27],[28 29 30],[32 33 34 35 36],[37 38 39 40 41],[42 43 44 45 46 47]};
pre_sleep_blocks = {[1],[1],[2],[1],[1],[1],[1],[1],[1],[1]};
post_sleep_blocks = {[],[3],[5],[5],[4],[3],[3],[5],[5],[6]};

pre_LFP_coh = cell(1,256);
post_LFP_coh = cell(1,256);

figure
hold on

for day = 2:10
    
    tmp_day_blocks = day_blocks{day};
    
    %%% LOAD BEH STATE AND START/STOP
    
    load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T200\sleep_activity\day_' num2str(day) '_block_1_sleep_activity.mat'])
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
    
    load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T200\sleep_activity\day_' num2str(day) '_block_2_sleep_activity.mat'])
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
    
    %%% FUNCTION START
    
    params.Fs= data.streams.Wave.fs;
    params.tapers=[3 5];
    params.fpass = [1 15];
    params.err = [2 0.05];
    params.trialave = 1;
    params.pad = 0;
    movingwin = [10 10];
    
    pre_coh = [];
    post_coh = [];
    
    pair_count = 1;
    
    for m1_chan = 17:32
        
        disp(['M1 chan ' num2str(m1_chan)])
        
        tmp_pre_m1_lfp = pre_lfp_norm(m1_chan,:);
        if day~=6
            tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
        end
        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
        
        tmp_post_m1_lfp = post_lfp_norm(m1_chan,:);
        tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop(1):post_video_start_stop(2));
        tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
        
        for dls_chan = 1:16
            
            tmp_pre_dls_lfp = pre_lfp_norm(dls_chan,:);
            if day~=6
                tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
            end            
            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_beh_state.nrem_interp==1);
            
            tmp_post_dls_lfp = post_lfp_norm(dls_chan,:);
            tmp_post_dls_lfp = tmp_post_dls_lfp(post_video_start_stop(1):post_video_start_stop(2));
            tmp_post_dls_lfp = tmp_post_dls_lfp(post_beh_state.nrem_interp==1);
            
            [tmp_coh,~,~,~,~,times,freqs,~,~,~]=cohgramc(tmp_pre_m1_lfp',tmp_pre_dls_lfp',movingwin,params);
            pre_LFP_coh{pair_count} = [pre_LFP_coh{pair_count}; mean(tmp_coh)];
            pre_coh = [pre_coh; mean(tmp_coh)];
            
            [tmp_coh,~,~,~,~,times,freqs,~,~,~]=cohgramc(tmp_post_m1_lfp',tmp_post_dls_lfp',movingwin,params);
            post_LFP_coh{pair_count} = [post_LFP_coh{pair_count}; mean(tmp_coh)];
            post_coh = [post_coh; mean(tmp_coh)];
            
            pair_count = pair_count + 1;
            
        end
    end
    
    plot(mean(pre_coh),'color',[day/10 0 1-(day/10)])
    plot(mean(post_coh),'color',[day/10 0 1-(day/10)])
    
end

%%% SAVE

    save_path = 'C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T200\sleep_activity';
    if ~exist([save_path '\'])
        mkdir([save_path '\'])
    end
    cd([save_path '\']);
    save('LFP_coherence.mat','pre_LFP_coh','post_LFP_coh','times','freqs','-v7.3');

%% T201

clc; clear all; close all;

day_blocks = {[1],[2 3],[4 5 6 7],[8 9 10],[11 12 13],[14 15 16 17],[18 19 20],[21 22 23 24 25 26],[27 28 29], ...
    [30 31 32],[33 34 35 36 37],[38 39 40],[41 42 43 44],[45 46 47 48],[49 50 51 52],[53 54 55 56],[57 58 59 60 61]};

pre_sleep_blocks = {[1],[1],[2],[1],[1],[2],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1]};
post_sleep_blocks = {[],[2],[4],[3],[3],[4],[3],[6],[3],[3],[5],[3],[4],[4],[4],[4],[5]};

pre_LFP_coh = cell(1,512);
post_LFP_coh = cell(1,512);

for day = 3:17
    
    tmp_day_blocks = day_blocks{day};
    
    %%% LOAD BEH STATE AND START/STOP
    
    load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T201\sleep_activity\day_' num2str(day) '_block_1_sleep_activity.mat'])
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
    
    load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T201\sleep_activity\day_' num2str(day) '_block_2_sleep_activity.mat'])
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
    
    %%% FUNCTION START
    
    params.Fs= data.streams.Wave.fs;
    params.tapers=[3 5];
    params.fpass = [1 15];
    params.err = [2 0.05];
    params.trialave = 1;
    params.pad = 0;
    movingwin = [10 10];
    
    pre_coh = [];
    post_coh = [];
    
    pair_count = 1;
    
    for m1_chan = 1:2:63
        
        disp(['M1 chan ' num2str(m1_chan)])
        
        tmp_pre_m1_lfp = pre_lfp_norm(m1_chan,:);
        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
        
        tmp_post_m1_lfp = post_lfp_norm(m1_chan,:);
        tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop(1):post_video_start_stop(2));
        tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
        
        for dls_chan = 34:2:64
            
            tmp_pre_dls_lfp = pre_lfp_norm(dls_chan,:);
            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_beh_state.nrem_interp==1);
            
            tmp_post_dls_lfp = post_lfp_norm(dls_chan,:);
            tmp_post_dls_lfp = tmp_post_dls_lfp(post_video_start_stop(1):post_video_start_stop(2));
            tmp_post_dls_lfp = tmp_post_dls_lfp(post_beh_state.nrem_interp==1);
            
            [tmp_coh,~,~,~,~,times,freqs,~,~,~]=cohgramc(tmp_pre_m1_lfp',tmp_pre_dls_lfp',movingwin,params);
            pre_LFP_coh{pair_count} = [pre_LFP_coh{pair_count}; mean(tmp_coh)];
            pre_coh = [pre_coh; mean(tmp_coh)];
            
            [tmp_coh,~,~,~,~,times,freqs,~,~,~]=cohgramc(tmp_post_m1_lfp',tmp_post_dls_lfp',movingwin,params);
            post_LFP_coh{pair_count} = [post_LFP_coh{pair_count}; mean(tmp_coh)];
            post_coh = [post_coh; mean(tmp_coh)];
            
            pair_count = pair_count + 1;
            
        end
    end
    
    plot(mean(pre_coh),'color',[day/17 0 1-(day/17)])
    plot(mean(post_coh),'color',[day/17 0 1-(day/17)])
    
end

%%% SAVE

    save_path = 'C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T201\sleep_activity';
    if ~exist([save_path '\'])
        mkdir([save_path '\'])
    end
    cd([save_path '\']);
    save('LFP_coherence.mat','pre_LFP_coh','post_LFP_coh','times','freqs','-v7.3');

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

pre_LFP_coh = cell(1,1024);
post_LFP_coh = cell(1,1024);

for day = 4:13
    
    tmp_sleep_blocks = day_blocks{day};
    
    %%% LOAD BEH STATE AND START/STOP
    
    load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\SLDD2\sleep_activity\day_' num2str(day) '_block_1_sleep_activity.mat'])
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
    
	load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\SLDD2\sleep_activity\day_' num2str(day) '_block_2_sleep_activity.mat'])
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
    
    %%% FUNCTION START
    
    pre_beh_state.nrem_interp = pre_beh_state.nrem_interp(1:length(pre_video_start_stop(1):pre_video_start_stop(2)));
    post_beh_state.nrem_interp = post_beh_state.nrem_interp(1:length(post_video_start_stop(1):post_video_start_stop(2)));
    
    params.Fs= data.RSn1.fs/20;
    params.tapers=[3 5];
    params.fpass = [1 15];
    params.err = [2 0.05];
    params.trialave = 1;
    params.pad = 0;
    movingwin = [10 10];
    
    pre_coh = [];
    post_coh = [];
    
    pair_count = 1;
    
    for m1_chan = [8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33]
        
        disp(['M1 chan ' num2str(m1_chan)])
        
        tmp_pre_m1_lfp = pre_lfp_norm(m1_chan,:);
        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
        tmp_pre_m1_lfp = tmp_pre_m1_lfp(pre_beh_state.nrem_interp==1);
        
        tmp_post_m1_lfp = post_lfp_norm(m1_chan,:);
        tmp_post_m1_lfp = tmp_post_m1_lfp(post_video_start_stop(1):post_video_start_stop(2));
        tmp_post_m1_lfp = tmp_post_m1_lfp(post_beh_state.nrem_interp==1);
        
        for dls_chan = [2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63]
            
            tmp_pre_dls_lfp = pre_lfp_norm(dls_chan,:);
            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_video_start_stop(1):pre_video_start_stop(2));
            tmp_pre_dls_lfp = tmp_pre_dls_lfp(pre_beh_state.nrem_interp==1);
            
            tmp_post_dls_lfp = post_lfp_norm(dls_chan,:);
            tmp_post_dls_lfp = tmp_post_dls_lfp(post_video_start_stop(1):post_video_start_stop(2));
            tmp_post_dls_lfp = tmp_post_dls_lfp(post_beh_state.nrem_interp==1);
            
            [tmp_coh,~,~,~,~,times,freqs,~,~,~]=cohgramc(tmp_pre_m1_lfp',tmp_pre_dls_lfp',movingwin,params);
            pre_LFP_coh{pair_count} = [pre_LFP_coh{pair_count}; mean(tmp_coh)];
            pre_coh = [pre_coh; mean(tmp_coh)];
            
            [tmp_coh,~,~,~,~,times,freqs,~,~,~]=cohgramc(tmp_post_m1_lfp',tmp_post_dls_lfp',movingwin,params);
            post_LFP_coh{pair_count} = [post_LFP_coh{pair_count}; mean(tmp_coh)];
            post_coh = [post_coh; mean(tmp_coh)];
            
            pair_count = pair_count + 1;
            
        end
    end
    
    figure;
    hold on;
    plot(mean(pre_coh),'color',[0 0 0])
    plot(mean(post_coh),'color',[1 0 0])
    pause(0.1)
    
end

%%% SAVE

    save_path = 'C:\Users\Stefan\Desktop\Lemke_etal_2021\data\SLDD2\sleep_activity';
    if ~exist([save_path '\'])
        mkdir([save_path '\'])
    end
    cd([save_path '\']);
    save('LFP_coherence.mat','pre_LFP_coh','post_LFP_coh','times','freqs','-v7.3');

%% T398

clc; clear all; close all;

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

pre_LFP_coh = cell(1,1406);
post_LFP_coh = cell(1,1406);

for day = [1:5 7:9]
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
        
        load(['C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T398\sleep_activity\day_' num2str(day) '_block_' num2str(block) '_sleep_activity.mat'])
        
        if day==7 && block==2
            video_start_stop = [10*(data.RSn1.fs/20) 7210*(data.RSn1.fs/20)];
        else
            wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\TDT_tanks\M1-DLS_256Probe-190222-144047\' sleep_blocks{day}{block}],'TYPE',[2]);
            video_start_stop = [min(wave.epocs.PC1_.onset)*(data.RSn1.fs/20) max(wave.epocs.PC1_.onset)*(data.RSn1.fs/20)];
            if day == 2 && block == 1
                video_start_stop(2) = video_start_stop(2)+7200*(data.RSn1.fs/20);
            end
        end
        
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
        
        beh_state.nrem_interp = beh_state.nrem_interp(1:length(video_start_stop(1):video_start_stop(2)));
        
        params.Fs= data.RSn1.fs/20;
        params.tapers=[3 5];
        params.fpass = [1 15];
        params.err = [2 0.05];
        params.trialave = 1;
        params.pad = 0;
        movingwin = [10 10];
        
        pre_coh = [];
        post_coh = [];
        
        pair_count = 1;
        
        for m1_chan = 1:size(m1_lfp_norm,1)
            
            disp(['M1 chan ' num2str(m1_chan)])
            
            tmp_m1_lfp = m1_lfp_norm(m1_chan,:);
            if day ~= 6
                tmp_m1_lfp = tmp_m1_lfp(video_start_stop(1):video_start_stop(2));
            end
            tmp_m1_lfp = tmp_m1_lfp(beh_state.nrem_interp==1);
            
            for dls_chan = 1:size(dls_lfp_norm,1)
                
                tmp_dls_lfp = dls_lfp_norm(dls_chan,:);
                if day ~= 6
                    tmp_dls_lfp = tmp_dls_lfp(video_start_stop(1):video_start_stop(2));
                end
                tmp_dls_lfp = tmp_dls_lfp(beh_state.nrem_interp==1);
                
                [tmp_coh,~,~,~,~,times,freqs,~,~,~]=cohgramc(tmp_m1_lfp',tmp_dls_lfp',movingwin,params);
                
                if block == 1
                    pre_LFP_coh{pair_count} = [pre_LFP_coh{pair_count}; mean(tmp_coh)];
                    pre_coh = [pre_coh; mean(tmp_coh)];
                else
                    post_LFP_coh{pair_count} = [post_LFP_coh{pair_count}; mean(tmp_coh)];
                    post_coh = [post_coh; mean(tmp_coh)];
                end
                
                pair_count = pair_count + 1;
                
            end
        end
        
        figure;
        hold on;
        plot(mean(pre_coh),'color',[0 0 0])
        plot(mean(post_coh),'color',[1 0 0])
        pause(0.1)
    end
end

%%% SAVE

    save_path = 'C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T398\sleep_activity';
    if ~exist([save_path '\'])
        mkdir([save_path '\'])
    end
    cd([save_path '\']);
    save('LFP_coherence.mat','pre_LFP_coh','post_LFP_coh','times','freqs','-v7.3');
