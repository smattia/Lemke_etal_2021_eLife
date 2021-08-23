function [beh_state, sleep_rhythms, mean_lfp] = detect_sleep(m1_lfp,dls_lfp,lfp_samp_rate, video, video_start_stop, sensitivity, day, block, varargin)
%   detects behavioral state and finds sleep rhythms - sleep spindles, slow
%   oscillations, and delta waves - in NREM sleep (SL 11/2019)
%
%   INPUTS:
%       m1_lfp = M1 LFP data (channels x timebins)
%       dls_lfp = DLS LFP data (channels x timebins)
%       lfp_samp_rate = LFP sampling rate
%       video = path to video
%       video_start_stop = start and stop time for video (seconds)
%
%   OPTIONAL INPUTS:
%       emg = EMG data (channels x timebins)
%       emg_samp_rate = EMG samping rate
%       save_path = path to save data & plot (will save if given)
%       remove artifcat LFP = 1 or [] (default: no removal)
%
%   OUTPUTS:
%       beh_state - struct of behavioral state classification
%       sleep_rhythms - struct of spindles, so, and delta
%       mean_lfp - mean lfp from M1 and DLS
%

%% ARTIFACT REJECT AND Z-SCORE LFP

    for chan = 1:size(m1_lfp,1) 
        artf = abs(zscore(m1_lfp(chan,:)))>6; 
        no_artf_mean = mean(m1_lfp(chan,~artf)); 
        no_artf_std = std(m1_lfp(chan,~artf));
        m1_lfp(chan,artf) = no_artf_mean;
        m1_lfp(chan,:) = (m1_lfp(chan,:)-no_artf_mean)/no_artf_std;
    end
    
    for chan = 1:size(dls_lfp,1) 
        artf = abs(zscore(dls_lfp(chan,:)))>6; 
        no_artf_mean = mean(dls_lfp(chan,~artf)); 
        no_artf_std = std(dls_lfp(chan,~artf));
        dls_lfp(chan,artf) = no_artf_mean;
        dls_lfp(chan,:) = (dls_lfp(chan,:)-no_artf_mean)/no_artf_std;
    end
    
    if ~isempty(video_start_stop) && exist(video)==2
        m1_lfp = double(mean(m1_lfp(:,round(video_start_stop(1)*lfp_samp_rate):round(video_start_stop(2)*lfp_samp_rate))));
        dls_lfp = double(mean(dls_lfp(:,round(video_start_stop(1)*lfp_samp_rate):round(video_start_stop(2)*lfp_samp_rate))));
    else
        m1_lfp = double(mean(m1_lfp));
        dls_lfp = double(mean(dls_lfp));
    end
    
%% CALCULATE LFP POWER

    params.Fs = lfp_samp_rate;
    params.tapers = [5 9];
    params.fpass = [.1 30];
    params.err = [2 0.05];
    params.trialave = 0;
    params.pad = 0;
    movingwin = [10 10];
    [lfp_power,t,f,~]=mtspecgramc(m1_lfp,movingwin,params);
    % 1-4 Hz
    [~,low] = min(abs(f-1));
    [~,high] = min(abs(f-4));
    delta_power = zscore(mean(lfp_power(:,low:high),2)); 
    
    [~,low1] = min(abs(f-5));
    [~,high1] = min(abs(f-10));
    [~,low2] = min(abs(f-2));
    [~,high2] = min(abs(f-15));
    theta_power = zscore(mean(lfp_power(:,low1:high1),2)./mean(lfp_power(:,low2:high2),2)); % 5-10 Hz / 2-15 Hz

%% CALCULATE MOVEMENT/EMG

    movement = zeros(1,length(delta_power));
    
    if exist(video)==2 && ~isempty(video_start_stop)
        vid = VideoReader(video);
        movement_tmp = zeros(1,vid.duration*vid.FrameRate);
        Frame1 = rgb2gray(readFrame(vid));
        count = 1;
        while hasFrame(vid)
            Frame2 = rgb2gray(readFrame(vid));
            movement_tmp(count) = nnz(imbinarize(Frame2-Frame1,sensitivity));
            Frame1 = Frame2;
            count = count + 1;
        end
        movement_tmp(zscore(movement_tmp)>6) = mean(movement_tmp);
        movement_tmp = zscore(movement_tmp);

        step_size = length(movement_tmp)/length(delta_power);
        for frame=0:length(delta_power)-1
            movement(frame+1) = mean(movement_tmp(floor(1+(step_size*frame)):floor(step_size+(step_size*frame))));
        end
    end
    
    if varargin{4} == 1 
        delta_power(movement>0)=-1;
        delta_power = zscore(delta_power);
    end

%% CLASSIFY SLEEP

    tmp_index_nrem = zeros(1,length(delta_power));
    tmp_index_nrem(delta_power'>0 & movement<=0) = 1;
    tmp_index_rem = zeros(1,length(delta_power));
    tmp_index_rem(theta_power'>0 & delta_power'<0 & movement<0) = 1;

    beh_state.nrem = zeros(1,length(delta_power));
    beh_state.rem = zeros(1,length(delta_power));
    beh_state.wake = zeros(1,length(delta_power));
    for n=1:length(delta_power)-5
        if sum(tmp_index_nrem(n:n+5))==6
            beh_state.nrem(n:n+5)=1;
        elseif sum(tmp_index_rem(n:n+5))==6
            beh_state.rem(n:n+5)=1;
        else
            beh_state.wake(n:n+5)=1;
        end
    end
    
    beh_state.delta_power = delta_power;
    
    beh_state.nrem_interp = interp1(1:length(beh_state.nrem),beh_state.nrem,linspace(1,length(beh_state.nrem),length(m1_lfp)));
    beh_state.rem_interp = interp1(1:length(beh_state.rem),beh_state.rem,linspace(1,length(beh_state.rem),length(m1_lfp)));
    beh_state.wake_interp = interp1(1:length(beh_state.wake),beh_state.wake,linspace(1,length(beh_state.wake),length(m1_lfp)));

%% DETECT NREM SLEEP RHYTHMS

    sleep_rhythms.spindles = detect_spindles({m1_lfp},lfp_samp_rate,{beh_state.nrem_interp});
    sleep_rhythms.so_delta = detect_slow_oscillations(m1_lfp',lfp_samp_rate,beh_state.nrem_interp');
    
%% PLOT

fig = figure('units','normalized','outerposition',[0 0 1 1]);

    subplot(4,2,1)
        hold on;
        imagesc(log(lfp_power(:,1:end)'));
        colorbar;
        for n=1:length(delta_power)
            if beh_state.nrem(n)==1
                fill([n n n+1 n+1], [1 size(lfp_power,2)-15 size(lfp_power,2)-15 1], [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2)
            end
        end
        for n=1:length(delta_power)
            if beh_state.rem(n)==1
                fill([n n n+1 n+1], [1 size(lfp_power,2)-15 size(lfp_power,2)-15 1], [0 0 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2)
            end
        end
        xlim([1 size(lfp_power,1)])
        ylim([1 size(lfp_power,2)-15])       

    subplot(4,2,3)
        hold on;
        plot(delta_power,'color','r','LineWidth',1);
        plot(theta_power,'color','b','LineWidth',1);
        plot(movement,'color','k','LineWidth',1);
        y_lim  = ylim;
        for n=1:length(delta_power)
            if beh_state.nrem(n)==1
                fill([n n n+1 n+1], [y_lim(1) y_lim(2) y_lim(2) y_lim(1)], [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3)
            end
        end
        for n=1:length(delta_power)
            if beh_state.rem(n)==1
                fill([n n n+1 n+1], [y_lim(1) y_lim(2) y_lim(2) y_lim(1)], [0 0 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3)
            end
        end
        xlim([1 size(lfp_power,1)])
        ylim([y_lim(1) y_lim(2)])
        if max([y_lim(1) y_lim(2)])>5
            ylim([-5 5])
        end
        ylabel('delta / theta / movement')

    subplot(4,2,[5 7])
        hold on;
        for n=1:length(beh_state.nrem)
            if beh_state.rem(n) == 1
                scatter3(delta_power(n),theta_power(n),movement(n),30,[0 0 1]);
            elseif beh_state.nrem(n) == 1
                scatter3(delta_power(n),theta_power(n),movement(n),30,[1 0 0]);
            else
                scatter3(delta_power(n),theta_power(n),movement(n),30,[0 0 0]);
            end
        end
        xlabel('delta')
        ylabel('theta')
        zlabel('movement')
        grid on

    subplot(4,2,2)
        hold on;
        m1_spindle = [];
        dls_spindle = [];
        for n_event = 1:length(sleep_rhythms.spindles{1,1}.pks)
            m1_spindle = [m1_spindle; m1_lfp(round((sleep_rhythms.spindles{1,1}.pks(n_event)*lfp_samp_rate)-lfp_samp_rate):round((sleep_rhythms.spindles{1,1}.pks(n_event)*lfp_samp_rate)+lfp_samp_rate))];
            dls_spindle = [dls_spindle; dls_lfp(round((sleep_rhythms.spindles{1,1}.pks(n_event)*lfp_samp_rate)-lfp_samp_rate):round((sleep_rhythms.spindles{1,1}.pks(n_event)*lfp_samp_rate)+lfp_samp_rate))];
        end
        plot(mean(m1_spindle,1),'color','r','LineWidth',2);
        plot(mean(dls_spindle,1),'color','b','LineWidth',2);
        title('spindles');
        xlim([1 2035])

    subplot(4,2,4)
        hold on;
        m1_SO = [];
        dls_SO = [];
        for n_event = 1:length(sleep_rhythms.so_delta.SO_up_states)
            m1_SO = [m1_SO; m1_lfp(round((sleep_rhythms.so_delta.SO_up_states(n_event)*lfp_samp_rate)-lfp_samp_rate):round((sleep_rhythms.so_delta.SO_up_states(n_event)*lfp_samp_rate)+lfp_samp_rate))];
            dls_SO = [dls_SO; dls_lfp(round((sleep_rhythms.so_delta.SO_up_states(n_event)*lfp_samp_rate)-lfp_samp_rate):round((sleep_rhythms.so_delta.SO_up_states(n_event)*lfp_samp_rate)+lfp_samp_rate))];
        end
        plot(mean(m1_SO,1),'color','r','LineWidth',2);
        plot(mean(dls_SO,1),'color','b','LineWidth',2);
        title('slow oscillations');
        xlim([1 2035])

    subplot(4,2,6)
        hold on;
        m1_spindle = [];
        dls_spindle = [];
        for n_event = 1:length(sleep_rhythms.so_delta.delta_up_states)
            m1_spindle = [m1_spindle; m1_lfp(round((sleep_rhythms.so_delta.delta_up_states(n_event)*lfp_samp_rate)-lfp_samp_rate):round((sleep_rhythms.so_delta.delta_up_states(n_event)*lfp_samp_rate)+lfp_samp_rate))];
            dls_spindle = [dls_spindle; dls_lfp(round((sleep_rhythms.so_delta.delta_up_states(n_event)*lfp_samp_rate)-lfp_samp_rate):round((sleep_rhythms.so_delta.delta_up_states(n_event)*lfp_samp_rate)+lfp_samp_rate))];
        end
        plot(mean(m1_spindle,1),'color','r','LineWidth',2);
        plot(mean(dls_spindle,1),'color','b','LineWidth',2);
        title('delta');
        xlim([1 2035])

    subplot(4,2,8)
        hold on;
        tbin = [1:lfp_samp_rate*60:length(beh_state.nrem_interp)];
        spindle_rate = [];
        so_rate = [];
        delta_rate = [];
        for n_bin = 1:length(tbin)-1
            if min(beh_state.nrem_interp(tbin(n_bin):tbin(n_bin+1)))==0; continue; end
            spindle_rate = [spindle_rate; sum(sleep_rhythms.spindles{1,1}.pks*lfp_samp_rate>tbin(n_bin) & sleep_rhythms.spindles{1,1}.pks*lfp_samp_rate<tbin(n_bin+1))];
            so_rate = [so_rate; sum(sleep_rhythms.so_delta.SO_up_states*lfp_samp_rate>tbin(n_bin) & sleep_rhythms.so_delta.SO_up_states*lfp_samp_rate<tbin(n_bin+1))];
            delta_rate = [delta_rate; sum(sleep_rhythms.so_delta.delta_up_states*lfp_samp_rate>tbin(n_bin) & sleep_rhythms.so_delta.delta_up_states*lfp_samp_rate<tbin(n_bin+1))];
        end
        plot(zscore(smooth(spindle_rate,10)),'color','r','LineWidth',1);
        plot(zscore(smooth(so_rate,10)),'color','b','LineWidth',1);
        plot(zscore(smooth(delta_rate,10)),'color','g','LineWidth',1);
        title('event rate')

%% SAVE
   
    if ~isempty(varargin{3})
        
        save_path = varargin{3};
        if ~exist(save_path)
            mkdir(save_path)
        end
        cd(save_path);
        saveas(fig,['day_' num2str(day) '_block_' num2str(block) '_sleep_activity.png'])

        mean_lfp.lfp_samp_rate = lfp_samp_rate;
        mean_lfp.m1_lfp = m1_lfp;
        mean_lfp.dls_lfp = dls_lfp;
        save(['day_' num2str(day) '_block_' num2str(block) '_sleep_activity.mat'],'mean_lfp','sleep_rhythms','beh_state','-v7.3')

    end
    
end