function spindles = detect_spindles(lfp,Fs,sleep_idx)
%
%   Detects spindles in LFP signal - contributions from Daniel 
%   Silversmith, Nikhilesh Natraj, and Jaekyung Kim
%
%   INPUTS:
%       lfp - cell array, each cell is a session containing
%           - a vector of lfp data (samples x 1) in sleep1
%
%   OPTIONAL NAME/VALUE PAIRS:
%       'Fs' - sampling rate of LFP (default is 2.44140625 / 24)
%       'sleep_idx' - cell array, each cell contains a logical vector 
%                   - (1-sleep,0-awake) corresponding to that session's lfp
%       'artifact_idx' - cell array, each cell contains a logical vector 
%                   - (1-artifact,0-fine) corresponding to that session's lfp
%       'fpass' - for filtering lfp signal (default = [10,14])
%       'PLOT' - [0|1]
%       'single_thresh' - [0|1], if 1 uses one threshold for all sessions
%       'avg_flag' - [0|1], if 1 use average lfp across channels
%       'bad_channels' - vector of bad channels (ignored)
% 
%   OUTPUTS:
%       spindles - output structure for sleep1 (defined below)
% 
%       spindles.
%           pks - time peak of filtered lfp signal for each thresh crossing
%           start - upward thresh crossing
%           finish - downward thresh crossing
%           amp - peak amplitude of envelope
%           dur - duration (secs)
%

%% INPUTS
% Fs = 2.44140625 / 24;
% sleep_idx = cell(1,3);
fpass = [10,14];
PLOT = 0;
artifact_idx = cell(1,3);
single_thresh = 0;
avg_flag = 0;
bad_channels = [];
% assignopts(who,varargin);
DEBUG = false;

%% SIZING INFO
NumSessions = length(lfp);
for s=1:NumSessions,
    N = length(lfp{s});
    dt = 1/Fs;
    time{s} = ((1:N)/Fs)' - dt;
end

%% REMOVE ARTIFACT FROM LFP 
for s=1:NumSessions,
    lfp{s}(artifact_idx{s}) = mean(lfp{s}(~artifact_idx{s}));
end

%% AVERAGE LFP SIGNAL
if avg_flag,
    for s=1:NumSessions,
        good = setdiff(1:size(lfp{s},2),bad_channels);
        lfp{s} = mean(zscore(lfp{s}(:,good)),2);
    end
end

%% GET SPINDLE POWER

% FILTER SIGNAL
[bhigh,ahigh] = butter(6,fpass(1)/(Fs/2),'high');
[blow,alow] = butter(10,fpass(2)/(Fs/2),'low');
for s=1:NumSessions,
    flfp{s} = filtfilt(bhigh,ahigh,lfp{s});
    flfp{s} = filtfilt(blow,alow,flfp{s});
end

%% ZSCORE LFP TO AWAKE ACTIVITY
for s=1:NumSessions,
%     artifact = zeros(size(sleep_idx{s}));
%     artifact(artifact_idx{s}) = 1;
%     artifact = artifact==1;
%     zlfp{s} = myzscore(flfp{s},sleep_idx{s}(:));
%     zlfp{s} = myzscore(flfp{s},[sleep_idx{s}(:)|artifact(:)]);
    zlfp{s} = flfp{s};
end

%% GET SIGNAL ENVELOPE
gwin = gausswin(round(200e-3*Fs));
gwin = gwin / sum(gwin);
for s=1:NumSessions,
    H{s} = abs(hilbert(zlfp{s}));
    H{s} = conv(H{s},gwin,'same'); % smooth w/ 100ms gaussian kernel 
%     H{s} = teager(zlfp{s});
    sleep_H{s} = nan(size(H{s}));
%     sleep_H{s}(sleep_idx{s}==1,:) = H{s}(sleep_idx{s}==1,:);
    sleep_H{s}(sleep_idx{s}==1) = H{s}(sleep_idx{s}==1);
end

%% SET THRESHOLD
if single_thresh,
    total_sleep_H = [];
    for s=1:NumSessions,
%         total_sleep_H = vertcat(1,total_sleep_H,H{s});
        total_sleep_H = vertcat(1,total_sleep_H,zlfp{s});
    end
    sigma = std(total_sleep_H);
    mu = mean(total_sleep_H);
    low = mu + 1.5*sigma; % 1.5
    high = mu + 2.5*sigma; % 2.5
    % can also try percentiles
%     low = prctile(total_sleep_H,85); 
%     high = prctile(total_sleep_H,95); 
    low_thresh = repmat(low,1,3);
    high_thresh = repmat(high,1,3);
else
    for s=1:NumSessions,
%         total_sleep = H{s};
        total_sleep = zlfp{s};
        sigma = std(total_sleep);
        mu = mean(total_sleep);
        low = mu + 1.5*sigma; % 1.5
        high = mu + 2.5*sigma; % 2.5
        % can also try percentiles
%         low = prctile(H{s},85); 
%         high = prctile(H{s},95); 
        low_thresh(s) = low;
        high_thresh(s) = high;
    end
end

%% FIND SPINDLES
for s=1:NumSessions,

    % EACH SESSION
    mag = sleep_H{s};
    low = low_thresh(s);
    high = high_thresh(s);
    
    % LIMIT SEARCH
    % only search within 1s of start and finish of recording,
    % also make sure signal is below threshold at boundaries
    idx0 = round(1*Fs);
    idx0 = idx0 + find(mag(idx0:end)<low,1,'first')-1;
    idxf = length(mag) - round(1*Fs)+1;
    idxf = idx0 + find(mag(idx0:idxf)<low,1,'last')-1;
    idx = idx0:idxf;

    % THRESHOLD CROSSINGS
    upidx  = idx0 + find(mag(idx)<high & mag(idx+1)>=high)-1;
    startidx = [];
    finishidx = [];
    for i=1:length(upidx),
        startidx(i) = find(H{s}(1:upidx(i))<low,1,'last');
        finishidx(i) = upidx(i)-1 + find(H{s}(upidx(i):end)<low,1,'first');
    end
    start = time{s}(startidx);
    finish = time{s}(finishidx);

    % COMBINE SPINDLES < 1 CYCLE (1/10Hz ~ 100ms)
    [startidx,finishidx,start,finish] = ...
        combine_spindles(startidx,finishidx,start,finish,.3);
    
    % REMOVE SPINDLES WITH DURATION < 400MS (< 3 CYCLES)
    dur = finish - start;
    idx = dur>.5 & dur<2.5;
    startidx = startidx(idx);
    finishidx = finishidx(idx);
    dur = dur(idx);
    
    % FIND TIME OF PEAK OF FILTERED LFP SIGNAL
    pksidx = zeros(size(startidx));
    for i=1:length(startidx),
        idx0 = startidx(i);
        idxf = finishidx(i);
        pksidx(i) = idx0 + find(flfp{s}(idx0:idxf)==max(flfp{s}(idx0:idxf)))-1;
    end

%     % DO MPP ALG AND FIND PEAK FREQ
%     for i=1:length(startidx),
%         idx = pksidx(i) + round(-2*Fs:2*Fs);
%         N = length(idx);
%         Z = zeros(1,2^nextpow2(N)-N);
%         sig = [zlfp{s}(idx);Z'];
%         [pwr,~,f] = mpp(sig,Fs);
%         fidx = f>6 & f<20;
%         spec = mean(pwr(fidx,:),2);
%         f = f(fidx);
%         maxfreq(i) = f(spec==max(spec));
% %         subplot(211)
% %         plot(lfp{s}(idx))
% %         subplot(212)
% %         plot(f(fidx),spec)
% %         waitforbuttonpress
%     end   
    
    % FIND PEAK OF ENVELOPE FILTERED LFP SIGNAL
    amp = zeros(size(startidx));
    for i=1:length(startidx),
        idx0 = startidx(i);
        idxf = finishidx(i);
        amp(i) = max(H{s}(idx0:idxf));
    end

    % CONVERT FROM INDICES TO TIME
    spindles{s}.pks = time{s}(pksidx);
    spindles{s}.start  = time{s}(startidx);
    spindles{s}.finish = time{s}(finishidx);
    spindles{s}.dur = dur;
    spindles{s}.amp = amp;
%     spindles{s}.pkFreq = maxfreq;
    spindlesidx{s} = pksidx;
end % session

%% plot
if PLOT,
    for s=1:NumSessions,
        low = low_thresh(s);
        high = high_thresh(s);
        
        figure(s)
        cc = get(gca,'ColorOrder');
        set(gcf,'Name',sprintf('Spindle_Oscillations'))
        set(gcf,'Position',[680,620,980,360])
        subplot(2,3,[1,4]), hold on
        [W,t] = triggered_lfp(lfp{s},Fs,spindles{s}.pks,[1,1]);
        shadedErrorBar(t,mean(W{1},2),std(W{1},[],2),{'-','Color',cc(1,:)},1)
        xlim([t(1),t(end)])
        hline([low,high])
        
        ax(1)=subplot(2,3,2:3); hold on
        sig = mean(lfp{s},2);
        plot(time{s},sig,'Color',cc(1,:))
        plot(time{s}(sleep_idx{s}),repmat(1*max(sig),sum(sleep_idx{s}),1),'.','Color',cc(2,:))
        
        ax(2)=subplot(2,3,5:6); hold on
        plot(time{s},zlfp{s},'Color',cc(1,:))
        plot(spindles{s}.pks,zlfp{s}(spindlesidx{s}),'k*')
        hline([low,high])
        
        linkaxes(ax,'x')
    end
end % plot
end % detect_spindles

function zlfp = myzscore(lfp,sleep_idx)
    % sleep1
    mu = mean(lfp(~sleep_idx,:));
    sigma = std(lfp(~sleep_idx,:));
    zlfp = bsxfun(@minus,lfp,mu);
    zlfp = bsxfun(@rdivide,zlfp,sigma);
    zlfp = mean(zlfp,2);
end


