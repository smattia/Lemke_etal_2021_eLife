function out = detect_slow_oscillations(lfp,Fs,sleep_idx)
%
%   Detects slow oscillations and delta waves in LFP signal - contributions
%   from Daniel Silversmith, Nikhilesh Natraj, and Jaekyung Kim
%
%   identifies times of down states during sleep block
%       1. averages activity across channels
%       2. filters average lfp signal for delta activity (.1-4Hz)
%       3. finds positive-to-negative zero crossings (zc) during sleep
%       4. for each pos-to-neg zc
%           4.1. get the prev peak (down state)
%           4.2. get the prev negative-to-positive zc
%           4.3. get the next negative-to-positive zc
%       5. returns slow oscillations that fit the following criteria
%           5.1. time btw neg-to-pos zcs is btw .8 and 2s
%           5.2. prev peak exceeds a threshold of 90% of all peaks
%       
%   Inputs:
%       lfp - lfp is a matrix (samples x channels) of lfp data from sleep
%       Fs - sample rate of lfp 
% 
%   Optional Name/Value Pairs:
%       'sleep_idx' - a logical vector 
%                   - (1-sleep,0-awake) corresponding to that session's lfp
%       'artifact_idx' - a logical vector 
%                   - (1-artifact,0-fine) corresponding to that session's lfp
%       'bad_channels' - vector of channels that are ignored (default = [])
%       'PLOT' - [0|1]
%       'sleep_classify' - [0|1], only for classified sleep (sleep_idx)
%       'mnl_parm' - [peak-thr trough-thr dur-min dur-max], manual parameter setting
% 
%   Outputs:
%       out
%           .down_states - time of down states
%           .zc - time of pos-to-neg down state
%           .peak - peak of z-scored delta band lfp signal
%           .trough - trough of z-scored delta band lfp signal

%%% PARAMETERS
    
    mnl_parm=[85 40 .15 .5]; 
    % SO = 85% treshold previous peak 
    %      40% threshold previous trough
    %      greater than .15 sec between peak and trough
    %      less than .5 sec between peak and trough 
    %
    % DELTA = bottom 85% previous peak 
    %         40% threshold previous trough
    %         greater than .10 sec between peak and trough
    
%%% ORGANIZE INPUTS 

    dt = 1/Fs;
    time = ((1:size(lfp,1))/Fs)' - dt;

%%% FILTER LFP FOR DELTA

    fpass = [.1,4];

    % highpass
    [b,a] = butter(2,fpass(1)/(Fs/2),'high');
    lfp_delta = filtfilt(b,a,lfp);

    % lowpass
    [b,a] = butter(4,fpass(2)/(Fs/2),'low');%original from Daniel
    %     [b,a] = butter(2,fpass(2)/(Fs/2),'low');%JK 04/18/19
    lfp_delta = filtfilt(b,a,lfp_delta);

%%% FIND POS-TO-NEG ZERO CROSSINGS

    idx = round(Fs):(size(lfp,1)-round(Fs));
    ptnzc = round(Fs)-1 + find(lfp_delta(idx)>=0 & lfp_delta(idx+1)<0 & sleep_idx(idx)==1);

%%% GET STATS ON ALL POS-TO-NEG ZERO CROSSINGS

    prev_peak_idx = zeros(size(ptnzc));
    prev_peak     = zeros(size(ptnzc));
    next_trough   = zeros(size(ptnzc));
    dur           = zeros(size(ptnzc));
    isasleep      = zeros(size(ptnzc));

    for i=2:length(ptnzc)-1

        prev_ntpzc          = ptnzc(i-1)-1 + find(lfp_delta(ptnzc(i-1):ptnzc(i  ))<0,1,'last');
        next_ntpzc          = ptnzc(i)  -1 + find(lfp_delta(ptnzc(i  ):ptnzc(i+1))<0,1,'last');
        [prev_peak(i),idx]  = max(lfp_delta(ptnzc(i-1):ptnzc(i)));
        prev_peak_idx(i)    = ptnzc(i-1)-1 + idx;
        [next_trough(i),idx]= min(lfp_delta(ptnzc(i):next_ntpzc));
        next_trough_idx(i)  = ptnzc(i)-1 + idx;
        dur(i)              = time(next_trough_idx(i) - prev_peak_idx(i));
        isasleep(i)         = sleep_idx(prev_peak_idx(i));

    end

%%% APPLY CRITERIA FOR SLOW OSCILLATIONS

    thresh_up = prctile(prev_peak,mnl_parm(1));
    thresh_dwn = prctile(next_trough,mnl_parm(2));
    thresh = [thresh_up,thresh_dwn];

    idx_SO = isasleep & prev_peak>=thresh_up & next_trough<thresh_dwn & dur>mnl_parm(3) & dur<mnl_parm(4);
    idx_DELTA = isasleep & prev_peak<thresh_up & next_trough<thresh_dwn & dur>(mnl_parm(3)-.05);

%%% OUTPUT

    out.SO_down_states = time(prev_peak_idx(idx_SO));
    out.SO_up_states = time(next_trough_idx(idx_SO));
    out.SO_zc = time(ptnzc(idx_SO));
    out.SO_peaks = prev_peak(idx_SO);
    out.SO_troughs = next_trough(idx_SO);
    out.SO_dur = dur(idx_SO);    

    out.delta_down_states = time(prev_peak_idx(idx_DELTA));
    out.delta_up_states = time(next_trough_idx(idx_DELTA));
    out.delta_zc = time(ptnzc(idx_DELTA));
    out.delta_peaks = prev_peak(idx_DELTA);
    out.delta_troughs = next_trough(idx_DELTA);
    out.delta_dur = dur(idx_DELTA);

end
