function [all_pairs, all_m1_spindle, all_dls_spindle, all_m1_SO, all_dls_SO, all_m1_delta, all_dls_delta] = CC_by_sleep_rhythm(pre_m1_spiking, pre_dls_spiking, post_m1_spiking, post_dls_spiking, pre_lfp, post_lfp, pre_video_start_stop, post_video_start_stop,pre_sleep_rhythms,post_sleep_rhythms);
%
    
%% FILTER LFP 

    [bhigh,ahigh] = butter(6,10/(pre_lfp.lfp_samp_rate/2),'high');
    [blow,alow] = butter(10,14/(pre_lfp.lfp_samp_rate/2),'low');

    pre_m1_spindle_lfp = filtfilt(bhigh,ahigh,pre_lfp.m1_lfp);
    pre_m1_spindle_lfp = filtfilt(blow,alow,pre_m1_spindle_lfp);
    post_m1_spindle_lfp = filtfilt(bhigh,ahigh,post_lfp.m1_lfp);
    post_m1_spindle_lfp = filtfilt(blow,alow,post_m1_spindle_lfp);

    if size(pre_m1_spindle_lfp,2)==1
        pre_m1_spindle_lfp = pre_m1_spindle_lfp';
    end
    if size(post_m1_spindle_lfp,2)==1
        post_m1_spindle_lfp = post_m1_spindle_lfp';
    end

    pre_dls_spindle_lfp = filtfilt(bhigh,ahigh,pre_lfp.dls_lfp);
    pre_dls_spindle_lfp = filtfilt(blow,alow,pre_dls_spindle_lfp);
    post_dls_spindle_lfp = filtfilt(bhigh,ahigh,post_lfp.dls_lfp);
    post_dls_spindle_lfp = filtfilt(blow,alow,post_dls_spindle_lfp);

    if size(pre_dls_spindle_lfp,2)==1
        pre_dls_spindle_lfp = pre_dls_spindle_lfp';
    end
    if size(post_dls_spindle_lfp,2)==1
        post_dls_spindle_lfp = post_dls_spindle_lfp';
    end

    [bhigh,ahigh] = butter(2,.1/(pre_lfp.lfp_samp_rate/2),'high');
    [blow,alow] = butter(4,4/(pre_lfp.lfp_samp_rate/2),'low');

    pre_m1_SO_lfp = filtfilt(bhigh,ahigh,pre_lfp.m1_lfp);
    pre_m1_SO_lfp = filtfilt(blow,alow,pre_m1_SO_lfp);
    post_m1_SO_lfp = filtfilt(bhigh,ahigh,post_lfp.m1_lfp);
    post_m1_SO_lfp = filtfilt(blow,alow,post_m1_SO_lfp);

    if size(pre_m1_SO_lfp,2)==1
        pre_m1_SO_lfp = pre_m1_SO_lfp';
    end
    if size(post_m1_SO_lfp,2)==1
        post_m1_SO_lfp = post_m1_SO_lfp';
    end

    pre_dls_SO_lfp = filtfilt(bhigh,ahigh,pre_lfp.dls_lfp);
    pre_dls_SO_lfp = filtfilt(blow,alow,pre_dls_SO_lfp);
    post_dls_SO_lfp = filtfilt(bhigh,ahigh,post_lfp.dls_lfp);
    post_dls_SO_lfp = filtfilt(blow,alow,post_dls_SO_lfp);

    if size(pre_dls_SO_lfp,2)==1
        pre_dls_SO_lfp = pre_dls_SO_lfp';
    end
    if size(post_dls_SO_lfp,2)==1
        post_dls_SO_lfp = post_dls_SO_lfp';
    end
    
%% M1 RHYTHM MODULATION & PHASE LOCKING
    
    unit_count = 1;
    
    for m1_chan = 1:size(pre_m1_spiking,1)
        for m1_unit = 1:size(pre_m1_spiking,2)
            if isempty(pre_m1_spiking{m1_chan,m1_unit}) || isempty(post_m1_spiking{m1_chan,m1_unit})
                continue
            end

            all_m1_spindle(unit_count).m1_chan_unit = [m1_chan m1_unit];
            all_m1_SO(unit_count).m1_chan_unit = [m1_chan m1_unit];
            all_m1_delta(unit_count).m1_chan_unit = [m1_chan m1_unit];
            
            %%% PRE %%%

                pre_m1_unit = pre_m1_spiking{m1_chan,m1_unit}(pre_m1_spiking{m1_chan,m1_unit}>pre_video_start_stop(1) & pre_m1_spiking{m1_chan,m1_unit}<pre_video_start_stop(2))-pre_video_start_stop(1);

                %%% SPINDLES

                    spindle_raster = [];
                    spindle_ms_raster = [];
                    spindle_mod = [];
                    spindle_phase = [];
                    spindle_lfp = [];

                    control_raster = [];
                    control_mod = [];
                    control_phase = [];
                    control_lfp = [];

                    for spindle = 1:length(pre_sleep_rhythms.spindles{1,1}.pks)

                        control_pk =pre_sleep_rhythms.spindles{1,1}.pks(spindle)+randi([-10 -5],1,1);
                        if control_pk<=0.5
                            continue
                        end
                           
                        spindle_lfp = [spindle_lfp; pre_m1_spindle_lfp(((pre_sleep_rhythms.spindles{1,1}.pks(spindle)-0.5)*pre_lfp.lfp_samp_rate):((pre_sleep_rhythms.spindles{1,1}.pks(spindle)+0.5)*pre_lfp.lfp_samp_rate))];
                        control_lfp = [control_lfp; pre_m1_spindle_lfp((control_pk-0.5)*pre_lfp.lfp_samp_rate:(control_pk+0.5)*pre_lfp.lfp_samp_rate)];

                        spindle_raster = [spindle_raster; histcounts(pre_m1_unit((pre_m1_unit>=pre_sleep_rhythms.spindles{1,1}.pks(spindle)-1) & (pre_m1_unit<=pre_sleep_rhythms.spindles{1,1}.pks(spindle)+1))-pre_sleep_rhythms.spindles{1,1}.pks(spindle),[-1:0.01:1])];
                        spindle_ms_raster = [spindle_ms_raster; histcounts(pre_m1_unit((pre_m1_unit>=pre_sleep_rhythms.spindles{1,1}.pks(spindle)-1) & (pre_m1_unit<=pre_sleep_rhythms.spindles{1,1}.pks(spindle)+1))-pre_sleep_rhythms.spindles{1,1}.pks(spindle),[-1:0.001:1])];
                        control_raster = [control_raster; histcounts(pre_m1_unit((pre_m1_unit>=control_pk-1) & (pre_m1_unit<=control_pk+1))-control_pk,[-1:0.01:1])];                

                        spindle_mod = [spindle_mod; histcounts(pre_m1_unit((pre_m1_unit>=pre_sleep_rhythms.spindles{1,1}.pks(spindle)-1) & (pre_m1_unit<=pre_sleep_rhythms.spindles{1,1}.pks(spindle)+1))-pre_sleep_rhythms.spindles{1,1}.pks(spindle),[-1:0.04:1])];
                        control_mod = [control_mod; histcounts(pre_m1_unit((pre_m1_unit>=control_pk-1) & (pre_m1_unit<=control_pk+1))-control_pk,[-1:0.04:1])];                

                        tmp_phase = angle(hilbert(pre_m1_spindle_lfp(((pre_sleep_rhythms.spindles{1,1}.pks(spindle)-0.5)*pre_lfp.lfp_samp_rate):((pre_sleep_rhythms.spindles{1,1}.pks(spindle)+0.5)*pre_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(pre_m1_unit((pre_m1_unit>=pre_sleep_rhythms.spindles{1,1}.pks(spindle)-0.5) & (pre_m1_unit<=pre_sleep_rhythms.spindles{1,1}.pks(spindle)+0.5))-pre_sleep_rhythms.spindles{1,1}.pks(spindle),[-0.5:1/pre_lfp.lfp_samp_rate:.5]);
                        spindle_phase = [spindle_phase tmp_phase(tmp_spikes==1)];

                        tmp_phase = angle(hilbert(pre_m1_spindle_lfp(((control_pk-0.5)*pre_lfp.lfp_samp_rate):((control_pk+0.5)*pre_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(pre_m1_unit((pre_m1_unit>=control_pk-0.5) & (pre_m1_unit<=control_pk+0.5))-control_pk,[-0.5:1/pre_lfp.lfp_samp_rate:.5]);
                        control_phase = [control_phase tmp_phase(tmp_spikes==1)];

                    end

                    %%% PLV

                        spindle_plv = 0;
                        for n=1:length(spindle_phase)
                            spindle_plv = spindle_plv + exp(1i*spindle_phase(n));
                        end
                        spindle_plv_angle = angle(spindle_plv);
                        spindle_plv_mag = abs(spindle_plv)/length(spindle_phase);
                        
                        all_m1_spindle(unit_count).pre_m1_phase_spindle = spindle_phase; 
                        all_m1_spindle(unit_count).pre_m1_plv_spindle = [spindle_plv_angle spindle_plv_mag]; 

                        control_plv = 0;
                        for n=1:length(control_phase)
                            control_plv = control_plv + exp(1i*control_phase(n));
                        end
                        control_plv_angle = angle(control_plv);
                        control_plv_mag = abs(control_plv)/length(control_phase);
                        
                        all_m1_spindle(unit_count).pre_m1_phase_spindle_control = control_phase; 
                        all_m1_spindle(unit_count).pre_m1_plv_spindle_control = [control_plv_angle control_plv_mag]; 

                    %%% MODULATION

                        all_m1_spindle(unit_count).pre_m1_spindle_maxFR = max(mean(spindle_mod(:,13:38)))*25;
                        all_m1_spindle(unit_count).pre_m1_spindle_minFR = min(mean(spindle_mod(:,13:38)))*25;

                        all_m1_spindle(unit_count).pre_m1_control_maxFR = max(mean(control_mod(:,13:38)))*25;
                        all_m1_spindle(unit_count).pre_m1_control_minFR = min(mean(control_mod(:,13:38)))*25;

                        all_m1_spindle(unit_count).pre_m1_spindle_raster = spindle_raster;
                        all_m1_spindle(unit_count).pre_m1_spindle_ms_raster = spindle_ms_raster;
                        all_m1_spindle(unit_count).pre_m1_spindle_control_raster = control_raster;

                        all_m1_spindle(unit_count).pre_m1_spindle_lfp = spindle_lfp;
                        all_m1_spindle(unit_count).pre_m1_spindle_control_lfp = control_lfp;
                        
%                     figure; 
% 
%                         subplot(4,2,[1 2]);
%                             hold on;
%                             plot(mean(spindle_raster),'color','r')
%                             plot(mean(control_raster),'color','k')
%                             xlim([1 200])
%                             title('spindle vs. control raster')
% 
%                         subplot(4,2,[3 4]);
%                             hold on;
%                             plot(mean(spindle_lfp),'r')
%                             plot(mean(control_lfp),'k')
%                             xlim([1 round(pre_lfp.lfp_samp_rate)])
% 
%                         subplot(4,2,5);
%                             hold on;
%                             plot(mean(spindle_mod),'color','k')
%                             [~, max_i] = max(mean(spindle_mod(:,13:38)));
%                             [~, min_i] = min(mean(spindle_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['spindle modulation: ' num2str(all_m1_spindle(unit_count).pre_m1_spindle_maxFR-all_m1_spindle(unit_count).pre_m1_spindle_minFR)])
% 
%                         subplot(4,2,6);
%                            hold on;
%                             plot(mean(control_mod),'color','k')
%                             [~, max_i] = max(mean(control_mod(:,13:38)));
%                             [~, min_i] = min(mean(control_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['control modulation: ' num2str(all_m1_spindle(unit_count).pre_m1_control_maxFR-all_m1_spindle(unit_count).pre_m1_control_minFR)])
% 
%                         subplot(4,2,7);
%                             polarhistogram(spindle_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 spindle_plv_angle],[0 spindle_plv_mag],'LineWidth',3)
%                             tmp_rlim = rlim;
% 
%                         subplot(4,2,8);
%                             polarhistogram(control_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 control_plv_angle],[0 control_plv_mag],'LineWidth',3)
%                             rlim(tmp_rlim)

                %%% SO

                    SO_raster = [];
                    SO_mod = [];
                    SO_phase = [];
                    SO_lfp = [];

                    control_raster = [];
                    control_mod = [];
                    control_phase = [];
                    control_lfp = [];

                    for SO = 1:length(pre_sleep_rhythms.so_delta.SO_up_states)

                        control_pk =pre_sleep_rhythms.so_delta.SO_up_states(SO)+randi([-10 -5],1,1);
                        if control_pk<=0.5
                            continue;
                        end
                        
                        SO_lfp = [SO_lfp; pre_m1_SO_lfp(((pre_sleep_rhythms.so_delta.SO_up_states(SO)-0.5)*pre_lfp.lfp_samp_rate):((pre_sleep_rhythms.so_delta.SO_up_states(SO)+0.5)*pre_lfp.lfp_samp_rate))];
                        control_lfp = [control_lfp; pre_m1_SO_lfp(((control_pk-0.5)*pre_lfp.lfp_samp_rate):((control_pk+0.5)*pre_lfp.lfp_samp_rate))];

                        SO_raster = [SO_raster; histcounts(pre_m1_unit((pre_m1_unit>=pre_sleep_rhythms.so_delta.SO_up_states(SO)-1) & (pre_m1_unit<=pre_sleep_rhythms.so_delta.SO_up_states(SO)+1))-pre_sleep_rhythms.so_delta.SO_up_states(SO),[-1:0.01:1])];
                        control_raster = [control_raster; histcounts(pre_m1_unit((pre_m1_unit>=control_pk-1) & (pre_m1_unit<=control_pk+1))-control_pk,[-1:0.01:1])];                

                        SO_mod = [SO_mod; histcounts(pre_m1_unit((pre_m1_unit>=pre_sleep_rhythms.so_delta.SO_up_states(SO)-1) & (pre_m1_unit<=pre_sleep_rhythms.so_delta.SO_up_states(SO)+1))-pre_sleep_rhythms.so_delta.SO_up_states(SO),[-1:0.04:1])];
                        control_mod = [control_mod; histcounts(pre_m1_unit((pre_m1_unit>=control_pk-1) & (pre_m1_unit<=control_pk+1))-control_pk,[-1:0.04:1])];                

                        tmp_phase = angle(hilbert(pre_m1_SO_lfp(((pre_sleep_rhythms.so_delta.SO_up_states(SO)-0.5)*pre_lfp.lfp_samp_rate):((pre_sleep_rhythms.so_delta.SO_up_states(SO)+0.5)*pre_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(pre_m1_unit((pre_m1_unit>=pre_sleep_rhythms.so_delta.SO_up_states(SO)-0.5) & (pre_m1_unit<=pre_sleep_rhythms.so_delta.SO_up_states(SO)+0.5))-pre_sleep_rhythms.so_delta.SO_up_states(SO),[-0.5:1/pre_lfp.lfp_samp_rate:.5]);
                        SO_phase = [SO_phase tmp_phase(tmp_spikes==1)];

                        tmp_phase = angle(hilbert(pre_m1_SO_lfp(((control_pk-0.5)*pre_lfp.lfp_samp_rate):((control_pk+0.5)*pre_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(pre_m1_unit((pre_m1_unit>=control_pk-0.5) & (pre_m1_unit<=control_pk+0.5))-control_pk,[-0.5:1/pre_lfp.lfp_samp_rate:.5]);
                        control_phase = [control_phase tmp_phase(tmp_spikes==1)];

                    end

                    %%% PLV

                        SO_plv = 0;
                        for n=1:length(SO_phase)
                            SO_plv = SO_plv + exp(1i*SO_phase(n));
                        end
                        SO_plv_angle = angle(SO_plv);
                        SO_plv_mag = abs(SO_plv)/length(SO_phase);
                        
                        all_m1_SO(unit_count).pre_m1_phase_SO = SO_phase; 
                        all_m1_SO(unit_count).pre_m1_plv_SO = [SO_plv_angle SO_plv_mag]; 

                        control_plv = 0;
                        for n=1:length(control_phase)
                            control_plv = control_plv + exp(1i*control_phase(n));
                        end
                        control_plv_angle = angle(control_plv);
                        control_plv_mag = abs(control_plv)/length(control_phase);
                        
                        all_m1_SO(unit_count).pre_m1_phase_SO_control = control_phase; 
                        all_m1_SO(unit_count).pre_m1_plv_SO_control = [control_plv_angle control_plv_mag]; 

                    %%% MODULATION

                        all_m1_SO(unit_count).pre_m1_SO_maxFR = max(mean(SO_mod(:,13:38)))*25;
                        all_m1_SO(unit_count).pre_m1_SO_minFR = min(mean(SO_mod(:,13:38)))*25;

                        all_m1_SO(unit_count).pre_m1_control_maxFR = max(mean(control_mod(:,13:38)))*25;
                        all_m1_SO(unit_count).pre_m1_control_minFR = min(mean(control_mod(:,13:38)))*25;

                        all_m1_SO(unit_count).pre_m1_SO_raster = SO_raster;
                        all_m1_SO(unit_count).pre_m1_SO_control_raster = control_raster;

                        all_m1_SO(unit_count).pre_m1_SO_lfp = SO_lfp;
                        all_m1_SO(unit_count).pre_m1_SO_control_lfp = control_lfp;
                        
%                     figure; 
% 
%                         subplot(4,2,[1 2]);
%                             hold on;
%                             plot(mean(SO_raster),'color','r')
%                             plot(mean(control_raster),'color','k')
%                             xlim([1 200])
%                             title('SO vs. control raster')
% 
%                         subplot(4,2,[3 4]);
%                             hold on;
%                             plot(mean(SO_lfp),'r')
%                             plot(mean(control_lfp),'k')
%                             xlim([1 round(pre_lfp.lfp_samp_rate)])
% 
%                         subplot(4,2,5);
%                             hold on;
%                             plot(mean(SO_mod),'color','k')
%                             [~, max_i] = max(mean(SO_mod(:,13:38)));
%                             [~, min_i] = min(mean(SO_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['SO modulation: ' num2str(all_m1_SO(unit_count).pre_m1_SO_maxFR-all_m1_SO(unit_count).pre_m1_SO_minFR)])
% 
%                         subplot(4,2,6);
%                            hold on;
%                             plot(mean(control_mod),'color','k')
%                             [~, max_i] = max(mean(control_mod(:,13:38)));
%                             [~, min_i] = min(mean(control_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['control modulation: ' num2str(all_m1_SO(unit_count).pre_m1_control_maxFR-all_m1_SO(unit_count).pre_m1_control_minFR)])
% 
%                         subplot(4,2,7);
%                             polarhistogram(SO_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 SO_plv_angle],[0 SO_plv_mag],'LineWidth',3)
%                             tmp_rlim = rlim;
% 
%                         subplot(4,2,8);
%                             polarhistogram(control_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 control_plv_angle],[0 control_plv_mag],'LineWidth',3)
%                             rlim(tmp_rlim)

                %%% delta

                    delta_raster = [];
                    delta_mod = [];
                    delta_phase = [];
                    delta_lfp = [];

                    control_raster = [];
                    control_mod = [];
                    control_phase = [];
                    control_lfp = [];

                    for delta = 10:length(pre_sleep_rhythms.so_delta.delta_up_states)

                        control_pk =pre_sleep_rhythms.so_delta.delta_up_states(delta)+randi([-10 -5],1,1);
                        if control_pk<=0.5
                            continue;
                        end
                        
                        delta_lfp = [delta_lfp; pre_m1_SO_lfp(((pre_sleep_rhythms.so_delta.delta_up_states(delta)-0.5)*pre_lfp.lfp_samp_rate):((pre_sleep_rhythms.so_delta.delta_up_states(delta)+0.5)*pre_lfp.lfp_samp_rate))];
                        control_lfp = [control_lfp; pre_m1_SO_lfp(((control_pk-0.5)*pre_lfp.lfp_samp_rate):((control_pk+0.5)*pre_lfp.lfp_samp_rate))];

                        delta_raster = [delta_raster; histcounts(pre_m1_unit((pre_m1_unit>=pre_sleep_rhythms.so_delta.delta_up_states(delta)-1) & (pre_m1_unit<=pre_sleep_rhythms.so_delta.delta_up_states(delta)+1))-pre_sleep_rhythms.so_delta.delta_up_states(delta),[-1:0.01:1])];
                        control_raster = [control_raster; histcounts(pre_m1_unit((pre_m1_unit>=control_pk-1) & (pre_m1_unit<=control_pk+1))-control_pk,[-1:0.01:1])];                

                        delta_mod = [delta_mod; histcounts(pre_m1_unit((pre_m1_unit>=pre_sleep_rhythms.so_delta.delta_up_states(delta)-1) & (pre_m1_unit<=pre_sleep_rhythms.so_delta.delta_up_states(delta)+1))-pre_sleep_rhythms.so_delta.delta_up_states(delta),[-1:0.04:1])];
                        control_mod = [control_mod; histcounts(pre_m1_unit((pre_m1_unit>=control_pk-1) & (pre_m1_unit<=control_pk+1))-control_pk,[-1:0.04:1])];                

                        tmp_phase = angle(hilbert(pre_m1_SO_lfp(((pre_sleep_rhythms.so_delta.delta_up_states(delta)-0.5)*pre_lfp.lfp_samp_rate):((pre_sleep_rhythms.so_delta.delta_up_states(delta)+0.5)*pre_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(pre_m1_unit((pre_m1_unit>=pre_sleep_rhythms.so_delta.delta_up_states(delta)-0.5) & (pre_m1_unit<=pre_sleep_rhythms.so_delta.delta_up_states(delta)+0.5))-pre_sleep_rhythms.so_delta.delta_up_states(delta),[-0.5:1/pre_lfp.lfp_samp_rate:.5]);
                        delta_phase = [delta_phase tmp_phase(tmp_spikes==1)];

                        tmp_phase = angle(hilbert(pre_m1_SO_lfp(((control_pk-0.5)*pre_lfp.lfp_samp_rate):((control_pk+0.5)*pre_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(pre_m1_unit((pre_m1_unit>=control_pk-0.5) & (pre_m1_unit<=control_pk+0.5))-control_pk,[-0.5:1/pre_lfp.lfp_samp_rate:.5]);
                        control_phase = [control_phase tmp_phase(tmp_spikes==1)];

                    end

                    %%% PLV

                        delta_plv = 0;
                        for n=1:length(delta_phase)
                            delta_plv = delta_plv + exp(1i*delta_phase(n));
                        end
                        delta_plv_angle = angle(delta_plv);
                        delta_plv_mag = abs(delta_plv)/length(delta_phase);
                        
                        all_m1_delta(unit_count).pre_m1_phase_delta = delta_phase; 
                        all_m1_delta(unit_count).pre_m1_plv_delta = [delta_plv_angle delta_plv_mag]; 

                        control_plv = 0;
                        for n=1:length(control_phase)
                            control_plv = control_plv + exp(1i*control_phase(n));
                        end
                        control_plv_angle = angle(control_plv);
                        control_plv_mag = abs(control_plv)/length(control_phase);
                        
                        all_m1_delta(unit_count).pre_m1_phase_delta_control = control_phase; 
                        all_m1_delta(unit_count).pre_m1_plv_delta_control = [control_plv_angle control_plv_mag]; 

                    %%% MODULATION

                        all_m1_delta(unit_count).pre_m1_delta_maxFR = max(mean(delta_mod(:,13:38)))*25;
                        all_m1_delta(unit_count).pre_m1_delta_minFR = min(mean(delta_mod(:,13:38)))*25;

                        all_m1_delta(unit_count).pre_m1_control_maxFR = max(mean(control_mod(:,13:38)))*25;
                        all_m1_delta(unit_count).pre_m1_control_minFR = min(mean(control_mod(:,13:38)))*25;

                        all_m1_delta(unit_count).pre_m1_delta_raster = delta_raster;
                        all_m1_delta(unit_count).pre_m1_delta_control_raster = control_raster;

                        all_m1_delta(unit_count).pre_m1_delta_lfp = delta_lfp;
                        all_m1_delta(unit_count).pre_m1_delta_control_lfp = control_lfp;
                        
%                     figure; 
% 
%                         subplot(4,2,[1 2]);
%                             hold on;
%                             plot(mean(delta_raster),'color','r')
%                             plot(mean(control_raster),'color','k')
%                             xlim([1 200])
%                             title('delta vs. control raster')
% 
%                         subplot(4,2,[3 4]);
%                             hold on;
%                             plot(mean(delta_lfp),'r')
%                             plot(mean(control_lfp),'k')
%                             xlim([1 round(pre_lfp.lfp_samp_rate)])
% 
%                         subplot(4,2,5);
%                             hold on;
%                             plot(mean(delta_mod),'color','k')
%                             [~, max_i] = max(mean(delta_mod(:,13:38)));
%                             [~, min_i] = min(mean(delta_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['delta modulation: ' num2str(all_m1_delta(unit_count).pre_m1_delta_maxFR-all_m1_delta(unit_count).pre_m1_delta_minFR)])
% 
%                         subplot(4,2,6);
%                            hold on;
%                             plot(mean(control_mod),'color','k')
%                             [~, max_i] = max(mean(control_mod(:,13:38)));
%                             [~, min_i] = min(mean(control_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['control modulation: ' num2str(all_m1_delta(unit_count).pre_m1_control_maxFR-all_m1_delta(unit_count).pre_m1_control_minFR)])
% 
%                         subplot(4,2,7);
%                             polarhistogram(delta_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 delta_plv_angle],[0 delta_plv_mag],'LineWidth',3)
%                             tmp_rlim = rlim;
% 
%                         subplot(4,2,8);
%                             polarhistogram(control_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 control_plv_angle],[0 control_plv_mag],'LineWidth',3)
%                             rlim(tmp_rlim)

            %%% POST %%%

                post_m1_unit = post_m1_spiking{m1_chan,m1_unit}(post_m1_spiking{m1_chan,m1_unit}>post_video_start_stop(1) & post_m1_spiking{m1_chan,m1_unit}<post_video_start_stop(2))-post_video_start_stop(1);

                %%% SPINDLES

                    spindle_raster = [];
                    spindle_ms_raster = [];
                    spindle_mod = [];
                    spindle_phase = [];
                    spindle_lfp = [];

                    control_raster = [];
                    control_mod = [];
                    control_phase = [];
                    control_lfp = [];

                    for spindle = 1:length(post_sleep_rhythms.spindles{1,1}.pks)

                        control_pk =post_sleep_rhythms.spindles{1,1}.pks(spindle)+randi([-10 -5],1,1);
                        if control_pk<=0.5
                            continue;
                        end
                        
                        spindle_lfp = [spindle_lfp; post_m1_spindle_lfp(((post_sleep_rhythms.spindles{1,1}.pks(spindle)-0.5)*post_lfp.lfp_samp_rate):((post_sleep_rhythms.spindles{1,1}.pks(spindle)+0.5)*post_lfp.lfp_samp_rate))];
                        control_lfp = [control_lfp; post_m1_spindle_lfp(((control_pk-0.5)*post_lfp.lfp_samp_rate):((control_pk+0.5)*post_lfp.lfp_samp_rate))];

                        spindle_raster = [spindle_raster; histcounts(post_m1_unit((post_m1_unit>=post_sleep_rhythms.spindles{1,1}.pks(spindle)-1) & (post_m1_unit<=post_sleep_rhythms.spindles{1,1}.pks(spindle)+1))-post_sleep_rhythms.spindles{1,1}.pks(spindle),[-1:0.01:1])];
                        spindle_ms_raster = [spindle_ms_raster; histcounts(post_m1_unit((post_m1_unit>=post_sleep_rhythms.spindles{1,1}.pks(spindle)-1) & (post_m1_unit<=post_sleep_rhythms.spindles{1,1}.pks(spindle)+1))-post_sleep_rhythms.spindles{1,1}.pks(spindle),[-1:0.001:1])];
                        control_raster = [control_raster; histcounts(post_m1_unit((post_m1_unit>=control_pk-1) & (post_m1_unit<=control_pk+1))-control_pk,[-1:0.01:1])];                

                        spindle_mod = [spindle_mod; histcounts(post_m1_unit((post_m1_unit>=post_sleep_rhythms.spindles{1,1}.pks(spindle)-1) & (post_m1_unit<=post_sleep_rhythms.spindles{1,1}.pks(spindle)+1))-post_sleep_rhythms.spindles{1,1}.pks(spindle),[-1:0.04:1])];
                        control_mod = [control_mod; histcounts(post_m1_unit((post_m1_unit>=control_pk-1) & (post_m1_unit<=control_pk+1))-control_pk,[-1:0.04:1])];                

                        tmp_phase = angle(hilbert(post_m1_spindle_lfp(((post_sleep_rhythms.spindles{1,1}.pks(spindle)-0.5)*post_lfp.lfp_samp_rate):((post_sleep_rhythms.spindles{1,1}.pks(spindle)+0.5)*post_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(post_m1_unit((post_m1_unit>=post_sleep_rhythms.spindles{1,1}.pks(spindle)-0.5) & (post_m1_unit<=post_sleep_rhythms.spindles{1,1}.pks(spindle)+0.5))-post_sleep_rhythms.spindles{1,1}.pks(spindle),[-0.5:1/post_lfp.lfp_samp_rate:.5]);
                        spindle_phase = [spindle_phase tmp_phase(tmp_spikes==1)];

                        tmp_phase = angle(hilbert(post_m1_spindle_lfp(((control_pk-0.5)*post_lfp.lfp_samp_rate):((control_pk+0.5)*post_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(post_m1_unit((post_m1_unit>=control_pk-0.5) & (post_m1_unit<=control_pk+0.5))-control_pk,[-0.5:1/post_lfp.lfp_samp_rate:.5]);
                        control_phase = [control_phase tmp_phase(tmp_spikes==1)];

                    end

                    %%% PLV

                        spindle_plv = 0;
                        for n=1:length(spindle_phase)
                            spindle_plv = spindle_plv + exp(1i*spindle_phase(n));
                        end
                        spindle_plv_angle = angle(spindle_plv);
                        spindle_plv_mag = abs(spindle_plv)/length(spindle_phase);
                        
                        all_m1_spindle(unit_count).post_m1_phase_spindle = spindle_phase; 
                        all_m1_spindle(unit_count).post_m1_plv_spindle = [spindle_plv_angle spindle_plv_mag]; 

                        control_plv = 0;
                        for n=1:length(control_phase)
                            control_plv = control_plv + exp(1i*control_phase(n));
                        end
                        control_plv_angle = angle(control_plv);
                        control_plv_mag = abs(control_plv)/length(control_phase);
                        
                        all_m1_spindle(unit_count).post_m1_phase_spindle_control = control_phase; 
                        all_m1_spindle(unit_count).post_m1_plv_spindle_control = [control_plv_angle control_plv_mag]; 

                    %%% MODULATION

                        all_m1_spindle(unit_count).post_m1_spindle_maxFR = max(mean(spindle_mod(:,13:38)))*25;
                        all_m1_spindle(unit_count).post_m1_spindle_minFR = min(mean(spindle_mod(:,13:38)))*25;

                        all_m1_spindle(unit_count).post_m1_control_maxFR = max(mean(control_mod(:,13:38)))*25;
                        all_m1_spindle(unit_count).post_m1_control_minFR = min(mean(control_mod(:,13:38)))*25;

                        all_m1_spindle(unit_count).post_m1_spindle_raster = spindle_raster;
                        all_m1_spindle(unit_count).post_m1_spindle_ms_raster = spindle_ms_raster;
                        all_m1_spindle(unit_count).post_m1_spindle_control_raster = control_raster;

                        all_m1_spindle(unit_count).post_m1_spindle_lfp = spindle_lfp;
                        all_m1_spindle(unit_count).post_m1_spindle_control_lfp = control_lfp;
                        
%                     figure; 
% 
%                         subplot(4,2,[1 2]);
%                             hold on;
%                             plot(mean(spindle_raster),'color','r')
%                             plot(mean(control_raster),'color','k')
%                             xlim([1 200])
%                             title('spindle vs. control raster')
% 
%                         subplot(4,2,[3 4]);
%                             hold on;
%                             plot(mean(spindle_lfp),'r')
%                             plot(mean(control_lfp),'k')
%                             xlim([1 round(post_lfp.lfp_samp_rate)])
% 
%                         subplot(4,2,5);
%                             hold on;
%                             plot(mean(spindle_mod),'color','k')
%                             [~, max_i] = max(mean(spindle_mod(:,13:38)));
%                             [~, min_i] = min(mean(spindle_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['spindle modulation: ' num2str(all_m1_spindle(unit_count).post_m1_spindle_maxFR-all_m1_spindle(unit_count).post_m1_spindle_minFR)])
% 
%                         subplot(4,2,6);
%                            hold on;
%                             plot(mean(control_mod),'color','k')
%                             [~, max_i] = max(mean(control_mod(:,13:38)));
%                             [~, min_i] = min(mean(control_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['control modulation: ' num2str(all_m1_spindle(unit_count).post_m1_control_maxFR-all_m1_spindle(unit_count).post_m1_control_minFR)])
% 
%                         subplot(4,2,7);
%                             polarhistogram(spindle_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 spindle_plv_angle],[0 spindle_plv_mag],'LineWidth',3)
%                             tmp_rlim = rlim;
% 
%                         subplot(4,2,8);
%                             polarhistogram(control_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 control_plv_angle],[0 control_plv_mag],'LineWidth',3)
%                             rlim(tmp_rlim)

                %%% SO

                    SO_raster = [];
                    SO_mod = [];
                    SO_phase = [];
                    SO_lfp = [];

                    control_raster = [];
                    control_mod = [];
                    control_phase = [];
                    control_lfp = [];

                    for SO = 1:length(post_sleep_rhythms.so_delta.SO_up_states)

                        control_pk =post_sleep_rhythms.so_delta.SO_up_states(SO)+randi([-10 -5],1,1);
                        if control_pk<=0.5
                            continue;
                        end
                        
                        SO_lfp = [SO_lfp; post_m1_SO_lfp(((post_sleep_rhythms.so_delta.SO_up_states(SO)-0.5)*post_lfp.lfp_samp_rate):((post_sleep_rhythms.so_delta.SO_up_states(SO)+0.5)*post_lfp.lfp_samp_rate))];
                        control_lfp = [control_lfp; post_m1_SO_lfp(((control_pk-0.5)*post_lfp.lfp_samp_rate):((control_pk+0.5)*post_lfp.lfp_samp_rate))];

                        SO_raster = [SO_raster; histcounts(post_m1_unit((post_m1_unit>=post_sleep_rhythms.so_delta.SO_up_states(SO)-1) & (post_m1_unit<=post_sleep_rhythms.so_delta.SO_up_states(SO)+1))-post_sleep_rhythms.so_delta.SO_up_states(SO),[-1:0.01:1])];
                        control_raster = [control_raster; histcounts(post_m1_unit((post_m1_unit>=control_pk-1) & (post_m1_unit<=control_pk+1))-control_pk,[-1:0.01:1])];                

                        SO_mod = [SO_mod; histcounts(post_m1_unit((post_m1_unit>=post_sleep_rhythms.so_delta.SO_up_states(SO)-1) & (post_m1_unit<=post_sleep_rhythms.so_delta.SO_up_states(SO)+1))-post_sleep_rhythms.so_delta.SO_up_states(SO),[-1:0.04:1])];
                        control_mod = [control_mod; histcounts(post_m1_unit((post_m1_unit>=control_pk-1) & (post_m1_unit<=control_pk+1))-control_pk,[-1:0.04:1])];                

                        tmp_phase = angle(hilbert(post_m1_SO_lfp(((post_sleep_rhythms.so_delta.SO_up_states(SO)-0.5)*post_lfp.lfp_samp_rate):((post_sleep_rhythms.so_delta.SO_up_states(SO)+0.5)*post_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(post_m1_unit((post_m1_unit>=post_sleep_rhythms.so_delta.SO_up_states(SO)-0.5) & (post_m1_unit<=post_sleep_rhythms.so_delta.SO_up_states(SO)+0.5))-post_sleep_rhythms.so_delta.SO_up_states(SO),[-0.5:1/post_lfp.lfp_samp_rate:.5]);
                        SO_phase = [SO_phase tmp_phase(tmp_spikes==1)];

                        tmp_phase = angle(hilbert(post_m1_SO_lfp(((control_pk-0.5)*post_lfp.lfp_samp_rate):((control_pk+0.5)*post_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(post_m1_unit((post_m1_unit>=control_pk-0.5) & (post_m1_unit<=control_pk+0.5))-control_pk,[-0.5:1/post_lfp.lfp_samp_rate:.5]);
                        control_phase = [control_phase tmp_phase(tmp_spikes==1)];

                    end

                    %%% PLV

                        SO_plv = 0;
                        for n=1:length(SO_phase)
                            SO_plv = SO_plv + exp(1i*SO_phase(n));
                        end
                        SO_plv_angle = angle(SO_plv);
                        SO_plv_mag = abs(SO_plv)/length(SO_phase);
                        
                        all_m1_SO(unit_count).post_m1_phase_SO = SO_phase;
                        all_m1_SO(unit_count).post_m1_plv_SO = [SO_plv_angle SO_plv_mag]; 

                        control_plv = 0;
                        for n=1:length(control_phase)
                            control_plv = control_plv + exp(1i*control_phase(n));
                        end
                        control_plv_angle = angle(control_plv);
                        control_plv_mag = abs(control_plv)/length(control_phase);
                        
                        all_m1_SO(unit_count).post_m1_phase_SO_control = control_phase; 
                        all_m1_SO(unit_count).post_m1_plv_SO_control = [control_plv_angle control_plv_mag]; 

                    %%% MODULATION

                        all_m1_SO(unit_count).post_m1_SO_maxFR = max(mean(SO_mod(:,13:38)))*25;
                        all_m1_SO(unit_count).post_m1_SO_minFR = min(mean(SO_mod(:,13:38)))*25;

                        all_m1_SO(unit_count).post_m1_control_maxFR = max(mean(control_mod(:,13:38)))*25;
                        all_m1_SO(unit_count).post_m1_control_minFR = min(mean(control_mod(:,13:38)))*25;

                        all_m1_SO(unit_count).post_m1_SO_raster = SO_raster;
                        all_m1_SO(unit_count).post_m1_SO_control_raster = control_raster;

                        all_m1_SO(unit_count).post_m1_SO_lfp = SO_lfp;
                        all_m1_SO(unit_count).post_m1_SO_control_lfp = control_lfp;
                        
%                     figure; 
% 
%                         subplot(4,2,[1 2]);
%                             hold on;
%                             plot(mean(SO_raster),'color','r')
%                             plot(mean(control_raster),'color','k')
%                             xlim([1 200])
%                             title('SO vs. control raster')
% 
%                         subplot(4,2,[3 4]);
%                             hold on;
%                             plot(mean(SO_lfp),'r')
%                             plot(mean(control_lfp),'k')
%                             xlim([1 round(post_lfp.lfp_samp_rate)])
% 
%                         subplot(4,2,5);
%                             hold on;
%                             plot(mean(SO_mod),'color','k')
%                             [~, max_i] = max(mean(SO_mod(:,13:38)));
%                             [~, min_i] = min(mean(SO_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['SO modulation: ' num2str(all_m1_SO(unit_count).post_m1_SO_maxFR-all_m1_SO(unit_count).post_m1_SO_minFR)])
% 
%                         subplot(4,2,6);
%                            hold on;
%                             plot(mean(control_mod),'color','k')
%                             [~, max_i] = max(mean(control_mod(:,13:38)));
%                             [~, min_i] = min(mean(control_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['control modulation: ' num2str(all_m1_SO(unit_count).post_m1_control_maxFR-all_m1_SO(unit_count).post_m1_control_minFR)])
% 
%                         subplot(4,2,7);
%                             polarhistogram(SO_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 SO_plv_angle],[0 SO_plv_mag],'LineWidth',3)
%                             tmp_rlim = rlim;
% 
%                         subplot(4,2,8);
%                             polarhistogram(control_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 control_plv_angle],[0 control_plv_mag],'LineWidth',3)
%                             rlim(tmp_rlim)

                %%% delta

                    delta_raster = [];
                    delta_mod = [];
                    delta_phase = [];
                    delta_lfp = [];

                    control_raster = [];
                    control_mod = [];
                    control_phase = [];
                    control_lfp = [];

                    for delta = 1:length(post_sleep_rhythms.so_delta.delta_up_states)

                        control_pk =post_sleep_rhythms.so_delta.delta_up_states(delta)+randi([-10 -5],1,1);
                        if control_pk<=0.5
                            continue;
                        end
                        
                        delta_lfp = [delta_lfp; post_m1_SO_lfp(((post_sleep_rhythms.so_delta.delta_up_states(delta)-0.5)*post_lfp.lfp_samp_rate):((post_sleep_rhythms.so_delta.delta_up_states(delta)+0.5)*post_lfp.lfp_samp_rate))];
                        control_lfp = [control_lfp; post_m1_SO_lfp(((control_pk-0.5)*post_lfp.lfp_samp_rate):((control_pk+0.5)*post_lfp.lfp_samp_rate))];

                        delta_raster = [delta_raster; histcounts(post_m1_unit((post_m1_unit>=post_sleep_rhythms.so_delta.delta_up_states(delta)-1) & (post_m1_unit<=post_sleep_rhythms.so_delta.delta_up_states(delta)+1))-post_sleep_rhythms.so_delta.delta_up_states(delta),[-1:0.01:1])];
                        control_raster = [control_raster; histcounts(post_m1_unit((post_m1_unit>=control_pk-1) & (post_m1_unit<=control_pk+1))-control_pk,[-1:0.01:1])];                

                        delta_mod = [delta_mod; histcounts(post_m1_unit((post_m1_unit>=post_sleep_rhythms.so_delta.delta_up_states(delta)-1) & (post_m1_unit<=post_sleep_rhythms.so_delta.delta_up_states(delta)+1))-post_sleep_rhythms.so_delta.delta_up_states(delta),[-1:0.04:1])];
                        control_mod = [control_mod; histcounts(post_m1_unit((post_m1_unit>=control_pk-1) & (post_m1_unit<=control_pk+1))-control_pk,[-1:0.04:1])];                

                        tmp_phase = angle(hilbert(post_m1_SO_lfp(((post_sleep_rhythms.so_delta.delta_up_states(delta)-0.5)*post_lfp.lfp_samp_rate):((post_sleep_rhythms.so_delta.delta_up_states(delta)+0.5)*post_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(post_m1_unit((post_m1_unit>=post_sleep_rhythms.so_delta.delta_up_states(delta)-0.5) & (post_m1_unit<=post_sleep_rhythms.so_delta.delta_up_states(delta)+0.5))-post_sleep_rhythms.so_delta.delta_up_states(delta),[-0.5:1/post_lfp.lfp_samp_rate:.5]);
                        delta_phase = [delta_phase tmp_phase(tmp_spikes==1)];

                        tmp_phase = angle(hilbert(post_m1_SO_lfp(((control_pk-0.5)*post_lfp.lfp_samp_rate):((control_pk+0.5)*post_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(post_m1_unit((post_m1_unit>=control_pk-0.5) & (post_m1_unit<=control_pk+0.5))-control_pk,[-0.5:1/post_lfp.lfp_samp_rate:.5]);
                        control_phase = [control_phase tmp_phase(tmp_spikes==1)];

                    end

                    %%% PLV

                        delta_plv = 0;
                        for n=1:length(delta_phase)
                            delta_plv = delta_plv + exp(1i*delta_phase(n));
                        end
                        delta_plv_angle = angle(delta_plv);
                        delta_plv_mag = abs(delta_plv)/length(delta_phase);
                        
                        all_m1_delta(unit_count).post_m1_phase_delta = delta_phase; 
                        all_m1_delta(unit_count).post_m1_plv_delta = [delta_plv_angle delta_plv_mag]; 

                        control_plv = 0;
                        for n=1:length(control_phase)
                            control_plv = control_plv + exp(1i*control_phase(n));
                        end
                        control_plv_angle = angle(control_plv);
                        control_plv_mag = abs(control_plv)/length(control_phase);
                        
                        all_m1_delta(unit_count).post_m1_phase_delta_control = control_phase; 
                        all_m1_delta(unit_count).post_m1_plv_delta_control = [control_plv_angle control_plv_mag]; 

                    %%% MODULATION

                        all_m1_delta(unit_count).post_m1_delta_maxFR = max(mean(delta_mod(:,13:38)))*25;
                        all_m1_delta(unit_count).post_m1_delta_minFR = min(mean(delta_mod(:,13:38)))*25;

                        all_m1_delta(unit_count).post_m1_control_maxFR = max(mean(control_mod(:,13:38)))*25;
                        all_m1_delta(unit_count).post_m1_control_minFR = min(mean(control_mod(:,13:38)))*25;

                        all_m1_delta(unit_count).post_m1_delta_raster = delta_raster;
                        all_m1_delta(unit_count).post_m1_delta_control_raster = control_raster;

                        all_m1_delta(unit_count).post_m1_delta_lfp = delta_lfp;
                        all_m1_delta(unit_count).post_m1_delta_control_lfp = control_lfp;
                        
%                     figure; 
% 
%                         subplot(4,2,[1 2]);
%                             hold on;
%                             plot(mean(delta_raster),'color','r')
%                             plot(mean(control_raster),'color','k')
%                             xlim([1 200])
%                             title('delta vs. control raster')
% 
%                         subplot(4,2,[3 4]);
%                             hold on;
%                             plot(mean(delta_lfp),'r')
%                             plot(mean(control_lfp),'k')
%                             xlim([1 round(post_lfp.lfp_samp_rate)])
% 
%                         subplot(4,2,5);
%                             hold on;
%                             plot(mean(delta_mod),'color','k')
%                             [~, max_i] = max(mean(delta_mod(:,13:38)));
%                             [~, min_i] = min(mean(delta_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['delta modulation: ' num2str(all_m1_delta(unit_count).post_m1_delta_maxFR-all_m1_delta(unit_count).post_m1_delta_minFR)])
% 
%                         subplot(4,2,6);
%                            hold on;
%                             plot(mean(control_mod),'color','k')
%                             [~, max_i] = max(mean(control_mod(:,13:38)));
%                             [~, min_i] = min(mean(control_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['control modulation: ' num2str(all_m1_delta(unit_count).post_m1_control_maxFR-all_m1_delta(unit_count).post_m1_control_minFR)])
% 
%                         subplot(4,2,7);
%                             polarhistogram(delta_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 delta_plv_angle],[0 delta_plv_mag],'LineWidth',3)
%                             tmp_rlim = rlim;
% 
%                         subplot(4,2,8);
%                             polarhistogram(control_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 control_plv_angle],[0 control_plv_mag],'LineWidth',3)
%                             rlim(tmp_rlim)

            unit_count = unit_count + 1;
            pause(0.1);
            
        end
    end
    
    %%% SPINDLE PLOT
    
        figure;

            subplot(3,2,1)
            
                spindle_tot = [0 0];
                control_tot = [0 0];
                for n=1:length(all_m1_spindle)
                    [x y] = pol2cart(all_m1_spindle(n).pre_m1_plv_spindle(1), all_m1_spindle(n).pre_m1_plv_spindle(2));
                    if isnan(x) || isnan(y); continue; end
                    spindle_tot = spindle_tot + [x y];
                    [x y] = pol2cart(all_m1_spindle(n).pre_m1_plv_spindle_control(1), all_m1_spindle(n).pre_m1_plv_spindle_control(2));
                    if isnan(x) || isnan(y); continue; end
                    control_tot = control_tot + [x y];
                end

                polarplot([0 all_m1_spindle(1).pre_m1_plv_spindle(1)],[0 all_m1_spindle(1).pre_m1_plv_spindle(2)],'color',[1 0.5 0.5],'LineWidth',1); 
                hold on;
                polarplot([0 all_m1_spindle(1).pre_m1_plv_spindle_control(1)],[0 all_m1_spindle(1).pre_m1_plv_spindle_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                for n=2:length(all_m1_spindle)
                    polarplot([0 all_m1_spindle(n).pre_m1_plv_spindle(1)],[0 all_m1_spindle(n).pre_m1_plv_spindle(2)],'color',[1 0.5 0.5],'LineWidth',1);
                    polarplot([0 all_m1_spindle(n).pre_m1_plv_spindle_control(1)],[0 all_m1_spindle(n).pre_m1_plv_spindle_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                end
                
                [theta rho] = cart2pol(spindle_tot(1), spindle_tot(2));
                polarplot([0 theta],[0 rho/length(all_m1_spindle)],'color','r','LineWidth',3);
                
                [theta rho] = cart2pol(control_tot(1), control_tot(2));
                polarplot([0 theta],[0 rho/length(all_m1_spindle)],'color','k','LineWidth',3);

            subplot(3,2,2)
            
                spindle_tot = [0 0];
                control_tot = [0 0];
                for n=1:length(all_m1_spindle)
                    [x y] = pol2cart(all_m1_spindle(n).post_m1_plv_spindle(1), all_m1_spindle(n).post_m1_plv_spindle(2));
                    if isnan(x) || isnan(y); continue; end
                    spindle_tot = spindle_tot + [x y];
                    [x y] = pol2cart(all_m1_spindle(n).post_m1_plv_spindle_control(1), all_m1_spindle(n).post_m1_plv_spindle_control(2));
                    if isnan(x) || isnan(y); continue; end
                    control_tot = control_tot + [x y];
                end

                polarplot([0 all_m1_spindle(1).post_m1_plv_spindle(1)],[0 all_m1_spindle(1).post_m1_plv_spindle(2)],'color',[1 0.5 0.5],'LineWidth',1); 
                hold on;
                polarplot([0 all_m1_spindle(1).post_m1_plv_spindle_control(1)],[0 all_m1_spindle(1).post_m1_plv_spindle_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                for n=2:length(all_m1_spindle)
                    polarplot([0 all_m1_spindle(n).post_m1_plv_spindle(1)],[0 all_m1_spindle(n).post_m1_plv_spindle(2)],'color',[1 0.5 0.5],'LineWidth',1);
                    polarplot([0 all_m1_spindle(n).post_m1_plv_spindle_control(1)],[0 all_m1_spindle(n).post_m1_plv_spindle_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                end
                
                [theta rho] = cart2pol(spindle_tot(1), spindle_tot(2));
                polarplot([0 theta],[0 rho/length(all_m1_spindle)],'color','r','LineWidth',3);
                
                [theta rho] = cart2pol(control_tot(1), control_tot(2));
                polarplot([0 theta],[0 rho/length(all_m1_spindle)],'color','k','LineWidth',3);


            subplot(3,2,3)
                hold on;
                tmp_mean_spindle = [];
                tmp_mean_control = [];
                for n=1:length(all_m1_spindle)
                    tmp_mean_spindle = [tmp_mean_spindle; mean(all_m1_spindle(n).pre_m1_spindle_raster)];
                    tmp_mean_control = [tmp_mean_control; mean(all_m1_spindle(n).pre_m1_spindle_control_raster)];
                end
                plot(mean(tmp_mean_spindle),'color',[1 0 0])
                plot(mean(tmp_mean_control),'color',[0.5 0.5 0.5])

            subplot(3,2,4)
                hold on;
                tmp_mean_spindle = [];
                tmp_mean_control = [];
                for n=1:length(all_m1_spindle)
                    tmp_mean_spindle = [tmp_mean_spindle; mean(all_m1_spindle(n).post_m1_spindle_raster)];
                    tmp_mean_control = [tmp_mean_control; mean(all_m1_spindle(n).post_m1_spindle_control_raster)];
                end
                plot(mean(tmp_mean_spindle),'color',[1 0 0])
                plot(mean(tmp_mean_control),'color',[0.5 0.5 0.5])

            subplot(3,2,5)
                hold on;
                tmp_spindle_mod = [];
                tmp_control_mod = [];
                for n=1:length(all_m1_spindle)
                    scatter(1,all_m1_spindle(n).pre_m1_spindle_maxFR-all_m1_spindle(n).pre_m1_spindle_minFR,[30],[1 0 0]);
                    scatter(2,all_m1_spindle(n).pre_m1_control_maxFR-all_m1_spindle(n).pre_m1_control_minFR,[30],[0.5 0.5 0.5]);
                    tmp_spindle_mod = [tmp_spindle_mod all_m1_spindle(n).pre_m1_spindle_maxFR-all_m1_spindle(n).pre_m1_spindle_minFR];
                    tmp_control_mod = [tmp_control_mod all_m1_spindle(n).pre_m1_control_maxFR-all_m1_spindle(n).pre_m1_control_minFR];
                end
                errorbar(1,mean(tmp_spindle_mod),std(tmp_spindle_mod)/sqrt(length(tmp_spindle_mod)),'color',[1 0 0],'LineWidth',2,'CapSize',20);
                errorbar(2,mean(tmp_control_mod),std(tmp_control_mod)/sqrt(length(tmp_control_mod)),'color',[0 0 0],'LineWidth',2,'CapSize',20);
                xlim([0 3])

            subplot(3,2,6)
                hold on;
                tmp_spindle_mod = [];
                tmp_control_mod = [];
                for n=1:length(all_m1_spindle)
                    scatter(1,all_m1_spindle(n).post_m1_spindle_maxFR-all_m1_spindle(n).post_m1_spindle_minFR,[30],[1 0 0]);
                    scatter(2,all_m1_spindle(n).post_m1_control_maxFR-all_m1_spindle(n).post_m1_control_minFR,[30],[0.5 0.5 0.5]);
                    tmp_spindle_mod = [tmp_spindle_mod all_m1_spindle(n).post_m1_spindle_maxFR-all_m1_spindle(n).post_m1_spindle_minFR];
                    tmp_control_mod = [tmp_control_mod all_m1_spindle(n).post_m1_control_maxFR-all_m1_spindle(n).post_m1_control_minFR];
                end
                errorbar(1,mean(tmp_spindle_mod),std(tmp_spindle_mod)/sqrt(length(tmp_spindle_mod)),'color',[1 0 0],'LineWidth',2,'CapSize',20);
                errorbar(2,mean(tmp_control_mod),std(tmp_control_mod)/sqrt(length(tmp_control_mod)),'color',[0 0 0],'LineWidth',2,'CapSize',20);
                xlim([0 3])  
            
    %%% SO PLOT
    
        figure;

            subplot(3,2,1)
            
                SO_tot = [0 0];
                control_tot = [0 0];
                for n=1:length(all_m1_SO)
                    [x y] = pol2cart(all_m1_SO(n).pre_m1_plv_SO(1), all_m1_SO(n).pre_m1_plv_SO(2));
                    if isnan(x) || isnan(y); continue; end
                    SO_tot = SO_tot + [x y];
                    [x y] = pol2cart(all_m1_SO(n).pre_m1_plv_SO_control(1), all_m1_SO(n).pre_m1_plv_SO_control(2));
                    if isnan(x) || isnan(y); continue; end
                    control_tot = control_tot + [x y];
                end

                polarplot([0 all_m1_SO(1).pre_m1_plv_SO(1)],[0 all_m1_SO(1).pre_m1_plv_SO(2)],'color',[1 0.5 0.5],'LineWidth',1); 
                hold on;
                polarplot([0 all_m1_SO(1).pre_m1_plv_SO_control(1)],[0 all_m1_SO(1).pre_m1_plv_SO_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                for n=2:length(all_m1_SO)
                    polarplot([0 all_m1_SO(n).pre_m1_plv_SO(1)],[0 all_m1_SO(n).pre_m1_plv_SO(2)],'color',[1 0.5 0.5],'LineWidth',1);
                    polarplot([0 all_m1_SO(n).pre_m1_plv_SO_control(1)],[0 all_m1_SO(n).pre_m1_plv_SO_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                end
                
                [theta rho] = cart2pol(SO_tot(1), SO_tot(2));
                polarplot([0 theta],[0 rho/length(all_m1_SO)],'color','r','LineWidth',3);
                
                [theta rho] = cart2pol(control_tot(1), control_tot(2));
                polarplot([0 theta],[0 rho/length(all_m1_SO)],'color','k','LineWidth',3);

            subplot(3,2,2)
            
                SO_tot = [0 0];
                control_tot = [0 0];
                for n=1:length(all_m1_SO)
                    [x y] = pol2cart(all_m1_SO(n).post_m1_plv_SO(1), all_m1_SO(n).post_m1_plv_SO(2));
                    if isnan(x) || isnan(y); continue; end
                    SO_tot = SO_tot + [x y];
                    [x y] = pol2cart(all_m1_SO(n).post_m1_plv_SO_control(1), all_m1_SO(n).post_m1_plv_SO_control(2));
                    if isnan(x) || isnan(y); continue; end
                    control_tot = control_tot + [x y];
                end

                polarplot([0 all_m1_SO(1).post_m1_plv_SO(1)],[0 all_m1_SO(1).post_m1_plv_SO(2)],'color',[1 0.5 0.5],'LineWidth',1); 
                hold on;
                polarplot([0 all_m1_SO(1).post_m1_plv_SO_control(1)],[0 all_m1_SO(1).post_m1_plv_SO_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                for n=2:length(all_m1_SO)
                    polarplot([0 all_m1_SO(n).post_m1_plv_SO(1)],[0 all_m1_SO(n).post_m1_plv_SO(2)],'color',[1 0.5 0.5],'LineWidth',1);
                    polarplot([0 all_m1_SO(n).post_m1_plv_SO_control(1)],[0 all_m1_SO(n).post_m1_plv_SO_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                end
                
                [theta rho] = cart2pol(SO_tot(1), SO_tot(2));
                polarplot([0 theta],[0 rho/length(all_m1_SO)],'color','r','LineWidth',3);
                
                [theta rho] = cart2pol(control_tot(1), control_tot(2));
                polarplot([0 theta],[0 rho/length(all_m1_SO)],'color','k','LineWidth',3);


            subplot(3,2,3)
                hold on;
                tmp_mean_SO = [];
                tmp_mean_control = [];
                for n=1:length(all_m1_SO)
                    tmp_mean_SO = [tmp_mean_SO; mean(all_m1_SO(n).pre_m1_SO_raster)];
                    tmp_mean_control = [tmp_mean_control; mean(all_m1_SO(n).pre_m1_SO_control_raster)];
                end
                plot(mean(tmp_mean_SO),'color',[1 0 0])
                plot(mean(tmp_mean_control),'color',[0.5 0.5 0.5])

            subplot(3,2,4)
                hold on;
                tmp_mean_SO = [];
                tmp_mean_control = [];
                for n=1:length(all_m1_SO)
                    tmp_mean_SO = [tmp_mean_SO; mean(all_m1_SO(n).post_m1_SO_raster)];
                    tmp_mean_control = [tmp_mean_control; mean(all_m1_SO(n).post_m1_SO_control_raster)];
                end
                plot(mean(tmp_mean_SO),'color',[1 0 0])
                plot(mean(tmp_mean_control),'color',[0.5 0.5 0.5])

            subplot(3,2,5)
                hold on;
                tmp_SO_mod = [];
                tmp_control_mod = [];
                for n=1:length(all_m1_SO)
                    scatter(1,all_m1_SO(n).pre_m1_SO_maxFR-all_m1_SO(n).pre_m1_SO_minFR,[30],[1 0 0]);
                    scatter(2,all_m1_SO(n).pre_m1_control_maxFR-all_m1_SO(n).pre_m1_control_minFR,[30],[0.5 0.5 0.5]);
                    tmp_SO_mod = [tmp_SO_mod all_m1_SO(n).pre_m1_SO_maxFR-all_m1_SO(n).pre_m1_SO_minFR];
                    tmp_control_mod = [tmp_control_mod all_m1_SO(n).pre_m1_control_maxFR-all_m1_SO(n).pre_m1_control_minFR];
                end
                errorbar(1,mean(tmp_SO_mod),std(tmp_SO_mod)/sqrt(length(tmp_SO_mod)),'color',[1 0 0],'LineWidth',2,'CapSize',20);
                errorbar(2,mean(tmp_control_mod),std(tmp_control_mod)/sqrt(length(tmp_control_mod)),'color',[0 0 0],'LineWidth',2,'CapSize',20);
                xlim([0 3])

            subplot(3,2,6)
                hold on;
                tmp_SO_mod = [];
                tmp_control_mod = [];
                for n=1:length(all_m1_SO)
                    scatter(1,all_m1_SO(n).post_m1_SO_maxFR-all_m1_SO(n).post_m1_SO_minFR,[30],[1 0 0]);
                    scatter(2,all_m1_SO(n).post_m1_control_maxFR-all_m1_SO(n).post_m1_control_minFR,[30],[0.5 0.5 0.5]);
                    tmp_SO_mod = [tmp_SO_mod all_m1_SO(n).post_m1_SO_maxFR-all_m1_SO(n).post_m1_SO_minFR];
                    tmp_control_mod = [tmp_control_mod all_m1_SO(n).post_m1_control_maxFR-all_m1_SO(n).post_m1_control_minFR];
                end
                errorbar(1,mean(tmp_SO_mod),std(tmp_SO_mod)/sqrt(length(tmp_SO_mod)),'color',[1 0 0],'LineWidth',2,'CapSize',20);
                errorbar(2,mean(tmp_control_mod),std(tmp_control_mod)/sqrt(length(tmp_control_mod)),'color',[0 0 0],'LineWidth',2,'CapSize',20);
                xlim([0 3]) 
                
    %%% DELTA PLOT
    
        figure;

            subplot(3,2,1)
            
                delta_tot = [0 0];
                control_tot = [0 0];
                for n=1:length(all_m1_delta)
                    [x y] = pol2cart(all_m1_delta(n).pre_m1_plv_delta(1), all_m1_delta(n).pre_m1_plv_delta(2));
                    if isnan(x) || isnan(y); continue; end
                    delta_tot = delta_tot + [x y];
                    [x y] = pol2cart(all_m1_delta(n).pre_m1_plv_delta_control(1), all_m1_delta(n).pre_m1_plv_delta_control(2));
                    if isnan(x) || isnan(y); continue; end
                    control_tot = control_tot + [x y];
                end

                polarplot([0 all_m1_delta(1).pre_m1_plv_delta(1)],[0 all_m1_delta(1).pre_m1_plv_delta(2)],'color',[1 0.5 0.5],'LineWidth',1); 
                hold on;
                polarplot([0 all_m1_delta(1).pre_m1_plv_delta_control(1)],[0 all_m1_delta(1).pre_m1_plv_delta_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                for n=2:length(all_m1_delta)
                    polarplot([0 all_m1_delta(n).pre_m1_plv_delta(1)],[0 all_m1_delta(n).pre_m1_plv_delta(2)],'color',[1 0.5 0.5],'LineWidth',1);
                    polarplot([0 all_m1_delta(n).pre_m1_plv_delta_control(1)],[0 all_m1_delta(n).pre_m1_plv_delta_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                end
                
                [theta rho] = cart2pol(delta_tot(1), delta_tot(2));
                polarplot([0 theta],[0 rho/length(all_m1_delta)],'color','r','LineWidth',3);
                
                [theta rho] = cart2pol(control_tot(1), control_tot(2));
                polarplot([0 theta],[0 rho/length(all_m1_delta)],'color','k','LineWidth',3);

            subplot(3,2,2)
            
                delta_tot = [0 0];
                control_tot = [0 0];
                for n=1:length(all_m1_delta)
                    [x y] = pol2cart(all_m1_delta(n).post_m1_plv_delta(1), all_m1_delta(n).post_m1_plv_delta(2));
                    delta_tot = delta_tot + [x y];
                    [x y] = pol2cart(all_m1_delta(n).post_m1_plv_delta_control(1), all_m1_delta(n).post_m1_plv_delta_control(2));
                    control_tot = control_tot + [x y];
                end

                polarplot([0 all_m1_delta(1).post_m1_plv_delta(1)],[0 all_m1_delta(1).post_m1_plv_delta(2)],'color',[1 0.5 0.5],'LineWidth',1); 
                hold on;
                polarplot([0 all_m1_delta(1).post_m1_plv_delta_control(1)],[0 all_m1_delta(1).post_m1_plv_delta_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                for n=2:length(all_m1_delta)
                    polarplot([0 all_m1_delta(n).post_m1_plv_delta(1)],[0 all_m1_delta(n).post_m1_plv_delta(2)],'color',[1 0.5 0.5],'LineWidth',1);
                    polarplot([0 all_m1_delta(n).post_m1_plv_delta_control(1)],[0 all_m1_delta(n).post_m1_plv_delta_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                end
                
                [theta rho] = cart2pol(delta_tot(1), delta_tot(2));
                polarplot([0 theta],[0 rho/length(all_m1_delta)],'color','r','LineWidth',3);
                
                [theta rho] = cart2pol(control_tot(1), control_tot(2));
                polarplot([0 theta],[0 rho/length(all_m1_delta)],'color','k','LineWidth',3);


            subplot(3,2,3)
                hold on;
                tmp_mean_delta = [];
                tmp_mean_control = [];
                for n=1:length(all_m1_delta)
                    tmp_mean_delta = [tmp_mean_delta; mean(all_m1_delta(n).pre_m1_delta_raster)];
                    tmp_mean_control = [tmp_mean_control; mean(all_m1_delta(n).pre_m1_delta_control_raster)];
                end
                plot(mean(tmp_mean_delta),'color',[1 0 0])
                plot(mean(tmp_mean_control),'color',[0.5 0.5 0.5])

            subplot(3,2,4)
                hold on;
                tmp_mean_delta = [];
                tmp_mean_control = [];
                for n=1:length(all_m1_delta)
                    tmp_mean_delta = [tmp_mean_delta; mean(all_m1_delta(n).post_m1_delta_raster)];
                    tmp_mean_control = [tmp_mean_control; mean(all_m1_delta(n).post_m1_delta_control_raster)];
                end
                plot(mean(tmp_mean_delta),'color',[1 0 0])
                plot(mean(tmp_mean_control),'color',[0.5 0.5 0.5])

            subplot(3,2,5)
                hold on;
                tmp_delta_mod = [];
                tmp_control_mod = [];
                for n=1:length(all_m1_delta)
                    scatter(1,all_m1_delta(n).pre_m1_delta_maxFR-all_m1_delta(n).pre_m1_delta_minFR,[30],[1 0 0]);
                    scatter(2,all_m1_delta(n).pre_m1_control_maxFR-all_m1_delta(n).pre_m1_control_minFR,[30],[0.5 0.5 0.5]);
                    tmp_delta_mod = [tmp_delta_mod all_m1_delta(n).pre_m1_delta_maxFR-all_m1_delta(n).pre_m1_delta_minFR];
                    tmp_control_mod = [tmp_control_mod all_m1_delta(n).pre_m1_control_maxFR-all_m1_delta(n).pre_m1_control_minFR];
                end
                errorbar(1,mean(tmp_delta_mod),std(tmp_delta_mod)/sqrt(length(tmp_delta_mod)),'color',[1 0 0],'LineWidth',2,'CapSize',20);
                errorbar(2,mean(tmp_control_mod),std(tmp_control_mod)/sqrt(length(tmp_control_mod)),'color',[0 0 0],'LineWidth',2,'CapSize',20);
                xlim([0 3])

            subplot(3,2,6)
                hold on;
                tmp_delta_mod = [];
                tmp_control_mod = [];
                for n=1:length(all_m1_delta)
                    scatter(1,all_m1_delta(n).post_m1_delta_maxFR-all_m1_delta(n).post_m1_delta_minFR,[30],[1 0 0]);
                    scatter(2,all_m1_delta(n).post_m1_control_maxFR-all_m1_delta(n).post_m1_control_minFR,[30],[0.5 0.5 0.5]);
                    tmp_delta_mod = [tmp_delta_mod all_m1_delta(n).post_m1_delta_maxFR-all_m1_delta(n).post_m1_delta_minFR];
                    tmp_control_mod = [tmp_control_mod all_m1_delta(n).post_m1_control_maxFR-all_m1_delta(n).post_m1_control_minFR];
                end
                errorbar(1,mean(tmp_delta_mod),std(tmp_delta_mod)/sqrt(length(tmp_delta_mod)),'color',[1 0 0],'LineWidth',2,'CapSize',20);
                errorbar(2,mean(tmp_control_mod),std(tmp_control_mod)/sqrt(length(tmp_control_mod)),'color',[0 0 0],'LineWidth',2,'CapSize',20);
                xlim([0 3]) 
             
%% DLS RHYTHM MODULATION & PHASE LOCKING
    
    unit_count = 1;
    
    for dls_chan = 1:size(pre_dls_spiking,1)
        for dls_unit = 1:size(pre_dls_spiking,2)
            if isempty(pre_dls_spiking{dls_chan,dls_unit}) || isempty(post_dls_spiking{dls_chan,dls_unit})
                continue
            end

            all_dls_spindle(unit_count).dls_chan_unit = [dls_chan dls_unit];
            all_dls_SO(unit_count).dls_chan_unit = [dls_chan dls_unit];
            all_dls_delta(unit_count).dls_chan_unit = [dls_chan dls_unit];
            
            %%% PRE %%%

                pre_dls_unit = pre_dls_spiking{dls_chan,dls_unit}(pre_dls_spiking{dls_chan,dls_unit}>pre_video_start_stop(1) & pre_dls_spiking{dls_chan,dls_unit}<pre_video_start_stop(2))-pre_video_start_stop(1);

                %%% SPINDLES

                    spindle_raster = [];
                    spindle_ms_raster = [];                    
                    spindle_mod = [];
                    spindle_phase = [];
                    spindle_lfp = [];

                    control_raster = [];
                    control_mod = [];
                    control_phase = [];
                    control_lfp = [];

                    for spindle = 1:length(pre_sleep_rhythms.spindles{1,1}.pks)

                        control_pk =pre_sleep_rhythms.spindles{1,1}.pks(spindle)+randi([-10 -5],1,1);
                        if control_pk<=0.5
                            continue;
                        end
                        
                        spindle_lfp = [spindle_lfp; pre_dls_spindle_lfp(((pre_sleep_rhythms.spindles{1,1}.pks(spindle)-0.5)*pre_lfp.lfp_samp_rate):((pre_sleep_rhythms.spindles{1,1}.pks(spindle)+0.5)*pre_lfp.lfp_samp_rate))];
                        control_lfp = [control_lfp; pre_dls_spindle_lfp(((control_pk-0.5)*pre_lfp.lfp_samp_rate):((control_pk+0.5)*pre_lfp.lfp_samp_rate))];

                        spindle_raster = [spindle_raster; histcounts(pre_dls_unit((pre_dls_unit>=pre_sleep_rhythms.spindles{1,1}.pks(spindle)-1) & (pre_dls_unit<=pre_sleep_rhythms.spindles{1,1}.pks(spindle)+1))-pre_sleep_rhythms.spindles{1,1}.pks(spindle),[-1:0.01:1])];
                        spindle_ms_raster = [spindle_ms_raster; histcounts(pre_dls_unit((pre_dls_unit>=pre_sleep_rhythms.spindles{1,1}.pks(spindle)-1) & (pre_dls_unit<=pre_sleep_rhythms.spindles{1,1}.pks(spindle)+1))-pre_sleep_rhythms.spindles{1,1}.pks(spindle),[-1:0.001:1])];
                        control_raster = [control_raster; histcounts(pre_dls_unit((pre_dls_unit>=control_pk-1) & (pre_dls_unit<=control_pk+1))-control_pk,[-1:0.01:1])];                

                        spindle_mod = [spindle_mod; histcounts(pre_dls_unit((pre_dls_unit>=pre_sleep_rhythms.spindles{1,1}.pks(spindle)-1) & (pre_dls_unit<=pre_sleep_rhythms.spindles{1,1}.pks(spindle)+1))-pre_sleep_rhythms.spindles{1,1}.pks(spindle),[-1:0.04:1])];
                        control_mod = [control_mod; histcounts(pre_dls_unit((pre_dls_unit>=control_pk-1) & (pre_dls_unit<=control_pk+1))-control_pk,[-1:0.04:1])];                

                        tmp_phase = angle(hilbert(pre_m1_spindle_lfp(((pre_sleep_rhythms.spindles{1,1}.pks(spindle)-0.5)*pre_lfp.lfp_samp_rate):((pre_sleep_rhythms.spindles{1,1}.pks(spindle)+0.5)*pre_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(pre_dls_unit((pre_dls_unit>=pre_sleep_rhythms.spindles{1,1}.pks(spindle)-0.5) & (pre_dls_unit<=pre_sleep_rhythms.spindles{1,1}.pks(spindle)+0.5))-pre_sleep_rhythms.spindles{1,1}.pks(spindle),[-0.5:1/pre_lfp.lfp_samp_rate:.5]);
                        spindle_phase = [spindle_phase tmp_phase(tmp_spikes==1)];

                        tmp_phase = angle(hilbert(pre_m1_spindle_lfp(((control_pk-0.5)*pre_lfp.lfp_samp_rate):((control_pk+0.5)*pre_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(pre_dls_unit((pre_dls_unit>=control_pk-0.5) & (pre_dls_unit<=control_pk+0.5))-control_pk,[-0.5:1/pre_lfp.lfp_samp_rate:.5]);
                        control_phase = [control_phase tmp_phase(tmp_spikes==1)];

                    end

                    %%% PLV

                        spindle_plv = 0;
                        for n=1:length(spindle_phase)
                            spindle_plv = spindle_plv + exp(1i*spindle_phase(n));
                        end
                        spindle_plv_angle = angle(spindle_plv);
                        spindle_plv_mag = abs(spindle_plv)/length(spindle_phase);
                        
                        all_dls_spindle(unit_count).pre_dls_phase_spindle = spindle_phase; 
                        all_dls_spindle(unit_count).pre_dls_plv_spindle = [spindle_plv_angle spindle_plv_mag]; 

                        control_plv = 0;
                        for n=1:length(control_phase)
                            control_plv = control_plv + exp(1i*control_phase(n));
                        end
                        control_plv_angle = angle(control_plv);
                        control_plv_mag = abs(control_plv)/length(control_phase);
                        
                        all_dls_spindle(unit_count).pre_dls_phase_spindle_control = control_phase; 
                        all_dls_spindle(unit_count).pre_dls_plv_spindle_control = [control_plv_angle control_plv_mag]; 

                    %%% MODULATION

                        all_dls_spindle(unit_count).pre_dls_spindle_maxFR = max(mean(spindle_mod(:,13:38)))*25;
                        all_dls_spindle(unit_count).pre_dls_spindle_minFR = min(mean(spindle_mod(:,13:38)))*25;

                        all_dls_spindle(unit_count).pre_dls_control_maxFR = max(mean(control_mod(:,13:38)))*25;
                        all_dls_spindle(unit_count).pre_dls_control_minFR = min(mean(control_mod(:,13:38)))*25;

                        all_dls_spindle(unit_count).pre_dls_spindle_raster = spindle_raster;
                        all_dls_spindle(unit_count).pre_dls_spindle_ms_raster = spindle_ms_raster;
                        all_dls_spindle(unit_count).pre_dls_spindle_control_raster = control_raster;

                        all_dls_spindle(unit_count).pre_dls_spindle_lfp = spindle_lfp;
                        all_dls_spindle(unit_count).pre_dls_spindle_control_lfp = control_lfp;
                        
%                     figure; 
% 
%                         subplot(4,2,[1 2]);
%                             hold on;
%                             plot(mean(spindle_raster),'color','r')
%                             plot(mean(control_raster),'color','k')
%                             xlim([1 200])
%                             title('spindle vs. control raster')
% 
%                         subplot(4,2,[3 4]);
%                             hold on;
%                             plot(mean(spindle_lfp),'r')
%                             plot(mean(control_lfp),'k')
%                             xlim([1 round(pre_lfp.lfp_samp_rate)])
% 
%                         subplot(4,2,5);
%                             hold on;
%                             plot(mean(spindle_mod),'color','k')
%                             [~, max_i] = max(mean(spindle_mod(:,13:38)));
%                             [~, min_i] = min(mean(spindle_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['spindle modulation: ' num2str(all_dls_spindle(unit_count).pre_dls_spindle_maxFR-all_dls_spindle(unit_count).pre_dls_spindle_minFR)])
% 
%                         subplot(4,2,6);
%                            hold on;
%                             plot(mean(control_mod),'color','k')
%                             [~, max_i] = max(mean(control_mod(:,13:38)));
%                             [~, min_i] = min(mean(control_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['control modulation: ' num2str(all_dls_spindle(unit_count).pre_dls_control_maxFR-all_dls_spindle(unit_count).pre_dls_control_minFR)])
% 
%                         subplot(4,2,7);
%                             polarhistogram(spindle_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 spindle_plv_angle],[0 spindle_plv_mag],'LineWidth',3)
%                             tmp_rlim = rlim;
% 
%                         subplot(4,2,8);
%                             polarhistogram(control_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 control_plv_angle],[0 control_plv_mag],'LineWidth',3)
%                             rlim(tmp_rlim)

                %%% SO

                    SO_raster = [];
                    SO_mod = [];
                    SO_phase = [];
                    SO_lfp = [];

                    control_raster = [];
                    control_mod = [];
                    control_phase = [];
                    control_lfp = [];

                    for SO = 1:length(pre_sleep_rhythms.so_delta.SO_up_states)

                        control_pk =pre_sleep_rhythms.so_delta.SO_up_states(SO)+randi([-10 -5],1,1);
                        if control_pk<=0.5
                            continue;
                        end
                        
                        SO_lfp = [SO_lfp; pre_dls_SO_lfp(((pre_sleep_rhythms.so_delta.SO_up_states(SO)-0.5)*pre_lfp.lfp_samp_rate):((pre_sleep_rhythms.so_delta.SO_up_states(SO)+0.5)*pre_lfp.lfp_samp_rate))];
                        control_lfp = [control_lfp; pre_dls_SO_lfp(((control_pk-0.5)*pre_lfp.lfp_samp_rate):((control_pk+0.5)*pre_lfp.lfp_samp_rate))];

                        SO_raster = [SO_raster; histcounts(pre_dls_unit((pre_dls_unit>=pre_sleep_rhythms.so_delta.SO_up_states(SO)-1) & (pre_dls_unit<=pre_sleep_rhythms.so_delta.SO_up_states(SO)+1))-pre_sleep_rhythms.so_delta.SO_up_states(SO),[-1:0.01:1])];
                        control_raster = [control_raster; histcounts(pre_dls_unit((pre_dls_unit>=control_pk-1) & (pre_dls_unit<=control_pk+1))-control_pk,[-1:0.01:1])];                

                        SO_mod = [SO_mod; histcounts(pre_dls_unit((pre_dls_unit>=pre_sleep_rhythms.so_delta.SO_up_states(SO)-1) & (pre_dls_unit<=pre_sleep_rhythms.so_delta.SO_up_states(SO)+1))-pre_sleep_rhythms.so_delta.SO_up_states(SO),[-1:0.04:1])];
                        control_mod = [control_mod; histcounts(pre_dls_unit((pre_dls_unit>=control_pk-1) & (pre_dls_unit<=control_pk+1))-control_pk,[-1:0.04:1])];                

                        tmp_phase = angle(hilbert(pre_m1_SO_lfp(((pre_sleep_rhythms.so_delta.SO_up_states(SO)-0.5)*pre_lfp.lfp_samp_rate):((pre_sleep_rhythms.so_delta.SO_up_states(SO)+0.5)*pre_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(pre_dls_unit((pre_dls_unit>=pre_sleep_rhythms.so_delta.SO_up_states(SO)-0.5) & (pre_dls_unit<=pre_sleep_rhythms.so_delta.SO_up_states(SO)+0.5))-pre_sleep_rhythms.so_delta.SO_up_states(SO),[-0.5:1/pre_lfp.lfp_samp_rate:.5]);
                        SO_phase = [SO_phase tmp_phase(tmp_spikes==1)];

                        tmp_phase = angle(hilbert(pre_m1_SO_lfp(((control_pk-0.5)*pre_lfp.lfp_samp_rate):((control_pk+0.5)*pre_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(pre_dls_unit((pre_dls_unit>=control_pk-0.5) & (pre_dls_unit<=control_pk+0.5))-control_pk,[-0.5:1/pre_lfp.lfp_samp_rate:.5]);
                        control_phase = [control_phase tmp_phase(tmp_spikes==1)];

                    end

                    %%% PLV

                        SO_plv = 0;
                        for n=1:length(SO_phase)
                            SO_plv = SO_plv + exp(1i*SO_phase(n));
                        end
                        SO_plv_angle = angle(SO_plv);
                        SO_plv_mag = abs(SO_plv)/length(SO_phase);
                        
                        all_dls_SO(unit_count).pre_dls_phase_SO = SO_phase; 
                        all_dls_SO(unit_count).pre_dls_plv_SO = [SO_plv_angle SO_plv_mag]; 

                        control_plv = 0;
                        for n=1:length(control_phase)
                            control_plv = control_plv + exp(1i*control_phase(n));
                        end
                        control_plv_angle = angle(control_plv);
                        control_plv_mag = abs(control_plv)/length(control_phase);
                        
                        all_dls_SO(unit_count).pre_dls_phase_SO_control = control_phase; 
                        all_dls_SO(unit_count).pre_dls_plv_SO_control = [control_plv_angle control_plv_mag]; 

                    %%% MODULATION

                        all_dls_SO(unit_count).pre_dls_SO_maxFR = max(mean(SO_mod(:,13:38)))*25;
                        all_dls_SO(unit_count).pre_dls_SO_minFR = min(mean(SO_mod(:,13:38)))*25;

                        all_dls_SO(unit_count).pre_dls_control_maxFR = max(mean(control_mod(:,13:38)))*25;
                        all_dls_SO(unit_count).pre_dls_control_minFR = min(mean(control_mod(:,13:38)))*25;

                        all_dls_SO(unit_count).pre_dls_SO_raster = SO_raster;
                        all_dls_SO(unit_count).pre_dls_SO_control_raster = control_raster;

                        all_dls_SO(unit_count).pre_dls_SO_lfp = SO_lfp;
                        all_dls_SO(unit_count).pre_dls_SO_control_lfp = control_lfp;
                        
%                     figure; 
% 
%                         subplot(4,2,[1 2]);
%                             hold on;
%                             plot(mean(SO_raster),'color','r')
%                             plot(mean(control_raster),'color','k')
%                             xlim([1 200])
%                             title('SO vs. control raster')
% 
%                         subplot(4,2,[3 4]);
%                             hold on;
%                             plot(mean(SO_lfp),'r')
%                             plot(mean(control_lfp),'k')
%                             xlim([1 round(pre_lfp.lfp_samp_rate)])
% 
%                         subplot(4,2,5);
%                             hold on;
%                             plot(mean(SO_mod),'color','k')
%                             [~, max_i] = max(mean(SO_mod(:,13:38)));
%                             [~, min_i] = min(mean(SO_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['SO modulation: ' num2str(all_dls_SO(unit_count).pre_dls_SO_maxFR-all_dls_SO(unit_count).pre_dls_SO_minFR)])
% 
%                         subplot(4,2,6);
%                            hold on;
%                             plot(mean(control_mod),'color','k')
%                             [~, max_i] = max(mean(control_mod(:,13:38)));
%                             [~, min_i] = min(mean(control_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['control modulation: ' num2str(all_dls_SO(unit_count).pre_dls_control_maxFR-all_dls_SO(unit_count).pre_dls_control_minFR)])
% 
%                         subplot(4,2,7);
%                             polarhistogram(SO_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 SO_plv_angle],[0 SO_plv_mag],'LineWidth',3)
%                             tmp_rlim = rlim;
% 
%                         subplot(4,2,8);
%                             polarhistogram(control_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 control_plv_angle],[0 control_plv_mag],'LineWidth',3)
%                             rlim(tmp_rlim)

                %%% delta

                    delta_raster = [];
                    delta_mod = [];
                    delta_phase = [];
                    delta_lfp = [];

                    control_raster = [];
                    control_mod = [];
                    control_phase = [];
                    control_lfp = [];

                    for delta = 1:length(pre_sleep_rhythms.so_delta.delta_up_states)

                        control_pk =pre_sleep_rhythms.so_delta.delta_up_states(delta)+randi([-10 -5],1,1);
                        if control_pk<=0.5
                            continue;
                        end
                        
                        delta_lfp = [delta_lfp; pre_dls_SO_lfp(((pre_sleep_rhythms.so_delta.delta_up_states(delta)-0.5)*pre_lfp.lfp_samp_rate):((pre_sleep_rhythms.so_delta.delta_up_states(delta)+0.5)*pre_lfp.lfp_samp_rate))];
                        control_lfp = [control_lfp; pre_dls_SO_lfp(((control_pk-0.5)*pre_lfp.lfp_samp_rate):((control_pk+0.5)*pre_lfp.lfp_samp_rate))];

                        delta_raster = [delta_raster; histcounts(pre_dls_unit((pre_dls_unit>=pre_sleep_rhythms.so_delta.delta_up_states(delta)-1) & (pre_dls_unit<=pre_sleep_rhythms.so_delta.delta_up_states(delta)+1))-pre_sleep_rhythms.so_delta.delta_up_states(delta),[-1:0.01:1])];
                        control_raster = [control_raster; histcounts(pre_dls_unit((pre_dls_unit>=control_pk-1) & (pre_dls_unit<=control_pk+1))-control_pk,[-1:0.01:1])];                

                        delta_mod = [delta_mod; histcounts(pre_dls_unit((pre_dls_unit>=pre_sleep_rhythms.so_delta.delta_up_states(delta)-1) & (pre_dls_unit<=pre_sleep_rhythms.so_delta.delta_up_states(delta)+1))-pre_sleep_rhythms.so_delta.delta_up_states(delta),[-1:0.04:1])];
                        control_mod = [control_mod; histcounts(pre_dls_unit((pre_dls_unit>=control_pk-1) & (pre_dls_unit<=control_pk+1))-control_pk,[-1:0.04:1])];                

                        tmp_phase = angle(hilbert(pre_m1_SO_lfp(((pre_sleep_rhythms.so_delta.delta_up_states(delta)-0.5)*pre_lfp.lfp_samp_rate):((pre_sleep_rhythms.so_delta.delta_up_states(delta)+0.5)*pre_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(pre_dls_unit((pre_dls_unit>=pre_sleep_rhythms.so_delta.delta_up_states(delta)-0.5) & (pre_dls_unit<=pre_sleep_rhythms.so_delta.delta_up_states(delta)+0.5))-pre_sleep_rhythms.so_delta.delta_up_states(delta),[-0.5:1/pre_lfp.lfp_samp_rate:.5]);
                        delta_phase = [delta_phase tmp_phase(tmp_spikes==1)];

                        tmp_phase = angle(hilbert(pre_m1_SO_lfp(((control_pk-0.5)*pre_lfp.lfp_samp_rate):((control_pk+0.5)*pre_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(pre_dls_unit((pre_dls_unit>=control_pk-0.5) & (pre_dls_unit<=control_pk+0.5))-control_pk,[-0.5:1/pre_lfp.lfp_samp_rate:.5]);
                        control_phase = [control_phase tmp_phase(tmp_spikes==1)];

                    end

                    %%% PLV

                        delta_plv = 0;
                        for n=1:length(delta_phase)
                            delta_plv = delta_plv + exp(1i*delta_phase(n));
                        end
                        delta_plv_angle = angle(delta_plv);
                        delta_plv_mag = abs(delta_plv)/length(delta_phase);
                        
                        all_dls_delta(unit_count).pre_dls_phase_delta = delta_phase; 
                        all_dls_delta(unit_count).pre_dls_plv_delta = [delta_plv_angle delta_plv_mag]; 

                        control_plv = 0;
                        for n=1:length(control_phase)
                            control_plv = control_plv + exp(1i*control_phase(n));
                        end
                        control_plv_angle = angle(control_plv);
                        control_plv_mag = abs(control_plv)/length(control_phase);
                        
                        all_dls_delta(unit_count).pre_dls_phase_delta_control = control_phase; 
                        all_dls_delta(unit_count).pre_dls_plv_delta_control = [control_plv_angle control_plv_mag]; 

                    %%% MODULATION

                        all_dls_delta(unit_count).pre_dls_delta_maxFR = max(mean(delta_mod(:,13:38)))*25;
                        all_dls_delta(unit_count).pre_dls_delta_minFR = min(mean(delta_mod(:,13:38)))*25;

                        all_dls_delta(unit_count).pre_dls_control_maxFR = max(mean(control_mod(:,13:38)))*25;
                        all_dls_delta(unit_count).pre_dls_control_minFR = min(mean(control_mod(:,13:38)))*25;

                        all_dls_delta(unit_count).pre_dls_delta_raster = delta_raster;
                        all_dls_delta(unit_count).pre_dls_delta_control_raster = control_raster;

                        all_dls_delta(unit_count).pre_dls_delta_lfp = delta_lfp;
                        all_dls_delta(unit_count).pre_dls_delta_control_lfp = control_lfp;
                        
%                     figure; 
% 
%                         subplot(4,2,[1 2]);
%                             hold on;
%                             plot(mean(delta_raster),'color','r')
%                             plot(mean(control_raster),'color','k')
%                             xlim([1 200])
%                             title('delta vs. control raster')
% 
%                         subplot(4,2,[3 4]);
%                             hold on;
%                             plot(mean(delta_lfp),'r')
%                             plot(mean(control_lfp),'k')
%                             xlim([1 round(pre_lfp.lfp_samp_rate)])
% 
%                         subplot(4,2,5);
%                             hold on;
%                             plot(mean(delta_mod),'color','k')
%                             [~, max_i] = max(mean(delta_mod(:,13:38)));
%                             [~, min_i] = min(mean(delta_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['delta modulation: ' num2str(all_dls_delta(unit_count).pre_dls_delta_maxFR-all_dls_delta(unit_count).pre_dls_delta_minFR)])
% 
%                         subplot(4,2,6);
%                            hold on;
%                             plot(mean(control_mod),'color','k')
%                             [~, max_i] = max(mean(control_mod(:,13:38)));
%                             [~, min_i] = min(mean(control_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['control modulation: ' num2str(all_dls_delta(unit_count).pre_dls_control_maxFR-all_dls_delta(unit_count).pre_dls_control_minFR)])
% 
%                         subplot(4,2,7);
%                             polarhistogram(delta_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 delta_plv_angle],[0 delta_plv_mag],'LineWidth',3)
%                             tmp_rlim = rlim;
% 
%                         subplot(4,2,8);
%                             polarhistogram(control_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 control_plv_angle],[0 control_plv_mag],'LineWidth',3)
%                             rlim(tmp_rlim)

            %%% POST %%%

                post_dls_unit = post_dls_spiking{dls_chan,dls_unit}(post_dls_spiking{dls_chan,dls_unit}>post_video_start_stop(1) & post_dls_spiking{dls_chan,dls_unit}<post_video_start_stop(2))-post_video_start_stop(1);

                %%% SPINDLES

                    spindle_raster = [];
                    spindle_ms_raster = [];
                    spindle_mod = [];
                    spindle_phase = [];
                    spindle_lfp = [];

                    control_raster = [];
                    control_mod = [];
                    control_phase = [];
                    control_lfp = [];

                    for spindle = 1:length(post_sleep_rhythms.spindles{1,1}.pks)

                        %%% FIND MAX PEAK IN CONTROL PERIOD (5-10 SECONDS PRIOR TO SPINDLE PEAK)

                        control_pk =post_sleep_rhythms.spindles{1,1}.pks(spindle)+randi([-10 -5],1,1);
                        if control_pk<=0.5
                            continue;
                        end
                        
                        spindle_lfp = [spindle_lfp; post_dls_spindle_lfp(((post_sleep_rhythms.spindles{1,1}.pks(spindle)-0.5)*post_lfp.lfp_samp_rate):((post_sleep_rhythms.spindles{1,1}.pks(spindle)+0.5)*post_lfp.lfp_samp_rate))];
                        control_lfp = [control_lfp; post_dls_spindle_lfp(((control_pk-0.5)*post_lfp.lfp_samp_rate):((control_pk+0.5)*post_lfp.lfp_samp_rate))];

                        spindle_raster = [spindle_raster; histcounts(post_dls_unit((post_dls_unit>=post_sleep_rhythms.spindles{1,1}.pks(spindle)-1) & (post_dls_unit<=post_sleep_rhythms.spindles{1,1}.pks(spindle)+1))-post_sleep_rhythms.spindles{1,1}.pks(spindle),[-1:0.01:1])];
                        spindle_ms_raster = [spindle_ms_raster; histcounts(post_dls_unit((post_dls_unit>=post_sleep_rhythms.spindles{1,1}.pks(spindle)-1) & (post_dls_unit<=post_sleep_rhythms.spindles{1,1}.pks(spindle)+1))-post_sleep_rhythms.spindles{1,1}.pks(spindle),[-1:0.001:1])];
                        control_raster = [control_raster; histcounts(post_dls_unit((post_dls_unit>=control_pk-1) & (post_dls_unit<=control_pk+1))-control_pk,[-1:0.01:1])];                

                        spindle_mod = [spindle_mod; histcounts(post_dls_unit((post_dls_unit>=post_sleep_rhythms.spindles{1,1}.pks(spindle)-1) & (post_dls_unit<=post_sleep_rhythms.spindles{1,1}.pks(spindle)+1))-post_sleep_rhythms.spindles{1,1}.pks(spindle),[-1:0.04:1])];
                        control_mod = [control_mod; histcounts(post_dls_unit((post_dls_unit>=control_pk-1) & (post_dls_unit<=control_pk+1))-control_pk,[-1:0.04:1])];                

                        tmp_phase = angle(hilbert(post_m1_spindle_lfp(((post_sleep_rhythms.spindles{1,1}.pks(spindle)-0.5)*post_lfp.lfp_samp_rate):((post_sleep_rhythms.spindles{1,1}.pks(spindle)+0.5)*post_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(post_dls_unit((post_dls_unit>=post_sleep_rhythms.spindles{1,1}.pks(spindle)-0.5) & (post_dls_unit<=post_sleep_rhythms.spindles{1,1}.pks(spindle)+0.5))-post_sleep_rhythms.spindles{1,1}.pks(spindle),[-0.5:1/post_lfp.lfp_samp_rate:.5]);
                        spindle_phase = [spindle_phase tmp_phase(tmp_spikes==1)];

                        tmp_phase = angle(hilbert(post_m1_spindle_lfp(((control_pk-0.5)*post_lfp.lfp_samp_rate):((control_pk+0.5)*post_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(post_dls_unit((post_dls_unit>=control_pk-0.5) & (post_dls_unit<=control_pk+0.5))-control_pk,[-0.5:1/post_lfp.lfp_samp_rate:.5]);
                        control_phase = [control_phase tmp_phase(tmp_spikes==1)];

                    end

                    %%% PLV

                        spindle_plv = 0;
                        for n=1:length(spindle_phase)
                            spindle_plv = spindle_plv + exp(1i*spindle_phase(n));
                        end
                        spindle_plv_angle = angle(spindle_plv);
                        spindle_plv_mag = abs(spindle_plv)/length(spindle_phase);
                        
                        all_dls_spindle(unit_count).post_dls_phase_spindle = spindle_phase; 
                        all_dls_spindle(unit_count).post_dls_plv_spindle = [spindle_plv_angle spindle_plv_mag]; 

                        control_plv = 0;
                        for n=1:length(control_phase)
                            control_plv = control_plv + exp(1i*control_phase(n));
                        end
                        control_plv_angle = angle(control_plv);
                        control_plv_mag = abs(control_plv)/length(control_phase);
                        
                        all_dls_spindle(unit_count).post_dls_phase_spindle_control = control_phase; 
                        all_dls_spindle(unit_count).post_dls_plv_spindle_control = [control_plv_angle control_plv_mag]; 

                    %%% MODULATION

                        all_dls_spindle(unit_count).post_dls_spindle_maxFR = max(mean(spindle_mod(:,13:38)))*25;
                        all_dls_spindle(unit_count).post_dls_spindle_minFR = min(mean(spindle_mod(:,13:38)))*25;

                        all_dls_spindle(unit_count).post_dls_control_maxFR = max(mean(control_mod(:,13:38)))*25;
                        all_dls_spindle(unit_count).post_dls_control_minFR = min(mean(control_mod(:,13:38)))*25;

                        all_dls_spindle(unit_count).post_dls_spindle_raster = spindle_raster;
                        all_dls_spindle(unit_count).post_dls_spindle_ms_raster = spindle_ms_raster;
                        all_dls_spindle(unit_count).post_dls_spindle_control_raster = control_raster;

                        all_dls_spindle(unit_count).post_dls_spindle_lfp = spindle_lfp;
                        all_dls_spindle(unit_count).post_dls_spindle_control_lfp = control_lfp;
                        
%                     figure; 
% 
%                         subplot(4,2,[1 2]);
%                             hold on;
%                             plot(mean(spindle_raster),'color','r')
%                             plot(mean(control_raster),'color','k')
%                             xlim([1 200])
%                             title('spindle vs. control raster')
% 
%                         subplot(4,2,[3 4]);
%                             hold on;
%                             plot(mean(spindle_lfp),'r')
%                             plot(mean(control_lfp),'k')
%                             xlim([1 round(post_lfp.lfp_samp_rate)])
% 
%                         subplot(4,2,5);
%                             hold on;
%                             plot(mean(spindle_mod),'color','k')
%                             [~, max_i] = max(mean(spindle_mod(:,13:38)));
%                             [~, min_i] = min(mean(spindle_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['spindle modulation: ' num2str(all_dls_spindle(unit_count).post_dls_spindle_maxFR-all_dls_spindle(unit_count).post_dls_spindle_minFR)])
% 
%                         subplot(4,2,6);
%                            hold on;
%                             plot(mean(control_mod),'color','k')
%                             [~, max_i] = max(mean(control_mod(:,13:38)));
%                             [~, min_i] = min(mean(control_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['control modulation: ' num2str(all_dls_spindle(unit_count).post_dls_control_maxFR-all_dls_spindle(unit_count).post_dls_control_minFR)])
% 
%                         subplot(4,2,7);
%                             polarhistogram(spindle_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 spindle_plv_angle],[0 spindle_plv_mag],'LineWidth',3)
%                             tmp_rlim = rlim;
% 
%                         subplot(4,2,8);
%                             polarhistogram(control_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 control_plv_angle],[0 control_plv_mag],'LineWidth',3)
%                             rlim(tmp_rlim)

                %%% SO

                    SO_raster = [];
                    SO_mod = [];
                    SO_phase = [];
                    SO_lfp = [];

                    control_raster = [];
                    control_mod = [];
                    control_phase = [];
                    control_lfp = [];

                    for SO = 1:length(post_sleep_rhythms.so_delta.SO_up_states)

                        control_pk =post_sleep_rhythms.so_delta.SO_up_states(SO)+randi([-10 -5],1,1);
                        if control_pk<=0.5
                            continue;
                        end
                        
                        SO_lfp = [SO_lfp; post_dls_SO_lfp(((post_sleep_rhythms.so_delta.SO_up_states(SO)-0.5)*post_lfp.lfp_samp_rate):((post_sleep_rhythms.so_delta.SO_up_states(SO)+0.5)*post_lfp.lfp_samp_rate))];
                        control_lfp = [control_lfp; post_dls_SO_lfp(((control_pk-0.5)*post_lfp.lfp_samp_rate):((control_pk+0.5)*post_lfp.lfp_samp_rate))];

                        SO_raster = [SO_raster; histcounts(post_dls_unit((post_dls_unit>=post_sleep_rhythms.so_delta.SO_up_states(SO)-1) & (post_dls_unit<=post_sleep_rhythms.so_delta.SO_up_states(SO)+1))-post_sleep_rhythms.so_delta.SO_up_states(SO),[-1:0.01:1])];
                        control_raster = [control_raster; histcounts(post_dls_unit((post_dls_unit>=control_pk-1) & (post_dls_unit<=control_pk+1))-control_pk,[-1:0.01:1])];                

                        SO_mod = [SO_mod; histcounts(post_dls_unit((post_dls_unit>=post_sleep_rhythms.so_delta.SO_up_states(SO)-1) & (post_dls_unit<=post_sleep_rhythms.so_delta.SO_up_states(SO)+1))-post_sleep_rhythms.so_delta.SO_up_states(SO),[-1:0.04:1])];
                        control_mod = [control_mod; histcounts(post_dls_unit((post_dls_unit>=control_pk-1) & (post_dls_unit<=control_pk+1))-control_pk,[-1:0.04:1])];                

                        tmp_phase = angle(hilbert(post_m1_SO_lfp(((post_sleep_rhythms.so_delta.SO_up_states(SO)-0.5)*post_lfp.lfp_samp_rate):((post_sleep_rhythms.so_delta.SO_up_states(SO)+0.5)*post_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(post_dls_unit((post_dls_unit>=post_sleep_rhythms.so_delta.SO_up_states(SO)-0.5) & (post_dls_unit<=post_sleep_rhythms.so_delta.SO_up_states(SO)+0.5))-post_sleep_rhythms.so_delta.SO_up_states(SO),[-0.5:1/post_lfp.lfp_samp_rate:.5]);
                        SO_phase = [SO_phase tmp_phase(tmp_spikes==1)];

                        tmp_phase = angle(hilbert(post_m1_SO_lfp(((control_pk-0.5)*post_lfp.lfp_samp_rate):((control_pk+0.5)*post_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(post_dls_unit((post_dls_unit>=control_pk-0.5) & (post_dls_unit<=control_pk+0.5))-control_pk,[-0.5:1/post_lfp.lfp_samp_rate:.5]);
                        control_phase = [control_phase tmp_phase(tmp_spikes==1)];

                    end

                    %%% PLV

                        SO_plv = 0;
                        for n=1:length(SO_phase)
                            SO_plv = SO_plv + exp(1i*SO_phase(n));
                        end
                        SO_plv_angle = angle(SO_plv);
                        SO_plv_mag = abs(SO_plv)/length(SO_phase);
                        
                        all_dls_SO(unit_count).post_dls_phase_SO = SO_phase; 
                        all_dls_SO(unit_count).post_dls_plv_SO = [SO_plv_angle SO_plv_mag]; 

                        control_plv = 0;
                        for n=1:length(control_phase)
                            control_plv = control_plv + exp(1i*control_phase(n));
                        end
                        control_plv_angle = angle(control_plv);
                        control_plv_mag = abs(control_plv)/length(control_phase);
                        
                        all_dls_SO(unit_count).post_dls_phase_SO_control = control_phase; 
                        all_dls_SO(unit_count).post_dls_plv_SO_control = [control_plv_angle control_plv_mag]; 

                    %%% MODULATION

                        all_dls_SO(unit_count).post_dls_SO_maxFR = max(mean(SO_mod(:,13:38)))*25;
                        all_dls_SO(unit_count).post_dls_SO_minFR = min(mean(SO_mod(:,13:38)))*25;

                        all_dls_SO(unit_count).post_dls_control_maxFR = max(mean(control_mod(:,13:38)))*25;
                        all_dls_SO(unit_count).post_dls_control_minFR = min(mean(control_mod(:,13:38)))*25;

                        all_dls_SO(unit_count).post_dls_SO_raster = SO_raster;
                        all_dls_SO(unit_count).post_dls_SO_control_raster = control_raster;

                        all_dls_SO(unit_count).post_dls_SO_lfp = SO_lfp;
                        all_dls_SO(unit_count).post_dls_SO_control_lfp = control_lfp;
                        
%                     figure; 
% 
%                         subplot(4,2,[1 2]);
%                             hold on;
%                             plot(mean(SO_raster),'color','r')
%                             plot(mean(control_raster),'color','k')
%                             xlim([1 200])
%                             title('SO vs. control raster')
% 
%                         subplot(4,2,[3 4]);
%                             hold on;
%                             plot(mean(SO_lfp),'r')
%                             plot(mean(control_lfp),'k')
%                             xlim([1 round(post_lfp.lfp_samp_rate)])
% 
%                         subplot(4,2,5);
%                             hold on;
%                             plot(mean(SO_mod),'color','k')
%                             [~, max_i] = max(mean(SO_mod(:,13:38)));
%                             [~, min_i] = min(mean(SO_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['SO modulation: ' num2str(all_dls_SO(unit_count).post_dls_SO_maxFR-all_dls_SO(unit_count).post_dls_SO_minFR)])
% 
%                         subplot(4,2,6);
%                            hold on;
%                             plot(mean(control_mod),'color','k')
%                             [~, max_i] = max(mean(control_mod(:,13:38)));
%                             [~, min_i] = min(mean(control_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['control modulation: ' num2str(all_dls_SO(unit_count).post_dls_control_maxFR-all_dls_SO(unit_count).post_dls_control_minFR)])
% 
%                         subplot(4,2,7);
%                             polarhistogram(SO_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 SO_plv_angle],[0 SO_plv_mag],'LineWidth',3)
%                             tmp_rlim = rlim;
% 
%                         subplot(4,2,8);
%                             polarhistogram(control_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 control_plv_angle],[0 control_plv_mag],'LineWidth',3)
%                             rlim(tmp_rlim)

                %%% delta

                    delta_raster = [];
                    delta_mod = [];
                    delta_phase = [];
                    delta_lfp = [];

                    control_raster = [];
                    control_mod = [];
                    control_phase = [];
                    control_lfp = [];

                    for delta = 1:length(post_sleep_rhythms.so_delta.delta_up_states)

                        control_pk =post_sleep_rhythms.so_delta.delta_up_states(delta)+randi([-10 -5],1,1);
                        if control_pk<=0.5
                            continue;
                        end
                        
                        delta_lfp = [delta_lfp; post_dls_SO_lfp(((post_sleep_rhythms.so_delta.delta_up_states(delta)-0.5)*post_lfp.lfp_samp_rate):((post_sleep_rhythms.so_delta.delta_up_states(delta)+0.5)*post_lfp.lfp_samp_rate))];
                        control_lfp = [control_lfp; post_dls_SO_lfp(((control_pk-0.5)*post_lfp.lfp_samp_rate):((control_pk+0.5)*post_lfp.lfp_samp_rate))];

                        delta_raster = [delta_raster; histcounts(post_dls_unit((post_dls_unit>=post_sleep_rhythms.so_delta.delta_up_states(delta)-1) & (post_dls_unit<=post_sleep_rhythms.so_delta.delta_up_states(delta)+1))-post_sleep_rhythms.so_delta.delta_up_states(delta),[-1:0.01:1])];
                        control_raster = [control_raster; histcounts(post_dls_unit((post_dls_unit>=control_pk-1) & (post_dls_unit<=control_pk+1))-control_pk,[-1:0.01:1])];                

                        delta_mod = [delta_mod; histcounts(post_dls_unit((post_dls_unit>=post_sleep_rhythms.so_delta.delta_up_states(delta)-1) & (post_dls_unit<=post_sleep_rhythms.so_delta.delta_up_states(delta)+1))-post_sleep_rhythms.so_delta.delta_up_states(delta),[-1:0.04:1])];
                        control_mod = [control_mod; histcounts(post_dls_unit((post_dls_unit>=control_pk-1) & (post_dls_unit<=control_pk+1))-control_pk,[-1:0.04:1])];                

                        tmp_phase = angle(hilbert(post_m1_SO_lfp(((post_sleep_rhythms.so_delta.delta_up_states(delta)-0.5)*post_lfp.lfp_samp_rate):((post_sleep_rhythms.so_delta.delta_up_states(delta)+0.5)*post_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(post_dls_unit((post_dls_unit>=post_sleep_rhythms.so_delta.delta_up_states(delta)-0.5) & (post_dls_unit<=post_sleep_rhythms.so_delta.delta_up_states(delta)+0.5))-post_sleep_rhythms.so_delta.delta_up_states(delta),[-0.5:1/post_lfp.lfp_samp_rate:.5]);
                        delta_phase = [delta_phase tmp_phase(tmp_spikes==1)];

                        tmp_phase = angle(hilbert(post_m1_SO_lfp(((control_pk-0.5)*post_lfp.lfp_samp_rate):((control_pk+0.5)*post_lfp.lfp_samp_rate)))); tmp_phase = tmp_phase(1:end-1);
                        tmp_spikes = histcounts(post_dls_unit((post_dls_unit>=control_pk-0.5) & (post_dls_unit<=control_pk+0.5))-control_pk,[-0.5:1/post_lfp.lfp_samp_rate:.5]);
                        control_phase = [control_phase tmp_phase(tmp_spikes==1)];

                    end

                    %%% PLV

                        delta_plv = 0;
                        for n=1:length(delta_phase)
                            delta_plv = delta_plv + exp(1i*delta_phase(n));
                        end
                        delta_plv_angle = angle(delta_plv);
                        delta_plv_mag = abs(delta_plv)/length(delta_phase);
                        
                        all_dls_delta(unit_count).post_dls_phase_delta = delta_phase; 
                        all_dls_delta(unit_count).post_dls_plv_delta = [delta_plv_angle delta_plv_mag]; 

                        control_plv = 0;
                        for n=1:length(control_phase)
                            control_plv = control_plv + exp(1i*control_phase(n));
                        end
                        control_plv_angle = angle(control_plv);
                        control_plv_mag = abs(control_plv)/length(control_phase);
                        
                        all_dls_delta(unit_count).post_dls_phase_delta_control = control_phase; 
                        all_dls_delta(unit_count).post_dls_plv_delta_control = [control_plv_angle control_plv_mag]; 

                    %%% MODULATION

                        all_dls_delta(unit_count).post_dls_delta_maxFR = max(mean(delta_mod(:,13:38)))*25;
                        all_dls_delta(unit_count).post_dls_delta_minFR = min(mean(delta_mod(:,13:38)))*25;

                        all_dls_delta(unit_count).post_dls_control_maxFR = max(mean(control_mod(:,13:38)))*25;
                        all_dls_delta(unit_count).post_dls_control_minFR = min(mean(control_mod(:,13:38)))*25;

                        all_dls_delta(unit_count).post_dls_delta_raster = delta_raster;
                        all_dls_delta(unit_count).post_dls_delta_control_raster = control_raster;

                        all_dls_delta(unit_count).post_dls_delta_lfp = delta_lfp;
                        all_dls_delta(unit_count).post_dls_delta_control_lfp = control_lfp;
                        
%                     figure; 
% 
%                         subplot(4,2,[1 2]);
%                             hold on;
%                             plot(mean(delta_raster),'color','r')
%                             plot(mean(control_raster),'color','k')
%                             xlim([1 200])
%                             title('delta vs. control raster')
% 
%                         subplot(4,2,[3 4]);
%                             hold on;
%                             plot(mean(delta_lfp),'r')
%                             plot(mean(control_lfp),'k')
%                             xlim([1 round(post_lfp.lfp_samp_rate)])
% 
%                         subplot(4,2,5);
%                             hold on;
%                             plot(mean(delta_mod),'color','k')
%                             [~, max_i] = max(mean(delta_mod(:,13:38)));
%                             [~, min_i] = min(mean(delta_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['delta modulation: ' num2str(all_dls_delta(unit_count).post_dls_delta_maxFR-all_dls_delta(unit_count).post_dls_delta_minFR)])
% 
%                         subplot(4,2,6);
%                            hold on;
%                             plot(mean(control_mod),'color','k')
%                             [~, max_i] = max(mean(control_mod(:,13:38)));
%                             [~, min_i] = min(mean(control_mod(:,13:38)));
%                             plot([max_i+12 max_i+12],[ylim],'color','r')
%                             plot([min_i+12 min_i+12],[ylim],'color','b')
%                             xlim([1 50])
%                             title(['control modulation: ' num2str(all_dls_delta(unit_count).post_dls_control_maxFR-all_dls_delta(unit_count).post_dls_control_minFR)])
% 
%                         subplot(4,2,7);
%                             polarhistogram(delta_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 delta_plv_angle],[0 delta_plv_mag],'LineWidth',3)
%                             tmp_rlim = rlim;
% 
%                         subplot(4,2,8);
%                             polarhistogram(control_phase,[-pi:pi/8:pi],'Normalization','Probability')
%                             hold on;
%                             polarplot([0 control_plv_angle],[0 control_plv_mag],'LineWidth',3)
%                             rlim(tmp_rlim)

            unit_count = unit_count + 1;
            pause(0.1);
            
        end
    end
   
    %%% SPINDLE PLOT
    
        figure;

            subplot(3,2,1)
            
                spindle_tot = [0 0];
                control_tot = [0 0];
                for n=1:length(all_dls_spindle)
                    [x y] = pol2cart(all_dls_spindle(n).pre_dls_plv_spindle(1), all_dls_spindle(n).pre_dls_plv_spindle(2));
                    if isnan(x) || isnan(y); continue; end
                    spindle_tot = spindle_tot + [x y];
                    [x y] = pol2cart(all_dls_spindle(n).pre_dls_plv_spindle_control(1), all_dls_spindle(n).pre_dls_plv_spindle_control(2));
                    if isnan(x) || isnan(y); continue; end
                    control_tot = control_tot + [x y];
                end

                polarplot([0 all_dls_spindle(1).pre_dls_plv_spindle(1)],[0 all_dls_spindle(1).pre_dls_plv_spindle(2)],'color',[1 0.5 0.5],'LineWidth',1); 
                hold on;
                polarplot([0 all_dls_spindle(1).pre_dls_plv_spindle_control(1)],[0 all_dls_spindle(1).pre_dls_plv_spindle_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                for n=2:length(all_dls_spindle)
                    polarplot([0 all_dls_spindle(n).pre_dls_plv_spindle(1)],[0 all_dls_spindle(n).pre_dls_plv_spindle(2)],'color',[1 0.5 0.5],'LineWidth',1);
                    polarplot([0 all_dls_spindle(n).pre_dls_plv_spindle_control(1)],[0 all_dls_spindle(n).pre_dls_plv_spindle_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                end
                
                [theta rho] = cart2pol(spindle_tot(1), spindle_tot(2));
                polarplot([0 theta],[0 rho/length(all_dls_spindle)],'color','r','LineWidth',3);
                
                [theta rho] = cart2pol(control_tot(1), control_tot(2));
                polarplot([0 theta],[0 rho/length(all_dls_spindle)],'color','k','LineWidth',3);

            subplot(3,2,2)
            
                spindle_tot = [0 0];
                control_tot = [0 0];
                for n=1:length(all_dls_spindle)
                    [x y] = pol2cart(all_dls_spindle(n).post_dls_plv_spindle(1), all_dls_spindle(n).post_dls_plv_spindle(2));
                    if isnan(x) || isnan(y); continue; end
                    spindle_tot = spindle_tot + [x y];
                    [x y] = pol2cart(all_dls_spindle(n).post_dls_plv_spindle_control(1), all_dls_spindle(n).post_dls_plv_spindle_control(2));
                    if isnan(x) || isnan(y); continue; end
                    control_tot = control_tot + [x y];
                end

                polarplot([0 all_dls_spindle(1).post_dls_plv_spindle(1)],[0 all_dls_spindle(1).post_dls_plv_spindle(2)],'color',[1 0.5 0.5],'LineWidth',1); 
                hold on;
                polarplot([0 all_dls_spindle(1).post_dls_plv_spindle_control(1)],[0 all_dls_spindle(1).post_dls_plv_spindle_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                for n=2:length(all_dls_spindle)
                    polarplot([0 all_dls_spindle(n).post_dls_plv_spindle(1)],[0 all_dls_spindle(n).post_dls_plv_spindle(2)],'color',[1 0.5 0.5],'LineWidth',1);
                    polarplot([0 all_dls_spindle(n).post_dls_plv_spindle_control(1)],[0 all_dls_spindle(n).post_dls_plv_spindle_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                end
                
                [theta rho] = cart2pol(spindle_tot(1), spindle_tot(2));
                polarplot([0 theta],[0 rho/length(all_dls_spindle)],'color','r','LineWidth',3);
                
                [theta rho] = cart2pol(control_tot(1), control_tot(2));
                polarplot([0 theta],[0 rho/length(all_dls_spindle)],'color','k','LineWidth',3);


            subplot(3,2,3)
                hold on;
                tmp_mean_spindle = [];
                tmp_mean_control = [];
                for n=1:length(all_dls_spindle)
                    tmp_mean_spindle = [tmp_mean_spindle; mean(all_dls_spindle(n).pre_dls_spindle_raster)];
                    tmp_mean_control = [tmp_mean_control; mean(all_dls_spindle(n).pre_dls_spindle_control_raster)];
                end
                plot(mean(tmp_mean_spindle),'color',[1 0 0])
                plot(mean(tmp_mean_control),'color',[0.5 0.5 0.5])

            subplot(3,2,4)
                hold on;
                tmp_mean_spindle = [];
                tmp_mean_control = [];
                for n=1:length(all_dls_spindle)
                    tmp_mean_spindle = [tmp_mean_spindle; mean(all_dls_spindle(n).post_dls_spindle_raster)];
                    tmp_mean_control = [tmp_mean_control; mean(all_dls_spindle(n).post_dls_spindle_control_raster)];
                end
                plot(mean(tmp_mean_spindle),'color',[1 0 0])
                plot(mean(tmp_mean_control),'color',[0.5 0.5 0.5])

            subplot(3,2,5)
                hold on;
                tmp_spindle_mod = [];
                tmp_control_mod = [];
                for n=1:length(all_dls_spindle)
                    scatter(1,all_dls_spindle(n).pre_dls_spindle_maxFR-all_dls_spindle(n).pre_dls_spindle_minFR,[30],[1 0 0]);
                    scatter(2,all_dls_spindle(n).pre_dls_control_maxFR-all_dls_spindle(n).pre_dls_control_minFR,[30],[0.5 0.5 0.5]);
                    tmp_spindle_mod = [tmp_spindle_mod all_dls_spindle(n).pre_dls_spindle_maxFR-all_dls_spindle(n).pre_dls_spindle_minFR];
                    tmp_control_mod = [tmp_control_mod all_dls_spindle(n).pre_dls_control_maxFR-all_dls_spindle(n).pre_dls_control_minFR];
                end
                errorbar(1,mean(tmp_spindle_mod),std(tmp_spindle_mod)/sqrt(length(tmp_spindle_mod)),'color',[1 0 0],'LineWidth',2,'CapSize',20);
                errorbar(2,mean(tmp_control_mod),std(tmp_control_mod)/sqrt(length(tmp_control_mod)),'color',[0 0 0],'LineWidth',2,'CapSize',20);
                xlim([0 3])

            subplot(3,2,6)
                hold on;
                tmp_spindle_mod = [];
                tmp_control_mod = [];
                for n=1:length(all_dls_spindle)
                    scatter(1,all_dls_spindle(n).post_dls_spindle_maxFR-all_dls_spindle(n).post_dls_spindle_minFR,[30],[1 0 0]);
                    scatter(2,all_dls_spindle(n).post_dls_control_maxFR-all_dls_spindle(n).post_dls_control_minFR,[30],[0.5 0.5 0.5]);
                    tmp_spindle_mod = [tmp_spindle_mod all_dls_spindle(n).post_dls_spindle_maxFR-all_dls_spindle(n).post_dls_spindle_minFR];
                    tmp_control_mod = [tmp_control_mod all_dls_spindle(n).post_dls_control_maxFR-all_dls_spindle(n).post_dls_control_minFR];
                end
                errorbar(1,mean(tmp_spindle_mod),std(tmp_spindle_mod)/sqrt(length(tmp_spindle_mod)),'color',[1 0 0],'LineWidth',2,'CapSize',20);
                errorbar(2,mean(tmp_control_mod),std(tmp_control_mod)/sqrt(length(tmp_control_mod)),'color',[0 0 0],'LineWidth',2,'CapSize',20);
                xlim([0 3])  
            
    %%% SO PLOT
    
        figure;

            subplot(3,2,1)
            
                SO_tot = [0 0];
                control_tot = [0 0];
                for n=1:length(all_dls_SO)
                    [x y] = pol2cart(all_dls_SO(n).pre_dls_plv_SO(1), all_dls_SO(n).pre_dls_plv_SO(2));
                    if isnan(x) || isnan(y); continue; end
                    SO_tot = SO_tot + [x y];
                    [x y] = pol2cart(all_dls_SO(n).pre_dls_plv_SO_control(1), all_dls_SO(n).pre_dls_plv_SO_control(2));
                    if isnan(x) || isnan(y); continue; end
                    control_tot = control_tot + [x y];
                end

                polarplot([0 all_dls_SO(1).pre_dls_plv_SO(1)],[0 all_dls_SO(1).pre_dls_plv_SO(2)],'color',[1 0.5 0.5],'LineWidth',1); 
                hold on;
                polarplot([0 all_dls_SO(1).pre_dls_plv_SO_control(1)],[0 all_dls_SO(1).pre_dls_plv_SO_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                for n=2:length(all_dls_SO)
                    polarplot([0 all_dls_SO(n).pre_dls_plv_SO(1)],[0 all_dls_SO(n).pre_dls_plv_SO(2)],'color',[1 0.5 0.5],'LineWidth',1);
                    polarplot([0 all_dls_SO(n).pre_dls_plv_SO_control(1)],[0 all_dls_SO(n).pre_dls_plv_SO_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                end
                
                [theta rho] = cart2pol(SO_tot(1), SO_tot(2));
                polarplot([0 theta],[0 rho/length(all_dls_SO)],'color','r','LineWidth',3);
                
                [theta rho] = cart2pol(control_tot(1), control_tot(2));
                polarplot([0 theta],[0 rho/length(all_dls_SO)],'color','k','LineWidth',3);

            subplot(3,2,2)
            
                SO_tot = [0 0];
                control_tot = [0 0];
                for n=1:length(all_dls_SO)
                    [x y] = pol2cart(all_dls_SO(n).post_dls_plv_SO(1), all_dls_SO(n).post_dls_plv_SO(2));
                    if isnan(x) || isnan(y); continue; end
                    SO_tot = SO_tot + [x y];
                    [x y] = pol2cart(all_dls_SO(n).post_dls_plv_SO_control(1), all_dls_SO(n).post_dls_plv_SO_control(2));
                    if isnan(x) || isnan(y); continue; end
                    control_tot = control_tot + [x y];
                end

                polarplot([0 all_dls_SO(1).post_dls_plv_SO(1)],[0 all_dls_SO(1).post_dls_plv_SO(2)],'color',[1 0.5 0.5],'LineWidth',1); 
                hold on;
                polarplot([0 all_dls_SO(1).post_dls_plv_SO_control(1)],[0 all_dls_SO(1).post_dls_plv_SO_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                for n=2:length(all_dls_SO)
                    polarplot([0 all_dls_SO(n).post_dls_plv_SO(1)],[0 all_dls_SO(n).post_dls_plv_SO(2)],'color',[1 0.5 0.5],'LineWidth',1);
                    polarplot([0 all_dls_SO(n).post_dls_plv_SO_control(1)],[0 all_dls_SO(n).post_dls_plv_SO_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                end
                
                [theta rho] = cart2pol(SO_tot(1), SO_tot(2));
                polarplot([0 theta],[0 rho/length(all_dls_SO)],'color','r','LineWidth',3);
                
                [theta rho] = cart2pol(control_tot(1), control_tot(2));
                polarplot([0 theta],[0 rho/length(all_dls_SO)],'color','k','LineWidth',3);


            subplot(3,2,3)
                hold on;
                tmp_mean_SO = [];
                tmp_mean_control = [];
                for n=1:length(all_dls_SO)
                    tmp_mean_SO = [tmp_mean_SO; mean(all_dls_SO(n).pre_dls_SO_raster)];
                    tmp_mean_control = [tmp_mean_control; mean(all_dls_SO(n).pre_dls_SO_control_raster)];
                end
                plot(mean(tmp_mean_SO),'color',[1 0 0])
                plot(mean(tmp_mean_control),'color',[0.5 0.5 0.5])

            subplot(3,2,4)
                hold on;
                tmp_mean_SO = [];
                tmp_mean_control = [];
                for n=1:length(all_dls_SO)
                    tmp_mean_SO = [tmp_mean_SO; mean(all_dls_SO(n).post_dls_SO_raster)];
                    tmp_mean_control = [tmp_mean_control; mean(all_dls_SO(n).post_dls_SO_control_raster)];
                end
                plot(mean(tmp_mean_SO),'color',[1 0 0])
                plot(mean(tmp_mean_control),'color',[0.5 0.5 0.5])

            subplot(3,2,5)
                hold on;
                tmp_SO_mod = [];
                tmp_control_mod = [];
                for n=1:length(all_dls_SO)
                    scatter(1,all_dls_SO(n).pre_dls_SO_maxFR-all_dls_SO(n).pre_dls_SO_minFR,[30],[1 0 0]);
                    scatter(2,all_dls_SO(n).pre_dls_control_maxFR-all_dls_SO(n).pre_dls_control_minFR,[30],[0.5 0.5 0.5]);
                    tmp_SO_mod = [tmp_SO_mod all_dls_SO(n).pre_dls_SO_maxFR-all_dls_SO(n).pre_dls_SO_minFR];
                    tmp_control_mod = [tmp_control_mod all_dls_SO(n).pre_dls_control_maxFR-all_dls_SO(n).pre_dls_control_minFR];
                end
                errorbar(1,mean(tmp_SO_mod),std(tmp_SO_mod)/sqrt(length(tmp_SO_mod)),'color',[1 0 0],'LineWidth',2,'CapSize',20);
                errorbar(2,mean(tmp_control_mod),std(tmp_control_mod)/sqrt(length(tmp_control_mod)),'color',[0 0 0],'LineWidth',2,'CapSize',20);
                xlim([0 3])

            subplot(3,2,6)
                hold on;
                tmp_SO_mod = [];
                tmp_control_mod = [];
                for n=1:length(all_dls_SO)
                    scatter(1,all_dls_SO(n).post_dls_SO_maxFR-all_dls_SO(n).post_dls_SO_minFR,[30],[1 0 0]);
                    scatter(2,all_dls_SO(n).post_dls_control_maxFR-all_dls_SO(n).post_dls_control_minFR,[30],[0.5 0.5 0.5]);
                    tmp_SO_mod = [tmp_SO_mod all_dls_SO(n).post_dls_SO_maxFR-all_dls_SO(n).post_dls_SO_minFR];
                    tmp_control_mod = [tmp_control_mod all_dls_SO(n).post_dls_control_maxFR-all_dls_SO(n).post_dls_control_minFR];
                end
                errorbar(1,mean(tmp_SO_mod),std(tmp_SO_mod)/sqrt(length(tmp_SO_mod)),'color',[1 0 0],'LineWidth',2,'CapSize',20);
                errorbar(2,mean(tmp_control_mod),std(tmp_control_mod)/sqrt(length(tmp_control_mod)),'color',[0 0 0],'LineWidth',2,'CapSize',20);
                xlim([0 3]) 
                
    %%% DELTA PLOT
    
        figure;

            subplot(3,2,1)
            
                delta_tot = [0 0];
                control_tot = [0 0];
                for n=1:length(all_dls_delta)
                    [x y] = pol2cart(all_dls_delta(n).pre_dls_plv_delta(1), all_dls_delta(n).pre_dls_plv_delta(2));
                    if isnan(x) || isnan(y); continue; end
                    delta_tot = delta_tot + [x y];
                    [x y] = pol2cart(all_dls_delta(n).pre_dls_plv_delta_control(1), all_dls_delta(n).pre_dls_plv_delta_control(2));
                    if isnan(x) || isnan(y); continue; end
                    control_tot = control_tot + [x y];
                end

                polarplot([0 all_dls_delta(1).pre_dls_plv_delta(1)],[0 all_dls_delta(1).pre_dls_plv_delta(2)],'color',[1 0.5 0.5],'LineWidth',1); 
                hold on;
                polarplot([0 all_dls_delta(1).pre_dls_plv_delta_control(1)],[0 all_dls_delta(1).pre_dls_plv_delta_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                for n=2:length(all_dls_delta)
                    polarplot([0 all_dls_delta(n).pre_dls_plv_delta(1)],[0 all_dls_delta(n).pre_dls_plv_delta(2)],'color',[1 0.5 0.5],'LineWidth',1);
                    polarplot([0 all_dls_delta(n).pre_dls_plv_delta_control(1)],[0 all_dls_delta(n).pre_dls_plv_delta_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                end
                
                [theta rho] = cart2pol(delta_tot(1), delta_tot(2));
                polarplot([0 theta],[0 rho/length(all_dls_delta)],'color','r','LineWidth',3);
                
                [theta rho] = cart2pol(control_tot(1), control_tot(2));
                polarplot([0 theta],[0 rho/length(all_dls_delta)],'color','k','LineWidth',3);

            subplot(3,2,2)
            
                delta_tot = [0 0];
                control_tot = [0 0];
                for n=1:length(all_dls_delta)
                    [x y] = pol2cart(all_dls_delta(n).post_dls_plv_delta(1), all_dls_delta(n).post_dls_plv_delta(2));
                    if isnan(x) || isnan(y); continue; end
                    delta_tot = delta_tot + [x y];
                    [x y] = pol2cart(all_dls_delta(n).post_dls_plv_delta_control(1), all_dls_delta(n).post_dls_plv_delta_control(2));
                    if isnan(x) || isnan(y); continue; end
                    control_tot = control_tot + [x y];
                end

                polarplot([0 all_dls_delta(1).post_dls_plv_delta(1)],[0 all_dls_delta(1).post_dls_plv_delta(2)],'color',[1 0.5 0.5],'LineWidth',1); 
                hold on;
                polarplot([0 all_dls_delta(1).post_dls_plv_delta_control(1)],[0 all_dls_delta(1).post_dls_plv_delta_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                for n=2:length(all_dls_delta)
                    polarplot([0 all_dls_delta(n).post_dls_plv_delta(1)],[0 all_dls_delta(n).post_dls_plv_delta(2)],'color',[1 0.5 0.5],'LineWidth',1);
                    polarplot([0 all_dls_delta(n).post_dls_plv_delta_control(1)],[0 all_dls_delta(n).post_dls_plv_delta_control(2)],'color',[0.5 0.5 0.5],'LineWidth',1);
                end
                
                [theta rho] = cart2pol(delta_tot(1), delta_tot(2));
                polarplot([0 theta],[0 rho/length(all_dls_delta)],'color','r','LineWidth',3);
                
                [theta rho] = cart2pol(control_tot(1), control_tot(2));
                polarplot([0 theta],[0 rho/length(all_dls_delta)],'color','k','LineWidth',3);


            subplot(3,2,3)
                hold on;
                tmp_mean_delta = [];
                tmp_mean_control = [];
                for n=1:length(all_dls_delta)
                    tmp_mean_delta = [tmp_mean_delta; mean(all_dls_delta(n).pre_dls_delta_raster)];
                    tmp_mean_control = [tmp_mean_control; mean(all_dls_delta(n).pre_dls_delta_control_raster)];
                end
                plot(mean(tmp_mean_delta),'color',[1 0 0])
                plot(mean(tmp_mean_control),'color',[0.5 0.5 0.5])

            subplot(3,2,4)
                hold on;
                tmp_mean_delta = [];
                tmp_mean_control = [];
                for n=1:length(all_dls_delta)
                    tmp_mean_delta = [tmp_mean_delta; mean(all_dls_delta(n).post_dls_delta_raster)];
                    tmp_mean_control = [tmp_mean_control; mean(all_dls_delta(n).post_dls_delta_control_raster)];
                end
                plot(mean(tmp_mean_delta),'color',[1 0 0])
                plot(mean(tmp_mean_control),'color',[0.5 0.5 0.5])

            subplot(3,2,5)
                hold on;
                tmp_delta_mod = [];
                tmp_control_mod = [];
                for n=1:length(all_dls_delta)
                    scatter(1,all_dls_delta(n).pre_dls_delta_maxFR-all_dls_delta(n).pre_dls_delta_minFR,[30],[1 0 0]);
                    scatter(2,all_dls_delta(n).pre_dls_control_maxFR-all_dls_delta(n).pre_dls_control_minFR,[30],[0.5 0.5 0.5]);
                    tmp_delta_mod = [tmp_delta_mod all_dls_delta(n).pre_dls_delta_maxFR-all_dls_delta(n).pre_dls_delta_minFR];
                    tmp_control_mod = [tmp_control_mod all_dls_delta(n).pre_dls_control_maxFR-all_dls_delta(n).pre_dls_control_minFR];
                end
                errorbar(1,mean(tmp_delta_mod),std(tmp_delta_mod)/sqrt(length(tmp_delta_mod)),'color',[1 0 0],'LineWidth',2,'CapSize',20);
                errorbar(2,mean(tmp_control_mod),std(tmp_control_mod)/sqrt(length(tmp_control_mod)),'color',[0 0 0],'LineWidth',2,'CapSize',20);
                xlim([0 3])

            subplot(3,2,6)
                hold on;
                tmp_delta_mod = [];
                tmp_control_mod = [];
                for n=1:length(all_dls_delta)
                    scatter(1,all_dls_delta(n).post_dls_delta_maxFR-all_dls_delta(n).post_dls_delta_minFR,[30],[1 0 0]);
                    scatter(2,all_dls_delta(n).post_dls_control_maxFR-all_dls_delta(n).post_dls_control_minFR,[30],[0.5 0.5 0.5]);
                    tmp_delta_mod = [tmp_delta_mod all_dls_delta(n).post_dls_delta_maxFR-all_dls_delta(n).post_dls_delta_minFR];
                    tmp_control_mod = [tmp_control_mod all_dls_delta(n).post_dls_control_maxFR-all_dls_delta(n).post_dls_control_minFR];
                end
                errorbar(1,mean(tmp_delta_mod),std(tmp_delta_mod)/sqrt(length(tmp_delta_mod)),'color',[1 0 0],'LineWidth',2,'CapSize',20);
                errorbar(2,mean(tmp_control_mod),std(tmp_control_mod)/sqrt(length(tmp_control_mod)),'color',[0 0 0],'LineWidth',2,'CapSize',20);
                xlim([0 3])     
                
%% CROSS CORRELATION
    
    pair_count = 1;
    
    for m1_chan = 1:size(pre_m1_spiking,1)
        for m1_unit = 1:size(pre_m1_spiking,2)
            if isempty(pre_m1_spiking{m1_chan,m1_unit}) || isempty(post_m1_spiking{m1_chan,m1_unit})
                continue
            end
            
            pre_m1_unit = pre_m1_spiking{m1_chan,m1_unit}(pre_m1_spiking{m1_chan,m1_unit}>pre_video_start_stop(1) & pre_m1_spiking{m1_chan,m1_unit}<pre_video_start_stop(2))-pre_video_start_stop(1);
            post_m1_unit = post_m1_spiking{m1_chan,m1_unit}(post_m1_spiking{m1_chan,m1_unit}>post_video_start_stop(1) & post_m1_spiking{m1_chan,m1_unit}<post_video_start_stop(2))-post_video_start_stop(1);
            
            for dls_chan = 1:size(pre_dls_spiking,1)
                for dls_unit = 1:size(pre_dls_spiking,2)
                    if isempty(pre_dls_spiking{dls_chan,dls_unit}) || isempty(post_dls_spiking{dls_chan,dls_unit})
                        continue
                    end
                    
                    disp(['M1 CHAN ' num2str(m1_chan) ' DLS CHAN ' num2str(dls_chan)]);
                    pre_dls_unit = pre_dls_spiking{dls_chan,dls_unit}(pre_dls_spiking{dls_chan,dls_unit}>pre_video_start_stop(1) & pre_dls_spiking{dls_chan,dls_unit}<pre_video_start_stop(2))-pre_video_start_stop(1);
                    post_dls_unit = post_dls_spiking{dls_chan,dls_unit}(post_dls_spiking{dls_chan,dls_unit}>post_video_start_stop(1) & post_dls_spiking{dls_chan,dls_unit}<post_video_start_stop(2))-post_video_start_stop(1);
                    
                    %%% PAIR INFORMATION
                    
                    all_pairs(pair_count).m1_chan_unit = [m1_chan m1_unit];
                    all_pairs(pair_count).dls_chan_unit = [dls_chan dls_unit];
                    
                    %%% PRE CONTROL

                        pre_control_m1 = [];
                        pre_control_dls = [];
                        for spindle = 1:length(pre_sleep_rhythms.spindles{1,1}.pks)
                            control_pk =pre_sleep_rhythms.spindles{1,1}.pks(spindle)+randi([-10 -5],1,1);
                            if control_pk<=0.5
                                continue;
                            end
                        
                            pre_control_m1 = [pre_control_m1; histcounts(pre_m1_unit((pre_m1_unit>=control_pk-.125) & (pre_m1_unit<=control_pk+.125))-control_pk,[-.125:0.001:.125])];
                            pre_control_dls = [pre_control_dls; histcounts(pre_dls_unit((pre_dls_unit>=control_pk-.125) & (pre_dls_unit<=control_pk+.125))-control_pk,[-.125:0.001:.125])];

                        end

                        real_cc =  xcorr(reshape(pre_control_m1',size(pre_control_m1,1)*size(pre_control_m1,2),1),reshape(pre_control_dls',size(pre_control_dls,1)*size(pre_control_dls,2),1),100,'normalized');
                        shuffle_cc = [];
                        for shuffles = 1:25
                            ind = randperm(size(pre_control_m1,1));
                            shuffle_cc = [shuffle_cc; xcorr(reshape(pre_control_m1',size(pre_control_m1,1)*size(pre_control_m1,2),1),reshape(pre_control_dls(ind,:)',size(pre_control_dls,1)*size(pre_control_dls,2),1),100,'normalized')'];
                        end
                        all_pairs(pair_count).pre_control_cc = real_cc;
                        all_pairs(pair_count).pre_control_cc_shuffle = mean(shuffle_cc);

                    %%% POST CONTROL

                        post_control_m1 = [];
                        post_control_dls = [];
                        for spindle = 1:length(post_sleep_rhythms.spindles{1,1}.pks)
                            control_pk =post_sleep_rhythms.spindles{1,1}.pks(spindle)+randi([-10 -5],1,1);
                            if control_pk<=0.5
                                continue;
                            end
                            post_control_m1 = [post_control_m1; histcounts(post_m1_unit((post_m1_unit>=control_pk-.125) & (post_m1_unit<=control_pk+.125))-control_pk,[-.125:0.001:.125])];
                            post_control_dls = [post_control_dls; histcounts(post_dls_unit((post_dls_unit>=control_pk-.125) & (post_dls_unit<=control_pk+.125))-control_pk,[-.125:0.001:.125])];

                        end

                        real_cc =  xcorr(reshape(post_control_m1',size(post_control_m1,1)*size(post_control_m1,2),1),reshape(post_control_dls',size(post_control_dls,1)*size(post_control_dls,2),1),100,'normalized');
                        shuffle_cc = [];
                        for shuffles = 1:25
                            ind = randperm(size(post_control_m1,1));
                            shuffle_cc = [shuffle_cc; xcorr(reshape(post_control_m1',size(post_control_m1,1)*size(post_control_m1,2),1),reshape(post_control_dls(ind,:)',size(post_control_dls,1)*size(post_control_dls,2),1),100,'normalized')'];
                        end
                        all_pairs(pair_count).post_control_cc = real_cc;
                        all_pairs(pair_count).post_control_cc_shuffle = mean(shuffle_cc);
                    
                    %%% PRE SPINDLE
                    
                        pre_spindle_m1 = [];
                        pre_spindle_dls = [];
                        for spindle = 1:length(pre_sleep_rhythms.spindles{1,1}.pks)
                            pre_spindle_m1 = [pre_spindle_m1; histcounts(pre_m1_unit((pre_m1_unit>=pre_sleep_rhythms.spindles{1,1}.pks(spindle)-.125) & (pre_m1_unit<=pre_sleep_rhythms.spindles{1,1}.pks(spindle)+.125))-pre_sleep_rhythms.spindles{1,1}.pks(spindle),[-.125:0.001:.125])];
                            pre_spindle_dls = [pre_spindle_dls; histcounts(pre_dls_unit((pre_dls_unit>=pre_sleep_rhythms.spindles{1,1}.pks(spindle)-.125) & (pre_dls_unit<=pre_sleep_rhythms.spindles{1,1}.pks(spindle)+.125))-pre_sleep_rhythms.spindles{1,1}.pks(spindle),[-.125:0.001:.125])];
                        end
                        
                        real_cc =  xcorr(reshape(pre_spindle_m1',size(pre_spindle_m1,1)*size(pre_spindle_m1,2),1),reshape(pre_spindle_dls',size(pre_spindle_dls,1)*size(pre_spindle_dls,2),1),100,'normalized');
                        shuffle_cc = [];
                        for shuffles = 1:25
                            ind = randperm(size(pre_spindle_m1,1));
                            shuffle_cc = [shuffle_cc; xcorr(reshape(pre_spindle_m1',size(pre_spindle_m1,1)*size(pre_spindle_m1,2),1),reshape(pre_spindle_dls(ind,:)',size(pre_spindle_dls,1)*size(pre_spindle_dls,2),1),100,'normalized')'];
                        end
                        all_pairs(pair_count).pre_spindle_cc = real_cc;
                        all_pairs(pair_count).pre_spindle_cc_shuffle = mean(shuffle_cc);
                        
                    %%% POST SPINDLE
                    
                        post_spindle_m1 = [];
                        post_spindle_dls = [];
                        for spindle = 1:length(post_sleep_rhythms.spindles{1,1}.pks)
                            post_spindle_m1 = [post_spindle_m1; histcounts(post_m1_unit((post_m1_unit>=post_sleep_rhythms.spindles{1,1}.pks(spindle)-.125) & (post_m1_unit<=post_sleep_rhythms.spindles{1,1}.pks(spindle)+.125))-post_sleep_rhythms.spindles{1,1}.pks(spindle),[-.125:0.001:.125])];
                            post_spindle_dls = [post_spindle_dls; histcounts(post_dls_unit((post_dls_unit>=post_sleep_rhythms.spindles{1,1}.pks(spindle)-.125) & (post_dls_unit<=post_sleep_rhythms.spindles{1,1}.pks(spindle)+.125))-post_sleep_rhythms.spindles{1,1}.pks(spindle),[-.125:0.001:.125])];
                        end
                        
                        real_cc =  xcorr(reshape(post_spindle_m1',size(post_spindle_m1,1)*size(post_spindle_m1,2),1),reshape(post_spindle_dls',size(post_spindle_dls,1)*size(post_spindle_dls,2),1),100,'normalized');
                        shuffle_cc = [];
                        for shuffles = 1:25
                            ind = randperm(size(post_spindle_m1,1));
                            shuffle_cc = [shuffle_cc; xcorr(reshape(post_spindle_m1',size(post_spindle_m1,1)*size(post_spindle_m1,2),1),reshape(post_spindle_dls(ind,:)',size(post_spindle_dls,1)*size(post_spindle_dls,2),1),100,'normalized')'];
                        end
                        all_pairs(pair_count).post_spindle_cc = real_cc;
                        all_pairs(pair_count).post_spindle_cc_shuffle = mean(shuffle_cc);

                    %%% PRE SO
                    
                        pre_so_m1 = [];
                        pre_so_dls = [];
                        for SO = 1:length(pre_sleep_rhythms.so_delta.SO_up_states)
                            pre_so_m1 = [pre_so_m1; histcounts(pre_m1_unit((pre_m1_unit>=pre_sleep_rhythms.so_delta.SO_up_states(SO)-.125) & (pre_m1_unit<=pre_sleep_rhythms.so_delta.SO_up_states(SO)+.125))-pre_sleep_rhythms.so_delta.SO_up_states(SO),[-.125:0.001:.125])];
                            pre_so_dls = [pre_so_dls; histcounts(pre_dls_unit((pre_dls_unit>=pre_sleep_rhythms.so_delta.SO_up_states(SO)-.125) & (pre_dls_unit<=pre_sleep_rhythms.so_delta.SO_up_states(SO)+.125))-pre_sleep_rhythms.so_delta.SO_up_states(SO),[-.125:0.001:.125])];
                        end
                        
                        real_cc =  xcorr(reshape(pre_so_m1',size(pre_so_m1,1)*size(pre_so_m1,2),1),reshape(pre_so_dls',size(pre_so_dls,1)*size(pre_so_dls,2),1),100,'normalized');
                        shuffle_cc = [];
                        for shuffles = 1:25
                            ind = randperm(size(pre_so_m1,1));
                            shuffle_cc = [shuffle_cc; xcorr(reshape(pre_so_m1',size(pre_so_m1,1)*size(pre_so_m1,2),1),reshape(pre_so_dls(ind,:)',size(pre_so_dls,1)*size(pre_so_dls,2),1),100,'normalized')'];
                        end
                        all_pairs(pair_count).pre_so_cc = real_cc;
                        all_pairs(pair_count).pre_so_cc_shuffle = mean(shuffle_cc);
                        
                    %%% POST SO
                    
                        post_so_m1 = [];
                        post_so_dls = [];
                        for SO = 1:length(post_sleep_rhythms.so_delta.SO_up_states)
                            post_so_m1 = [post_so_m1; histcounts(post_m1_unit((post_m1_unit>=post_sleep_rhythms.so_delta.SO_up_states(SO)-.125) & (post_m1_unit<=post_sleep_rhythms.so_delta.SO_up_states(SO)+.125))-post_sleep_rhythms.so_delta.SO_up_states(SO),[-.125:0.001:.125])];
                            post_so_dls = [post_so_dls; histcounts(post_dls_unit((post_dls_unit>=post_sleep_rhythms.so_delta.SO_up_states(SO)-.125) & (post_dls_unit<=post_sleep_rhythms.so_delta.SO_up_states(SO)+.125))-post_sleep_rhythms.so_delta.SO_up_states(SO),[-.125:0.001:.125])];
                        end
                        
                        real_cc =  xcorr(reshape(post_so_m1',size(post_so_m1,1)*size(post_so_m1,2),1),reshape(post_so_dls',size(post_so_dls,1)*size(post_so_dls,2),1),100,'normalized');
                        shuffle_cc = [];
                        for shuffles = 1:25
                            ind = randperm(size(post_so_m1,1));
                            shuffle_cc = [shuffle_cc; xcorr(reshape(post_so_m1',size(post_so_m1,1)*size(post_so_m1,2),1),reshape(post_so_dls(ind,:)',size(post_so_dls,1)*size(post_so_dls,2),1),100,'normalized')'];
                        end
                        all_pairs(pair_count).post_so_cc = real_cc;
                        all_pairs(pair_count).post_so_cc_shuffle = mean(shuffle_cc);
                        
                    %%% PRE DELTA
                    
                        pre_delta_m1 = [];
                        pre_delta_dls = [];
                        for delta_ind = 1:length(pre_sleep_rhythms.so_delta.delta_up_states)
                            pre_delta_m1 = [pre_delta_m1; histcounts(pre_m1_unit((pre_m1_unit>=pre_sleep_rhythms.so_delta.delta_up_states(delta_ind)-.125) & (pre_m1_unit<=pre_sleep_rhythms.so_delta.delta_up_states(delta_ind)+.125))-pre_sleep_rhythms.so_delta.delta_up_states(delta_ind),[-.125:0.001:.125])];
                            pre_delta_dls = [pre_delta_dls; histcounts(pre_dls_unit((pre_dls_unit>=pre_sleep_rhythms.so_delta.delta_up_states(delta_ind)-.125) & (pre_dls_unit<=pre_sleep_rhythms.so_delta.delta_up_states(delta_ind)+.125))-pre_sleep_rhythms.so_delta.delta_up_states(delta_ind),[-.125:0.001:.125])];
                        end
                        
                        real_cc =  xcorr(reshape(pre_delta_m1',size(pre_delta_m1,1)*size(pre_delta_m1,2),1),reshape(pre_delta_dls',size(pre_delta_dls,1)*size(pre_delta_dls,2),1),100,'normalized');
                        shuffle_cc = [];
                        for shuffles = 1:25
                            ind = randperm(size(pre_delta_m1,1));
                            shuffle_cc = [shuffle_cc; xcorr(reshape(pre_delta_m1',size(pre_delta_m1,1)*size(pre_delta_m1,2),1),reshape(pre_delta_dls(ind,:)',size(pre_delta_dls,1)*size(pre_delta_dls,2),1),100,'normalized')'];
                        end
                        all_pairs(pair_count).pre_delta_cc = real_cc;
                        all_pairs(pair_count).pre_delta_cc_shuffle = mean(shuffle_cc);
                        
                    %%% POST DELTA
                    
                        post_delta_m1 = [];
                        post_delta_dls = [];
                        for delta_ind = 1:length(post_sleep_rhythms.so_delta.delta_up_states)
                            post_delta_m1 = [post_delta_m1; histcounts(post_m1_unit((post_m1_unit>=post_sleep_rhythms.so_delta.delta_up_states(delta_ind)-.125) & (post_m1_unit<=post_sleep_rhythms.so_delta.delta_up_states(delta_ind)+.125))-post_sleep_rhythms.so_delta.delta_up_states(delta_ind),[-.125:0.001:.125])];
                            post_delta_dls = [post_delta_dls; histcounts(post_dls_unit((post_dls_unit>=post_sleep_rhythms.so_delta.delta_up_states(delta_ind)-.125) & (post_dls_unit<=post_sleep_rhythms.so_delta.delta_up_states(delta_ind)+.125))-post_sleep_rhythms.so_delta.delta_up_states(delta_ind),[-.125:0.001:.125])];
                        end
                        
                        real_cc =  xcorr(reshape(post_delta_m1',size(post_delta_m1,1)*size(post_delta_m1,2),1),reshape(post_delta_dls',size(post_delta_dls,1)*size(post_delta_dls,2),1),100,'normalized');
                        shuffle_cc = [];
                        for shuffles = 1:25
                            ind = randperm(size(post_delta_m1,1));
                            shuffle_cc = [shuffle_cc; xcorr(reshape(post_delta_m1',size(post_delta_m1,1)*size(post_delta_m1,2),1),reshape(post_delta_dls(ind,:)',size(post_delta_dls,1)*size(post_delta_dls,2),1),100,'normalized')'];
                        end
                        all_pairs(pair_count).post_delta_cc = real_cc;
                        all_pairs(pair_count).post_delta_cc_shuffle = mean(shuffle_cc);
                        
                        pair_count = pair_count + 1;
                end
            end
            
        end
    end
    
%     tmp_real = [];
%     tmp_shuffle = [];
%     for n=1:length(all_pairs)
%         tmp_real = [tmp_real; all_pairs(n).pre_spindle_cc'];
%         tmp_shuffle = [tmp_shuffle; all_pairs(n).pre_spindle_cc_shuffle];
%     end
%     figure;
%         subplot(2,1,1)
%             hold on;
%             plot(mean(tmp_real));
%             plot(mean(tmp_shuffle));
%         subplot(2,1,2)
%             plot(mean(tmp_real)-mean(tmp_shuffle));
%     
%     tmp_real = [];
%     tmp_shuffle = [];
%     for n=1:length(all_pairs)
%         tmp_real = [tmp_real; all_pairs(n).post_spindle_cc'];
%         tmp_shuffle = [tmp_shuffle; all_pairs(n).post_spindle_cc_shuffle];
%     end
%     figure;
%         subplot(2,1,1)
%             hold on;
%             plot(mean(tmp_real));
%             plot(mean(tmp_shuffle));
%         subplot(2,1,2)
%             plot(mean(tmp_real)-mean(tmp_shuffle));
%             
%     tmp_real = [];
%     tmp_shuffle = [];
%     for n=1:length(all_pairs)
%         tmp_real = [tmp_real; all_pairs(n).pre_so_cc'];
%         tmp_shuffle = [tmp_shuffle; all_pairs(n).pre_so_cc_shuffle];
%     end
%     figure;
%         subplot(2,1,1)
%             hold on;
%             plot(mean(tmp_real));
%             plot(mean(tmp_shuffle));
%         subplot(2,1,2)
%             plot(mean(tmp_real)-mean(tmp_shuffle));
%     
%     tmp_real = [];
%     tmp_shuffle = [];
%     for n=1:length(all_pairs)
%         tmp_real = [tmp_real; all_pairs(n).post_so_cc'];
%         tmp_shuffle = [tmp_shuffle; all_pairs(n).post_so_cc_shuffle];
%     end
%     figure;
%         subplot(2,1,1)
%             hold on;
%             plot(mean(tmp_real));
%             plot(mean(tmp_shuffle));
%         subplot(2,1,2)
%             plot(mean(tmp_real)-mean(tmp_shuffle));
%             
%     tmp_real = [];
%     tmp_shuffle = [];
%     for n=1:length(all_pairs)
%         tmp_real = [tmp_real; all_pairs(n).pre_delta_cc'];
%         tmp_shuffle = [tmp_shuffle; all_pairs(n).pre_delta_cc_shuffle];
%     end
%     figure;
%         subplot(2,1,1)
%             hold on;
%             plot(mean(tmp_real));
%             plot(mean(tmp_shuffle));
%         subplot(2,1,2)
%             plot(mean(tmp_real)-mean(tmp_shuffle));
%     
%     tmp_real = [];
%     tmp_shuffle = [];
%     for n=1:length(all_pairs)
%         tmp_real = [tmp_real; all_pairs(n).post_delta_cc'];
%         tmp_shuffle = [tmp_shuffle; all_pairs(n).post_delta_cc_shuffle];
%     end
%     figure;
%         subplot(2,1,1)
%             hold on;
%             plot(mean(tmp_real));
%             plot(mean(tmp_shuffle));
%         subplot(2,1,2)
%             plot(mean(tmp_real)-mean(tmp_shuffle));
            
end
