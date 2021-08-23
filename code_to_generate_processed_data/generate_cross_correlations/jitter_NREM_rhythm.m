function [all_pairs] = jitter_NREM_rhythm(pre_m1_spiking, pre_dls_spiking, post_m1_spiking, post_dls_spiking, pre_video_start_stop, post_video_start_stop,pre_sleep_rhythms,post_sleep_rhythms,sig)

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
                
                %%% CROSS CORRELATIONS
                
                if sig(pair_count)>.99
                    
                    %%% PRE SPINDLE
                    
                    pre_spindle_m1 = [];
                    pre_spindle_dls = [];
                    for spindle = 1:length(pre_sleep_rhythms.spindles{1,1}.pks)
                        pre_spindle_m1 = [pre_spindle_m1; histcounts(pre_m1_unit((pre_m1_unit>=pre_sleep_rhythms.spindles{1,1}.pks(spindle)-.25) & (pre_m1_unit<=pre_sleep_rhythms.spindles{1,1}.pks(spindle)+.25))-pre_sleep_rhythms.spindles{1,1}.pks(spindle),[-.25:0.001:.25])];
                        pre_spindle_dls = [pre_spindle_dls; histcounts(pre_dls_unit((pre_dls_unit>=pre_sleep_rhythms.spindles{1,1}.pks(spindle)-.25) & (pre_dls_unit<=pre_sleep_rhythms.spindles{1,1}.pks(spindle)+.25))-pre_sleep_rhythms.spindles{1,1}.pks(spindle),[-.25:0.001:.25])];
                    end
                    
                    real_cc =  xcorr(reshape(pre_spindle_m1',size(pre_spindle_m1,1)*size(pre_spindle_m1,2),1),reshape(pre_spindle_dls',size(pre_spindle_dls,1)*size(pre_spindle_dls,2),1),100,'normalized');
                    all_pairs(pair_count).pre_spindle_cc = real_cc;
                    
                    shuffle_cc = zeros(100,201);
                    parfor shuffle = 1:100
                        pre_spindle_dls_shufffle = zeros(size(pre_spindle_dls,1),size(pre_spindle_dls,2));
                        for n = 1:size(pre_spindle_m1,1)
                            tmp_dls_spikes = find(pre_spindle_dls(n,:))+randi([-25 25],1,length(find(pre_spindle_dls(n,:))));
                            tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<size(pre_spindle_dls,2) & tmp_dls_spikes>0);
                            pre_spindle_dls_shufffle(n,tmp_dls_spikes) = 1;
                        end
                        shuffle_cc(shuffle,:) = xcorr(reshape(pre_spindle_m1',size(pre_spindle_m1,1)*size(pre_spindle_m1,2),1),reshape(pre_spindle_dls_shufffle',size(pre_spindle_dls_shufffle,1)*size(pre_spindle_dls_shufffle,2),1),100,'normalized')';
                    end
                    all_pairs(pair_count).pre_spindle_cc_shuffle = shuffle_cc;
                    
                    %%% POST SPINDLE
                    
                    post_spindle_m1 = [];
                    post_spindle_dls = [];
                    for spindle = 1:length(post_sleep_rhythms.spindles{1,1}.pks)
                        post_spindle_m1 = [post_spindle_m1; histcounts(post_m1_unit((post_m1_unit>=post_sleep_rhythms.spindles{1,1}.pks(spindle)-.25) & (post_m1_unit<=post_sleep_rhythms.spindles{1,1}.pks(spindle)+.25))-post_sleep_rhythms.spindles{1,1}.pks(spindle),[-.25:0.001:.25])];
                        post_spindle_dls = [post_spindle_dls; histcounts(post_dls_unit((post_dls_unit>=post_sleep_rhythms.spindles{1,1}.pks(spindle)-.25) & (post_dls_unit<=post_sleep_rhythms.spindles{1,1}.pks(spindle)+.25))-post_sleep_rhythms.spindles{1,1}.pks(spindle),[-.25:0.001:.25])];
                    end
                    
                    real_cc =  xcorr(reshape(post_spindle_m1',size(post_spindle_m1,1)*size(post_spindle_m1,2),1),reshape(post_spindle_dls',size(post_spindle_dls,1)*size(post_spindle_dls,2),1),100,'normalized');
                    all_pairs(pair_count).post_spindle_cc = real_cc;
                    
                    shuffle_cc = zeros(100,201);
                    parfor shuffle = 1:100
                        post_spindle_dls_shufffle = zeros(size(post_spindle_dls,1),size(post_spindle_dls,2));
                        for n = 1:size(post_spindle_m1,1)
                            tmp_dls_spikes = find(post_spindle_dls(n,:))+randi([-25 25],1,length(find(post_spindle_dls(n,:))));
                            tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<size(post_spindle_dls,2) & tmp_dls_spikes>0);
                            post_spindle_dls_shufffle(n,tmp_dls_spikes) = 1;
                        end
                        shuffle_cc(shuffle,:) = xcorr(reshape(post_spindle_m1',size(post_spindle_m1,1)*size(post_spindle_m1,2),1),reshape(post_spindle_dls_shufffle',size(post_spindle_dls_shufffle,1)*size(post_spindle_dls_shufffle,2),1),100,'normalized')';
                    end
                    all_pairs(pair_count).post_spindle_cc_shuffle = shuffle_cc;
                    
                    %%% PRE SO
                    
                    pre_SO_m1 = [];
                    pre_SO_dls = [];
                    for SO = 1:length(pre_sleep_rhythms.so_delta.SO_up_states)
                        pre_SO_m1 = [pre_SO_m1; histcounts(pre_m1_unit((pre_m1_unit>=pre_sleep_rhythms.so_delta.SO_up_states(SO)-.25) & (pre_m1_unit<=pre_sleep_rhythms.so_delta.SO_up_states(SO)+.25))-pre_sleep_rhythms.so_delta.SO_up_states(SO),[-.25:0.001:.25])];
                        pre_SO_dls = [pre_SO_dls; histcounts(pre_dls_unit((pre_dls_unit>=pre_sleep_rhythms.so_delta.SO_up_states(SO)-.25) & (pre_dls_unit<=pre_sleep_rhythms.so_delta.SO_up_states(SO)+.25))-pre_sleep_rhythms.so_delta.SO_up_states(SO),[-.25:0.001:.25])];
                    end
                    
                    real_cc =  xcorr(reshape(pre_SO_m1',size(pre_SO_m1,1)*size(pre_SO_m1,2),1),reshape(pre_SO_dls',size(pre_SO_dls,1)*size(pre_SO_dls,2),1),100,'normalized');
                    all_pairs(pair_count).pre_SO_cc = real_cc;
                    
                    shuffle_cc = zeros(100,201);
                    parfor shuffle = 1:100
                        pre_SO_dls_shufffle = zeros(size(pre_SO_dls,1),size(pre_SO_dls,2));
                        for n = 1:size(pre_SO_m1,1)
                            tmp_dls_spikes = find(pre_SO_dls(n,:))+randi([-25 25],1,length(find(pre_SO_dls(n,:))));
                            tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<size(pre_SO_dls,2) & tmp_dls_spikes>0);
                            pre_SO_dls_shufffle(n,tmp_dls_spikes) = 1;
                        end
                        shuffle_cc(shuffle,:) = xcorr(reshape(pre_SO_m1',size(pre_SO_m1,1)*size(pre_SO_m1,2),1),reshape(pre_SO_dls_shufffle',size(pre_SO_dls_shufffle,1)*size(pre_SO_dls_shufffle,2),1),100,'normalized')';
                    end
                    all_pairs(pair_count).pre_SO_cc_shuffle = shuffle_cc;
                    
                    %%% POST SO
                    
                    post_SO_m1 = [];
                    post_SO_dls = [];
                    for SO = 1:length(post_sleep_rhythms.so_delta.SO_up_states)
                        post_SO_m1 = [post_SO_m1; histcounts(post_m1_unit((post_m1_unit>=post_sleep_rhythms.so_delta.SO_up_states(SO)-.25) & (post_m1_unit<=post_sleep_rhythms.so_delta.SO_up_states(SO)+.25))-post_sleep_rhythms.so_delta.SO_up_states(SO),[-.25:0.001:.25])];
                        post_SO_dls = [post_SO_dls; histcounts(post_dls_unit((post_dls_unit>=post_sleep_rhythms.so_delta.SO_up_states(SO)-.25) & (post_dls_unit<=post_sleep_rhythms.so_delta.SO_up_states(SO)+.25))-post_sleep_rhythms.so_delta.SO_up_states(SO),[-.25:0.001:.25])];
                    end
                    
                    real_cc =  xcorr(reshape(post_SO_m1',size(post_SO_m1,1)*size(post_SO_m1,2),1),reshape(post_SO_dls',size(post_SO_dls,1)*size(post_SO_dls,2),1),100,'normalized');
                    all_pairs(pair_count).post_SO_cc = real_cc;
                    
                    shuffle_cc = zeros(100,201);
                    parfor shuffle = 1:100
                        post_SO_dls_shufffle = zeros(size(post_SO_dls,1),size(post_SO_dls,2));
                        for n = 1:size(post_SO_m1,1)
                            tmp_dls_spikes = find(post_SO_dls(n,:))+randi([-25 25],1,length(find(post_SO_dls(n,:))));
                            tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<size(post_SO_dls,2) & tmp_dls_spikes>0);
                            post_SO_dls_shufffle(n,tmp_dls_spikes) = 1;
                        end
                        shuffle_cc(shuffle,:) = xcorr(reshape(post_SO_m1',size(post_SO_m1,1)*size(post_SO_m1,2),1),reshape(post_SO_dls_shufffle',size(post_SO_dls_shufffle,1)*size(post_SO_dls_shufffle,2),1),100,'normalized')';
                    end
                    all_pairs(pair_count).post_SO_cc_shuffle = shuffle_cc;
                    
                    %%% PRE delta
                    
                    pre_delta_m1 = [];
                    pre_delta_dls = [];
                    for delta = 1:length(pre_sleep_rhythms.so_delta.delta_up_states)
                        pre_delta_m1 = [pre_delta_m1; histcounts(pre_m1_unit((pre_m1_unit>=pre_sleep_rhythms.so_delta.delta_up_states(delta)-.25) & (pre_m1_unit<=pre_sleep_rhythms.so_delta.delta_up_states(delta)+.25))-pre_sleep_rhythms.so_delta.delta_up_states(delta),[-.25:0.001:.25])];
                        pre_delta_dls = [pre_delta_dls; histcounts(pre_dls_unit((pre_dls_unit>=pre_sleep_rhythms.so_delta.delta_up_states(delta)-.25) & (pre_dls_unit<=pre_sleep_rhythms.so_delta.delta_up_states(delta)+.25))-pre_sleep_rhythms.so_delta.delta_up_states(delta),[-.25:0.001:.25])];
                    end
                    
                    real_cc =  xcorr(reshape(pre_delta_m1',size(pre_delta_m1,1)*size(pre_delta_m1,2),1),reshape(pre_delta_dls',size(pre_delta_dls,1)*size(pre_delta_dls,2),1),100,'normalized');
                    all_pairs(pair_count).pre_delta_cc = real_cc;
                    
                    shuffle_cc = zeros(100,201);
                    parfor shuffle = 1:100
                        pre_delta_dls_shufffle = zeros(size(pre_delta_dls,1),size(pre_delta_dls,2));
                        for n = 1:size(pre_delta_m1,1)
                            tmp_dls_spikes = find(pre_delta_dls(n,:))+randi([-25 25],1,length(find(pre_delta_dls(n,:))));
                            tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<size(pre_delta_dls,2) & tmp_dls_spikes>0);
                            pre_delta_dls_shufffle(n,tmp_dls_spikes) = 1;
                        end
                        shuffle_cc(shuffle,:) = xcorr(reshape(pre_delta_m1',size(pre_delta_m1,1)*size(pre_delta_m1,2),1),reshape(pre_delta_dls_shufffle',size(pre_delta_dls_shufffle,1)*size(pre_delta_dls_shufffle,2),1),100,'normalized')';
                    end
                    all_pairs(pair_count).pre_delta_cc_shuffle = shuffle_cc;
                    
                    %%% POST delta
                    
                    post_delta_m1 = [];
                    post_delta_dls = [];
                    for delta = 1:length(post_sleep_rhythms.so_delta.delta_up_states)
                        post_delta_m1 = [post_delta_m1; histcounts(post_m1_unit((post_m1_unit>=post_sleep_rhythms.so_delta.delta_up_states(delta)-.25) & (post_m1_unit<=post_sleep_rhythms.so_delta.delta_up_states(delta)+.25))-post_sleep_rhythms.so_delta.delta_up_states(delta),[-.25:0.001:.25])];
                        post_delta_dls = [post_delta_dls; histcounts(post_dls_unit((post_dls_unit>=post_sleep_rhythms.so_delta.delta_up_states(delta)-.25) & (post_dls_unit<=post_sleep_rhythms.so_delta.delta_up_states(delta)+.25))-post_sleep_rhythms.so_delta.delta_up_states(delta),[-.25:0.001:.25])];
                    end
                    
                    real_cc =  xcorr(reshape(post_delta_m1',size(post_delta_m1,1)*size(post_delta_m1,2),1),reshape(post_delta_dls',size(post_delta_dls,1)*size(post_delta_dls,2),1),100,'normalized');
                    all_pairs(pair_count).post_delta_cc = real_cc;
                    
                    shuffle_cc = zeros(100,201);
                    parfor shuffle = 1:100
                        post_delta_dls_shufffle = zeros(size(post_delta_dls,1),size(post_delta_dls,2));
                        for n = 1:size(post_delta_m1,1)
                            tmp_dls_spikes = find(post_delta_dls(n,:))+randi([-25 25],1,length(find(post_delta_dls(n,:))));
                            tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<size(post_delta_dls,2) & tmp_dls_spikes>0);
                            post_delta_dls_shufffle(n,tmp_dls_spikes) = 1;
                        end
                        shuffle_cc(shuffle,:) = xcorr(reshape(post_delta_m1',size(post_delta_m1,1)*size(post_delta_m1,2),1),reshape(post_delta_dls_shufffle',size(post_delta_dls_shufffle,1)*size(post_delta_dls_shufffle,2),1),100,'normalized')';
                    end
                    all_pairs(pair_count).post_delta_cc_shuffle = shuffle_cc;
                    
                end
                
                pair_count = pair_count + 1;
                
            end
        end
        
    end
end

end
