function [all_pairs] = spindle_CC_before_after(pre_m1_spiking, pre_dls_spiking, post_m1_spiking, post_dls_spiking, pre_video_start_stop, post_video_start_stop,pre_sleep_rhythms,post_sleep_rhythms,sig)

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
                    
                    %%% PRE
                    
                    pre_spindle_cc = zeros(length(pre_sleep_rhythms.spindles{1,1}.pks),6);
                    pre_spindle_raster_m1 = zeros(length(pre_sleep_rhythms.spindles{1,1}.pks),30000);
                    pre_spindle_raster_dls = zeros(length(pre_sleep_rhythms.spindles{1,1}.pks),30000);
                    
                    for spindle = 1:length(pre_sleep_rhythms.spindles{1,1}.pks)
                        
                        pre_spindle_raster_m1(spindle,:) = histcounts(pre_m1_unit((pre_m1_unit>=pre_sleep_rhythms.spindles{1,1}.pks(spindle)-150) & (pre_m1_unit<=pre_sleep_rhythms.spindles{1,1}.pks(spindle)+150))-pre_sleep_rhythms.spindles{1,1}.pks(spindle),[-150:0.01:150]);
                        pre_spindle_raster_dls(spindle,:) = histcounts(pre_dls_unit((pre_dls_unit>=pre_sleep_rhythms.spindles{1,1}.pks(spindle)-150) & (pre_dls_unit<=pre_sleep_rhythms.spindles{1,1}.pks(spindle)+150))-pre_sleep_rhythms.spindles{1,1}.pks(spindle),[-150:0.01:150]);
                        
                        bin_count=1;
                        for bins = [-91 -61 -31 1 31 61]
                            
                            tmp_m1 = histcounts(pre_m1_unit((pre_m1_unit>=pre_sleep_rhythms.spindles{1,1}.pks(spindle)+bins) & (pre_m1_unit<=pre_sleep_rhythms.spindles{1,1}.pks(spindle)+(bins+30)))-pre_sleep_rhythms.spindles{1,1}.pks(spindle),[bins:0.001:(bins+30)]);
                            tmp_dls = histcounts(pre_dls_unit((pre_dls_unit>=pre_sleep_rhythms.spindles{1,1}.pks(spindle)+bins) & (pre_dls_unit<=pre_sleep_rhythms.spindles{1,1}.pks(spindle)+(bins+30)))-pre_sleep_rhythms.spindles{1,1}.pks(spindle),[bins:0.001:(bins+30)]);
                            
                            real_CC = xcorr(tmp_m1,tmp_dls,10,'normalized')';
                            
                            shuffle_CC = zeros(100,21);
                            parfor shuffle = 1:100
                                tmp_pre_dls_shuffle = zeros(length(tmp_dls),1);
                                tmp_dls_spikes = find(tmp_dls)+randi([-25 25],1,length(find(tmp_dls)));
                                tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<length(tmp_dls) & tmp_dls_spikes>0);
                                tmp_pre_dls_shuffle(tmp_dls_spikes) = 1;
                                tmp_CC = xcorr(tmp_m1,tmp_pre_dls_shuffle,10,'normalized')';
                                shuffle_CC(shuffle,:) = tmp_CC;
                            end
                            
                            tmp_cc = repmat(real_CC',100,1)-shuffle_CC;
                            pre_spindle_cc(spindle,bin_count) = mean(mean(tmp_cc(:,1:10)));
                            bin_count=bin_count+1;
                            
                        end
                    end
                    
                    all_pairs(pair_count).pre_spindle_cc = pre_spindle_cc;
                    all_pairs(pair_count).pre_spindle_raster_m1 = pre_spindle_raster_m1;
                    all_pairs(pair_count).pre_spindle_raster_dls = pre_spindle_raster_dls;
                    
                    %%% POST
                    
                    post_spindle_cc = zeros(length(post_sleep_rhythms.spindles{1,1}.pks),6);
                    post_spindle_raster_m1 = zeros(length(post_sleep_rhythms.spindles{1,1}.pks),30000);
                    post_spindle_raster_dls = zeros(length(post_sleep_rhythms.spindles{1,1}.pks),30000);
                    
                    for spindle = 1:length(post_sleep_rhythms.spindles{1,1}.pks)
                        
                        post_spindle_raster_m1(spindle,:) = histcounts(post_m1_unit((post_m1_unit>=post_sleep_rhythms.spindles{1,1}.pks(spindle)-150) & (post_m1_unit<=post_sleep_rhythms.spindles{1,1}.pks(spindle)+150))-post_sleep_rhythms.spindles{1,1}.pks(spindle),[-150:0.01:150]);
                        post_spindle_raster_dls(spindle,:) = histcounts(post_dls_unit((post_dls_unit>=post_sleep_rhythms.spindles{1,1}.pks(spindle)-150) & (post_dls_unit<=post_sleep_rhythms.spindles{1,1}.pks(spindle)+150))-post_sleep_rhythms.spindles{1,1}.pks(spindle),[-150:0.01:150]);
                        
                        bin_count=1;
                        for bins = [-91 -61 -31 1 31 61]
                            
                            tmp_m1 = histcounts(post_m1_unit((post_m1_unit>=post_sleep_rhythms.spindles{1,1}.pks(spindle)+bins) & (post_m1_unit<=post_sleep_rhythms.spindles{1,1}.pks(spindle)+(bins+30)))-post_sleep_rhythms.spindles{1,1}.pks(spindle),[bins:0.001:(bins+30)]);
                            tmp_dls = histcounts(post_dls_unit((post_dls_unit>=post_sleep_rhythms.spindles{1,1}.pks(spindle)+bins) & (post_dls_unit<=post_sleep_rhythms.spindles{1,1}.pks(spindle)+(bins+30)))-post_sleep_rhythms.spindles{1,1}.pks(spindle),[bins:0.001:(bins+30)]);
                            
                            real_CC = xcorr(tmp_m1,tmp_dls,10,'normalized')';
                            
                            shuffle_CC = zeros(100,21);
                            parfor shuffle = 1:100
                                tmp_post_dls_shuffle = zeros(length(tmp_dls),1);
                                tmp_dls_spikes = find(tmp_dls)+randi([-25 25],1,length(find(tmp_dls)));
                                tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<length(tmp_dls) & tmp_dls_spikes>0);
                                tmp_post_dls_shuffle(tmp_dls_spikes) = 1;
                                tmp_CC = xcorr(tmp_m1,tmp_post_dls_shuffle,10,'normalized')';
                                shuffle_CC(shuffle,:) = tmp_CC;
                            end
                            
                            tmp_cc = repmat(real_CC',100,1)-shuffle_CC;
                            post_spindle_cc(spindle,bin_count) = mean(mean(tmp_cc(:,1:10)));
                            bin_count=bin_count+1;
                            
                        end
                    end
                    
                    all_pairs(pair_count).post_spindle_cc = post_spindle_cc;
                    all_pairs(pair_count).post_spindle_raster_m1 = post_spindle_raster_m1;
                    all_pairs(pair_count).post_spindle_raster_dls = post_spindle_raster_dls;
                    
                end
                pair_count = pair_count + 1;
            end
        end
        
    end
end

end
