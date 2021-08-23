function [all_pairs] = jitter_test(pre_m1_spiking, pre_dls_spiking, post_m1_spiking, post_dls_spiking, pre_beh_state, post_beh_state, pre_video_start_stop, post_video_start_stop)
%

%%% DEFINE BEHAVIORAL STATE

pre_nrem = interp1(1:length(pre_beh_state.nrem),pre_beh_state.nrem,linspace(1,length(pre_beh_state.nrem),length([0:0.001:pre_video_start_stop(2)-pre_video_start_stop(1)])));
post_nrem = interp1(1:length(post_beh_state.nrem),post_beh_state.nrem,linspace(1,length(post_beh_state.nrem),length([0:0.001:post_video_start_stop(2)-post_video_start_stop(1)])));

%%% ALL PAIRS

pair_count = 1;

for m1_chan = 1:size(pre_m1_spiking,1)
    for m1_unit = 1:size(pre_m1_spiking,2)
        if isempty(pre_m1_spiking{m1_chan,m1_unit}) || isempty(post_m1_spiking{m1_chan,m1_unit})
            continue
        end
        
        pre_m1_unit = pre_m1_spiking{m1_chan,m1_unit}(pre_m1_spiking{m1_chan,m1_unit}>pre_video_start_stop(1) & pre_m1_spiking{m1_chan,m1_unit}<pre_video_start_stop(2))-pre_video_start_stop(1);
        pre_m1_unit = histc(pre_m1_unit,[0:0.001:pre_video_start_stop(2)-pre_video_start_stop(1)]);
        post_m1_unit = post_m1_spiking{m1_chan,m1_unit}(post_m1_spiking{m1_chan,m1_unit}>post_video_start_stop(1) & post_m1_spiking{m1_chan,m1_unit}<post_video_start_stop(2))-post_video_start_stop(1);
        post_m1_unit = histc(post_m1_unit,[0:0.001:post_video_start_stop(2)-post_video_start_stop(1)]);
        
        for dls_chan = 1:size(pre_dls_spiking,1)
            for dls_unit = 1:size(pre_dls_spiking,2)
                if isempty(pre_dls_spiking{dls_chan,dls_unit}) || isempty(post_dls_spiking{dls_chan,dls_unit})
                    continue
                end
                
                disp(['M1 CHAN ' num2str(m1_chan) ' UNIT ' num2str(m1_unit) ' AND DLS CHAN ' num2str(dls_chan) ' UNIT ' num2str(dls_unit)]);
                
                pre_dls_unit = pre_dls_spiking{dls_chan,dls_unit}(pre_dls_spiking{dls_chan,dls_unit}>pre_video_start_stop(1) & pre_dls_spiking{dls_chan,dls_unit}<pre_video_start_stop(2))-pre_video_start_stop(1);
                pre_dls_unit = histc(pre_dls_unit,[0:0.001:pre_video_start_stop(2)-pre_video_start_stop(1)]);
                post_dls_unit = post_dls_spiking{dls_chan,dls_unit}(post_dls_spiking{dls_chan,dls_unit}>post_video_start_stop(1) & post_dls_spiking{dls_chan,dls_unit}<post_video_start_stop(2))-post_video_start_stop(1);
                post_dls_unit = histc(post_dls_unit,[0:0.001:post_video_start_stop(2)-post_video_start_stop(1)]);

                %%% PAIR INFORMATION

                    all_pairs(pair_count).m1_chan_unit = [m1_chan m1_unit];
                    all_pairs(pair_count).dls_chan_unit = [dls_chan dls_unit];

                %%% JITTER TEST FOR SIGNIFICANCE

                    if sum(pre_nrem==1)<1000*60*5
                        tmp_pre_m1 = pre_m1_unit(pre_nrem==1);
                        tmp_pre_dls = pre_dls_unit(pre_nrem==1);
                    else
                        pre_index = find(pre_nrem==1);
                        pre_index = pre_index(1:1000*60*5);
                        tmp_pre_m1 = pre_m1_unit(pre_index);
                        tmp_pre_dls = pre_dls_unit(pre_index);
                    end

                    if sum(post_nrem==1)<1000*60*5
                        tmp_post_m1 = post_m1_unit(post_nrem==1);
                        tmp_post_dls = post_dls_unit(post_nrem==1);
                    else
                        post_index = find(post_nrem==1);
                        post_index = post_index(1:1000*60*5);
                        tmp_post_m1 = post_m1_unit(post_index);
                        tmp_post_dls = post_dls_unit(post_index);
                    end

                    tmp_m1 = [tmp_pre_m1' tmp_post_m1'];
                    tmp_dls = [tmp_pre_dls' tmp_post_dls'];

                    real_CC = xcorr(tmp_m1,tmp_dls,10,'normalized')';
                    shuffle_CC = zeros(1000,21);
                    parfor shuffle = 1:1000
                        tmp_pre_dls_shuffle = zeros(length(tmp_dls),1);
                        tmp_dls_spikes = find(tmp_dls)+randi([-25 25],1,length(find(tmp_dls)));
                        tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<length(tmp_dls) & tmp_dls_spikes>0);
                        tmp_pre_dls_shuffle(tmp_dls_spikes) = 1;
                        tmp_CC = xcorr(tmp_m1,tmp_pre_dls_shuffle,10,'normalized')';
                        shuffle_CC(shuffle,:) = tmp_CC;
                    end

                    all_pairs(pair_count).real_CC = real_CC;                    
                    all_pairs(pair_count).shuffle_CC = shuffle_CC;
                    all_pairs(pair_count).p_val = sum(mean(real_CC(1:10))>mean(shuffle_CC(:,1:10),2))/size(shuffle_CC,1);
                    pair_count = pair_count + 1;
                
            end
        end
    end
end

end

