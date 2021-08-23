function [all_pairs] = jitter_behavioral_state(pre_m1_spiking, pre_dls_spiking, post_m1_spiking, post_dls_spiking, pre_beh_state, post_beh_state, pre_video_start_stop, post_video_start_stop,sig)
%

%%% DEFINE BEHAVIORAL STATE

    pre_nrem = interp1(1:length(pre_beh_state.nrem),pre_beh_state.nrem,linspace(1,length(pre_beh_state.nrem),length([0:0.001:pre_video_start_stop(2)-pre_video_start_stop(1)])));
    pre_rem = interp1(1:length(pre_beh_state.rem),pre_beh_state.rem,linspace(1,length(pre_beh_state.rem),length([0:0.001:pre_video_start_stop(2)-pre_video_start_stop(1)])));
    pre_wake = interp1(1:length(pre_beh_state.wake),pre_beh_state.wake,linspace(1,length(pre_beh_state.wake),length([0:0.001:pre_video_start_stop(2)-pre_video_start_stop(1)])));

    post_nrem = interp1(1:length(post_beh_state.nrem),post_beh_state.nrem,linspace(1,length(post_beh_state.nrem),length([0:0.001:post_video_start_stop(2)-post_video_start_stop(1)])));
    post_rem = interp1(1:length(post_beh_state.rem),post_beh_state.rem,linspace(1,length(post_beh_state.rem),length([0:0.001:post_video_start_stop(2)-post_video_start_stop(1)])));
    post_wake = interp1(1:length(post_beh_state.wake),post_beh_state.wake,linspace(1,length(post_beh_state.wake),length([0:0.001:post_video_start_stop(2)-post_video_start_stop(1)])));

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

                    %%% CROSS CORRELATION BY BEHAVIORAL STATE

                    if sig(pair_count)>.99

                        %%% PRE NREM

                            pre_nrem_1 = find(pre_nrem==1);
                            pre_nrem_1 = pre_nrem_1(pre_nrem_1<floor(length(pre_nrem)/2));
                            tmp_pre_m1 = pre_m1_unit(pre_nrem_1);
                            tmp_pre_dls = pre_dls_unit(pre_nrem_1);

                            pre_CC_nrem = xcorr(tmp_pre_m1,tmp_pre_dls,100,'normalized')';
                            all_pairs(pair_count).pre_CC_nrem_1 = pre_CC_nrem;

                            shuffle_pre_CC_nrem = zeros(100,201);
                            parfor shuffle = 1:100
                                tmp_pre_dls_shuffle = zeros(length(tmp_pre_dls),1);
                                tmp_dls_spikes = find(tmp_pre_dls)+randi([-25 25],1,length(find(tmp_pre_dls)))';
                                tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<length(tmp_pre_dls) & tmp_dls_spikes>0);
                                tmp_pre_dls_shuffle(tmp_dls_spikes) = 1;
                                shuffle_pre_CC_nrem(shuffle,:) = xcorr(tmp_pre_m1,tmp_pre_dls_shuffle,100,'normalized')';
                            end
                            all_pairs(pair_count).shuffle_pre_CC_nrem_1 = shuffle_pre_CC_nrem;

                            pre_nrem_2 = find(pre_nrem==1);
                            pre_nrem_2 = pre_nrem_2(pre_nrem_2>ceil(length(pre_nrem)/2));
                            tmp_pre_m1 = pre_m1_unit(pre_nrem_2);
                            tmp_pre_dls = pre_dls_unit(pre_nrem_2);

                            pre_CC_nrem = xcorr(tmp_pre_m1,tmp_pre_dls,100,'normalized')';
                            all_pairs(pair_count).pre_CC_nrem_2 = pre_CC_nrem;

                            shuffle_pre_CC_nrem = zeros(100,201);
                            parfor shuffle = 1:100
                                tmp_pre_dls_shuffle = zeros(length(tmp_pre_dls),1);
                                tmp_dls_spikes = find(tmp_pre_dls)+randi([-25 25],1,length(find(tmp_pre_dls)))';
                                tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<length(tmp_pre_dls) & tmp_dls_spikes>0);
                                tmp_pre_dls_shuffle(tmp_dls_spikes) = 1;
                                shuffle_pre_CC_nrem(shuffle,:) = xcorr(tmp_pre_m1,tmp_pre_dls_shuffle,100,'normalized')';
                            end
                            all_pairs(pair_count).shuffle_pre_CC_nrem_2 = shuffle_pre_CC_nrem;

                        %%% PRE REM

                            pre_rem_1 = find(pre_rem==1);
                            pre_rem_1 = pre_rem_1(pre_rem_1<floor(length(pre_rem)/2));
                            tmp_pre_m1 = pre_m1_unit(pre_rem_1);
                            tmp_pre_dls = pre_dls_unit(pre_rem_1);

                            pre_CC_rem = xcorr(tmp_pre_m1,tmp_pre_dls,100,'normalized')';
                            all_pairs(pair_count).pre_CC_rem_1 = pre_CC_rem;

                            shuffle_pre_CC_rem = zeros(100,201);
                            parfor shuffle = 1:100
                                tmp_pre_dls_shuffle = zeros(length(tmp_pre_dls),1);
                                tmp_dls_spikes = find(tmp_pre_dls)+randi([-25 25],1,length(find(tmp_pre_dls)))';
                                tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<length(tmp_pre_dls) & tmp_dls_spikes>0);
                                tmp_pre_dls_shuffle(tmp_dls_spikes) = 1;
                                shuffle_pre_CC_rem(shuffle,:) = xcorr(tmp_pre_m1,tmp_pre_dls_shuffle,100,'normalized')';
                            end
                            all_pairs(pair_count).shuffle_pre_CC_rem_1 = shuffle_pre_CC_rem;

                            pre_rem_2 = find(pre_rem==1);
                            pre_rem_2 = pre_rem_2(pre_rem_2>ceil(length(pre_rem)/2));
                            tmp_pre_m1 = pre_m1_unit(pre_rem_2);
                            tmp_pre_dls = pre_dls_unit(pre_rem_2);

                            pre_CC_rem = xcorr(tmp_pre_m1,tmp_pre_dls,100,'normalized')';
                            all_pairs(pair_count).pre_CC_rem_2 = pre_CC_rem;

                            shuffle_pre_CC_rem = zeros(100,201);
                            parfor shuffle = 1:100
                                tmp_pre_dls_shuffle = zeros(length(tmp_pre_dls),1);
                                tmp_dls_spikes = find(tmp_pre_dls)+randi([-25 25],1,length(find(tmp_pre_dls)))';
                                tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<length(tmp_pre_dls) & tmp_dls_spikes>0);
                                tmp_pre_dls_shuffle(tmp_dls_spikes) = 1;
                                shuffle_pre_CC_rem(shuffle,:) = xcorr(tmp_pre_m1,tmp_pre_dls_shuffle,100,'normalized')';
                            end
                            all_pairs(pair_count).shuffle_pre_CC_rem_2 = shuffle_pre_CC_rem;

                        %%% PRE WAKE

                            pre_wake_1 = find(pre_wake==1);
                            pre_wake_1 = pre_wake_1(pre_wake_1<floor(length(pre_wake)/2));
                            tmp_pre_m1 = pre_m1_unit(pre_wake_1);
                            tmp_pre_dls = pre_dls_unit(pre_wake_1);

                            pre_CC_wake = xcorr(tmp_pre_m1,tmp_pre_dls,100,'normalized')';
                            all_pairs(pair_count).pre_CC_wake_1 = pre_CC_wake;

                            shuffle_pre_CC_wake = zeros(100,201);
                            parfor shuffle = 1:100
                                tmp_pre_dls_shuffle = zeros(length(tmp_pre_dls),1);
                                tmp_dls_spikes = find(tmp_pre_dls)+randi([-25 25],1,length(find(tmp_pre_dls)))';
                                tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<length(tmp_pre_dls) & tmp_dls_spikes>0);
                                tmp_pre_dls_shuffle(tmp_dls_spikes) = 1;
                                shuffle_pre_CC_wake(shuffle,:) = xcorr(tmp_pre_m1,tmp_pre_dls_shuffle,100,'normalized')';
                            end
                            all_pairs(pair_count).shuffle_pre_CC_wake_1 = shuffle_pre_CC_wake;

                            pre_wake_2 = find(pre_wake==1);
                            pre_wake_2 = pre_wake_2(pre_wake_2>ceil(length(pre_wake)/2));
                            tmp_pre_m1 = pre_m1_unit(pre_wake_2);
                            tmp_pre_dls = pre_dls_unit(pre_wake_2);

                            pre_CC_wake = xcorr(tmp_pre_m1,tmp_pre_dls,100,'normalized')';
                            all_pairs(pair_count).pre_CC_wake_2 = pre_CC_wake;

                            shuffle_pre_CC_wake = zeros(100,201);
                            parfor shuffle = 1:100
                                tmp_pre_dls_shuffle = zeros(length(tmp_pre_dls),1);
                                tmp_dls_spikes = find(tmp_pre_dls)+randi([-25 25],1,length(find(tmp_pre_dls)))';
                                tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<length(tmp_pre_dls) & tmp_dls_spikes>0);
                                tmp_pre_dls_shuffle(tmp_dls_spikes) = 1;
                                shuffle_pre_CC_wake(shuffle,:) = xcorr(tmp_pre_m1,tmp_pre_dls_shuffle,100,'normalized')';
                            end
                            all_pairs(pair_count).shuffle_pre_CC_wake_2 = shuffle_pre_CC_wake;

                        %%% POST NREM

                            post_nrem_1 = find(post_nrem==1);
                            post_nrem_1 = post_nrem_1(post_nrem_1<floor(length(post_nrem)/2));
                            tmp_post_m1 = post_m1_unit(post_nrem_1);
                            tmp_post_dls = post_dls_unit(post_nrem_1);

                            post_CC_nrem = xcorr(tmp_post_m1,tmp_post_dls,100,'normalized')';
                            all_pairs(pair_count).post_CC_nrem_1 = post_CC_nrem;

                            shuffle_post_CC_nrem = zeros(100,201);
                            parfor shuffle = 1:100
                                tmp_post_dls_shuffle = zeros(length(tmp_post_dls),1);
                                tmp_dls_spikes = find(tmp_post_dls)+randi([-25 25],1,length(find(tmp_post_dls)))';
                                tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<length(tmp_post_dls) & tmp_dls_spikes>0);
                                tmp_post_dls_shuffle(tmp_dls_spikes) = 1;
                                shuffle_post_CC_nrem(shuffle,:) = xcorr(tmp_post_m1,tmp_post_dls_shuffle,100,'normalized')';
                            end
                            all_pairs(pair_count).shuffle_post_CC_nrem_1 = shuffle_post_CC_nrem;

                            post_nrem_2 = find(post_nrem==1);
                            post_nrem_2 = post_nrem_2(post_nrem_2>ceil(length(post_nrem)/2));
                            tmp_post_m1 = post_m1_unit(post_nrem_2);
                            tmp_post_dls = post_dls_unit(post_nrem_2);

                            post_CC_nrem = xcorr(tmp_post_m1,tmp_post_dls,100,'normalized')';
                            all_pairs(pair_count).post_CC_nrem_2 = post_CC_nrem;

                            shuffle_post_CC_nrem = zeros(100,201);
                            parfor shuffle = 1:100
                                tmp_post_dls_shuffle = zeros(length(tmp_post_dls),1);
                                tmp_dls_spikes = find(tmp_post_dls)+randi([-25 25],1,length(find(tmp_post_dls)))';
                                tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<length(tmp_post_dls) & tmp_dls_spikes>0);
                                tmp_post_dls_shuffle(tmp_dls_spikes) = 1;
                                shuffle_post_CC_nrem(shuffle,:) = xcorr(tmp_post_m1,tmp_post_dls_shuffle,100,'normalized')';
                            end
                            all_pairs(pair_count).shuffle_post_CC_nrem_2 = shuffle_post_CC_nrem;

                        %%% POST REM

                            post_rem_1 = find(post_rem==1);
                            post_rem_1 = post_rem_1(post_rem_1<floor(length(post_rem)/2));
                            tmp_post_m1 = post_m1_unit(post_rem_1);
                            tmp_post_dls = post_dls_unit(post_rem_1);

                            post_CC_rem = xcorr(tmp_post_m1,tmp_post_dls,100,'normalized')';
                            all_pairs(pair_count).post_CC_rem_1 = post_CC_rem;

                            shuffle_post_CC_rem = zeros(100,201);
                            parfor shuffle = 1:100
                                tmp_post_dls_shuffle = zeros(length(tmp_post_dls),1);
                                tmp_dls_spikes = find(tmp_post_dls)+randi([-25 25],1,length(find(tmp_post_dls)))';
                                tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<length(tmp_post_dls) & tmp_dls_spikes>0);
                                tmp_post_dls_shuffle(tmp_dls_spikes) = 1;
                                shuffle_post_CC_rem(shuffle,:) = xcorr(tmp_post_m1,tmp_post_dls_shuffle,100,'normalized')';
                            end
                            all_pairs(pair_count).shuffle_post_CC_rem_1 = shuffle_post_CC_rem;

                            post_rem_2 = find(post_rem==1);
                            post_rem_2 = post_rem_2(post_rem_2>ceil(length(post_rem)/2));
                            tmp_post_m1 = post_m1_unit(post_rem_2);
                            tmp_post_dls = post_dls_unit(post_rem_2);

                            post_CC_rem = xcorr(tmp_post_m1,tmp_post_dls,100,'normalized')';
                            all_pairs(pair_count).post_CC_rem_2 = post_CC_rem;

                            shuffle_post_CC_rem = zeros(100,201);
                            parfor shuffle = 1:100
                                tmp_post_dls_shuffle = zeros(length(tmp_post_dls),1);
                                tmp_dls_spikes = find(tmp_post_dls)+randi([-25 25],1,length(find(tmp_post_dls)))';
                                tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<length(tmp_post_dls) & tmp_dls_spikes>0);
                                tmp_post_dls_shuffle(tmp_dls_spikes) = 1;
                                shuffle_post_CC_rem(shuffle,:) = xcorr(tmp_post_m1,tmp_post_dls_shuffle,100,'normalized')';
                            end
                            all_pairs(pair_count).shuffle_post_CC_rem_2 = shuffle_post_CC_rem;

                        %%% POST WAKE

                            post_wake_1 = find(post_wake==1);
                            post_wake_1 = post_wake_1(post_wake_1<floor(length(post_wake)/2));
                            tmp_post_m1 = post_m1_unit(post_wake_1);
                            tmp_post_dls = post_dls_unit(post_wake_1);

                            post_CC_wake = xcorr(tmp_post_m1,tmp_post_dls,100,'normalized')';
                            all_pairs(pair_count).post_CC_wake_1 = post_CC_wake;

                            shuffle_post_CC_wake = zeros(100,201);
                            parfor shuffle = 1:100
                                tmp_post_dls_shuffle = zeros(length(tmp_post_dls),1);
                                tmp_dls_spikes = find(tmp_post_dls)+randi([-25 25],1,length(find(tmp_post_dls)))';
                                tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<length(tmp_post_dls) & tmp_dls_spikes>0);
                                tmp_post_dls_shuffle(tmp_dls_spikes) = 1;
                                shuffle_post_CC_wake(shuffle,:) = xcorr(tmp_post_m1,tmp_post_dls_shuffle,100,'normalized')';
                            end
                            all_pairs(pair_count).shuffle_post_CC_wake_1 = shuffle_post_CC_wake;

                            post_wake_2 = find(post_wake==1);
                            post_wake_2 = post_wake_2(post_wake_2>ceil(length(post_wake)/2));
                            tmp_post_m1 = post_m1_unit(post_wake_2);
                            tmp_post_dls = post_dls_unit(post_wake_2);

                            post_CC_wake = xcorr(tmp_post_m1,tmp_post_dls,100,'normalized')';
                            all_pairs(pair_count).post_CC_wake_2 = post_CC_wake;

                            shuffle_post_CC_wake = zeros(100,201);
                            parfor shuffle = 1:100
                                tmp_post_dls_shuffle = zeros(length(tmp_post_dls),1);
                                tmp_dls_spikes = find(tmp_post_dls)+randi([-25 25],1,length(find(tmp_post_dls)))';
                                tmp_dls_spikes = tmp_dls_spikes(tmp_dls_spikes<length(tmp_post_dls) & tmp_dls_spikes>0);
                                tmp_post_dls_shuffle(tmp_dls_spikes) = 1;
                                shuffle_post_CC_wake(shuffle,:) = xcorr(tmp_post_m1,tmp_post_dls_shuffle,100,'normalized')';
                            end
                            all_pairs(pair_count).shuffle_post_CC_wake_2 = shuffle_post_CC_wake;

                        %%% PLOT

                            figure;

                                subplot(2,3,1); hold on;
                                plot(prctile(shuffle_pre_CC_nrem,99),'r')
                                plot(prctile(shuffle_pre_CC_nrem,1),'r')
                                plot(pre_CC_nrem,'b')
                                plot([101 101],[ylim],'k')
                                xlim([1 201])
                                title(['PRE NREM M1 CHAN ' num2str(m1_chan) ' UNIT ' num2str(m1_unit) ' DLS CHAN ' num2str(dls_chan) ' UNIT ' num2str(dls_unit)])

                                subplot(2,3,2); hold on;
                                plot(prctile(shuffle_pre_CC_rem,99),'r')
                                plot(prctile(shuffle_pre_CC_rem,1),'r')
                                plot(pre_CC_rem,'b')
                                plot([101 101],[ylim],'k')
                                xlim([1 201])
                                title(['PRE REM M1 CHAN ' num2str(m1_chan) ' UNIT ' num2str(m1_unit) ' DLS CHAN ' num2str(dls_chan) ' UNIT ' num2str(dls_unit)])

                                subplot(2,3,3); hold on;
                                plot(prctile(shuffle_pre_CC_wake,99),'r')
                                plot(prctile(shuffle_pre_CC_wake,1),'r')
                                plot(pre_CC_wake,'b')
                                plot([101 101],[ylim],'k')
                                xlim([1 201])
                                title(['PRE WAKE M1 CHAN ' num2str(m1_chan) ' UNIT ' num2str(m1_unit) ' DLS CHAN ' num2str(dls_chan) ' UNIT ' num2str(dls_unit)])

                                subplot(2,3,4); hold on;
                                plot(prctile(shuffle_post_CC_nrem,99),'r')
                                plot(prctile(shuffle_post_CC_nrem,1),'r')
                                plot(post_CC_nrem,'b')
                                plot([101 101],[ylim],'k')
                                xlim([1 201])
                                title(['POST NREM M1 CHAN ' num2str(m1_chan) ' UNIT ' num2str(m1_unit) ' DLS CHAN ' num2str(dls_chan) ' UNIT ' num2str(dls_unit)])

                                subplot(2,3,5); hold on;
                                plot(prctile(shuffle_post_CC_rem,99),'r')
                                plot(prctile(shuffle_post_CC_rem,1),'r')
                                plot(post_CC_rem,'b')
                                plot([101 101],[ylim],'k')
                                xlim([1 201])
                                title(['POST REM M1 CHAN ' num2str(m1_chan) ' UNIT ' num2str(m1_unit) ' DLS CHAN ' num2str(dls_chan) ' UNIT ' num2str(dls_unit)])

                                subplot(2,3,6); hold on;
                                plot(prctile(shuffle_post_CC_wake,99),'r')
                                plot(prctile(shuffle_post_CC_wake,1),'r')
                                plot(post_CC_wake,'b')
                                plot([101 101],[ylim],'k')
                                xlim([1 201])
                                title(['POST WAKE M1 CHAN ' num2str(m1_chan) ' UNIT ' num2str(m1_unit) ' DLS CHAN ' num2str(dls_chan) ' UNIT ' num2str(dls_unit)])

                        pause(0.1)

                    end

                    %%% FR BY BEH STATE

                    all_pairs(pair_count).pre_M1_FR_nrem = sum(pre_m1_unit(pre_nrem==1))/(sum(pre_nrem==1)/1000);
                    all_pairs(pair_count).pre_M1_FR_rem = sum(pre_m1_unit(pre_rem==1))/(sum(pre_rem==1)/1000);
                    all_pairs(pair_count).pre_M1_FR_wake = sum(pre_m1_unit(pre_wake==1))/(sum(pre_wake==1)/1000);

                    all_pairs(pair_count).post_M1_FR_nrem = sum(post_m1_unit(post_nrem==1))/(sum(post_nrem==1)/1000);
                    all_pairs(pair_count).post_M1_FR_rem = sum(post_m1_unit(post_rem==1))/(sum(post_rem==1)/1000);
                    all_pairs(pair_count).post_M1_FR_wake = sum(post_m1_unit(post_wake==1))/(sum(post_wake==1)/1000);

                    all_pairs(pair_count).pre_DLS_FR_nrem = sum(pre_dls_unit(pre_nrem==1))/(sum(pre_nrem==1)/1000);
                    all_pairs(pair_count).pre_DLS_FR_rem = sum(pre_dls_unit(pre_rem==1))/(sum(pre_rem==1)/1000);
                    all_pairs(pair_count).pre_DLS_FR_wake = sum(pre_dls_unit(pre_wake==1))/(sum(pre_wake==1)/1000);

                    all_pairs(pair_count).post_DLS_FR_nrem = sum(post_dls_unit(post_nrem==1))/(sum(post_nrem==1)/1000);
                    all_pairs(pair_count).post_DLS_FR_rem = sum(post_dls_unit(post_rem==1))/(sum(post_rem==1)/1000);
                    all_pairs(pair_count).post_DLS_FR_wake = sum(post_dls_unit(post_wake==1))/(sum(post_wake==1)/1000);

                    pair_count = pair_count + 1;

                end
            end
        end
    end
    
end

