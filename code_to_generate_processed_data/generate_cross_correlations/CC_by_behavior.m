function [all_pairs] = CC_by_behavior(pre_m1_spiking, pre_dls_spiking, post_m1_spiking, post_dls_spiking, pre_beh_state, post_beh_state, pre_video_start_stop, post_video_start_stop,m1_spiking_width,dls_spiking_width,m1_spiking_peak2valley,dls_spiking_peak2valley);
    %   
    
    %%% DEFINE BEHAVIORAL STATE

        pre_nrem = interp1(1:length(pre_beh_state.nrem),pre_beh_state.nrem,linspace(1,length(pre_beh_state.nrem),length([0:0.001:pre_video_start_stop(2)-pre_video_start_stop(1)])));
%         pre_nrem_starts = find(diff(pre_nrem==1)>0.5);
%         pre_nrem_ends = find(diff(pre_nrem==1)<-0.5);
%         if pre_nrem(1)==1; pre_nrem_starts = [1 pre_nrem_starts]; end
%         if pre_nrem(end)==1; pre_nrem_ends = [pre_nrem_ends length(pre_nrem)]; end

        pre_rem = interp1(1:length(pre_beh_state.rem),pre_beh_state.rem,linspace(1,length(pre_beh_state.rem),length([0:0.001:pre_video_start_stop(2)-pre_video_start_stop(1)])));
%         pre_rem_starts = find(diff(pre_rem==1)>0.5);
%         pre_rem_ends = find(diff(pre_rem==1)<-0.5);
%         if pre_rem(1)==1; pre_rem_starts = [1 pre_rem_starts]; end
%         if pre_rem(end)==1; pre_rem_ends = [pre_rem_ends length(pre_rem)]; end

        pre_wake = interp1(1:length(pre_beh_state.wake),pre_beh_state.wake,linspace(1,length(pre_beh_state.wake),length([0:0.001:pre_video_start_stop(2)-pre_video_start_stop(1)])));
%         pre_wake_starts = find(diff(pre_wake==1)>0.5);
%         pre_wake_ends = find(diff(pre_wake==1)<-0.5);
%         if pre_wake(1)==1; pre_wake_starts = [1 pre_wake_starts]; end
%         if pre_wake(end)==1; pre_wake_ends = [pre_wake_ends length(pre_wake)]; end

        post_nrem = interp1(1:length(post_beh_state.nrem),post_beh_state.nrem,linspace(1,length(post_beh_state.nrem),length([0:0.001:post_video_start_stop(2)-post_video_start_stop(1)])));
%         post_nrem_starts = find(diff(post_nrem==1)>0.5);
%         post_nrem_ends = find(diff(post_nrem==1)<-0.5);
%         if post_nrem(1)==1; post_nrem_starts = [1 post_nrem_starts]; end
%         if post_nrem(end)==1; post_nrem_ends = [post_nrem_ends length(post_nrem)]; end

        post_rem = interp1(1:length(post_beh_state.rem),post_beh_state.rem,linspace(1,length(post_beh_state.rem),length([0:0.001:post_video_start_stop(2)-post_video_start_stop(1)])));
%         post_rem_starts = find(diff(post_rem==1)>0.5);
%         post_rem_ends = find(diff(post_rem==1)<-0.5);
%         if post_rem(1)==1; post_rem_starts = [1 post_rem_starts]; end
%         if post_rem(end)==1; post_rem_ends = [post_rem_ends length(post_rem)]; end

        post_wake = interp1(1:length(post_beh_state.wake),post_beh_state.wake,linspace(1,length(post_beh_state.wake),length([0:0.001:post_video_start_stop(2)-post_video_start_stop(1)])));
%         post_wake_starts = find(diff(post_wake==1)>0.5);
%         post_wake_ends = find(diff(post_wake==1)<-0.5);
%         if post_wake(1)==1; post_wake_starts = [1 post_wake_starts]; end
%         if post_wake(end)==1; post_wake_ends = [post_wake_ends length(post_wake)]; end

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
                    
                    %%% DETERMINE IF CONNECTED
                    
%                         m1_all_hist = [pre_m1_unit' post_m1_unit'];
%                         dls_all_hist = [pre_dls_unit' post_dls_unit'];
% 
%                         [CC, ~, CI] = crosscorr(m1_all_hist, dls_all_hist,'NumLags',100,'NumSTD',2);
% 
%                         all_pairs(pair_count).CC = CC;
%                         all_pairs(pair_count).CI = CI;
% 
%                         if mean(CC(102:111))>CI(1)
%                             all_pairs(pair_count).sig = 1;
%                         else
%                             all_pairs(pair_count).sig = 0;
%                         end
% 
%                         figure;
%                             hold on
%                             plot(CC)
%                             plot([0 200],[CI(1) CI(1)],'color','k')
%                             plot([101 101],[ylim],'color','k')
%                             xlim([1 200])
%                             title(['M1 chan ' num2str(m1_chan) ' unit ' num2str(m1_unit) ' DLS chan ' num2str(dls_chan) ' unit ' num2str(dls_unit) ' - SIG: ' num2str(all_pairs(pair_count).sig)])
                            
                    %%% OVERALL CC BY BEHAVIORAL STATE

                        pre_CC_nrem = xcorr(pre_m1_unit(pre_nrem==1),pre_dls_unit(pre_nrem==1),200,'normalized');
                        all_pairs(pair_count).pre_CC_nrem = pre_CC_nrem;

                        pre_CC_rem = xcorr(pre_m1_unit(pre_rem==1),pre_dls_unit(pre_rem==1),200,'normalized');
                        all_pairs(pair_count).pre_CC_rem = pre_CC_rem;

                        pre_CC_wake = xcorr(pre_m1_unit(pre_wake==1),pre_dls_unit(pre_wake==1),200,'normalized');
                        all_pairs(pair_count).pre_CC_wake = pre_CC_wake;

                        post_CC_nrem = xcorr(post_m1_unit(post_nrem==1),post_dls_unit(post_nrem==1),200,'normalized');
                        all_pairs(pair_count).post_CC_nrem = post_CC_nrem;

                        post_CC_rem = xcorr(post_m1_unit(post_rem==1),post_dls_unit(post_rem==1),200,'normalized');
                        all_pairs(pair_count).post_CC_rem = post_CC_rem;

                        post_CC_wake = xcorr(post_m1_unit(post_wake==1),post_dls_unit(post_wake==1),200,'normalized');
                        all_pairs(pair_count).post_CC_wake = post_CC_wake;
                    
                    %%% CC 1ST/2ND HALF BY BEHAVIORAL STATE
                    
                        pre_nrem_1 = find(pre_nrem==1);
                        pre_nrem_1 = pre_nrem_1(pre_nrem_1<floor(length(pre_nrem)/2));
                        all_pairs(pair_count).pre_nrem_CC_1 = xcorr(pre_m1_unit(pre_nrem_1),pre_dls_unit(pre_nrem_1),300,'normalized');

                        pre_nrem_2 = find(pre_nrem==1);
                        pre_nrem_2 = pre_nrem_2(pre_nrem_2>ceil(length(pre_nrem)/2));
                        all_pairs(pair_count).pre_nrem_CC_2 = xcorr(pre_m1_unit(pre_nrem_2),pre_dls_unit(pre_nrem_2),300,'normalized');

                        pre_rem_1 = find(pre_rem==1);
                        pre_rem_1 = pre_rem_1(pre_rem_1<floor(length(pre_rem)/2));
                        all_pairs(pair_count).pre_rem_CC_1 = xcorr(pre_m1_unit(pre_rem_1),pre_dls_unit(pre_rem_1),100,'normalized');

                        pre_rem_2 = find(pre_rem==1);
                        pre_rem_2 = pre_rem_2(pre_rem_2>ceil(length(pre_rem)/2));
                        all_pairs(pair_count).pre_rem_CC_2 = xcorr(pre_m1_unit(pre_rem_2),pre_dls_unit(pre_rem_2),100,'normalized');

                        pre_wake_1 = find(pre_wake==1);
                        pre_wake_1 = pre_wake_1(pre_wake_1<floor(length(pre_wake)/2));
                        all_pairs(pair_count).pre_wake_CC_1 = xcorr(pre_m1_unit(pre_wake_1),pre_dls_unit(pre_wake_1),100,'normalized');

                        pre_wake_2 = find(pre_wake==1);
                        pre_wake_2 = pre_wake_2(pre_wake_2>ceil(length(pre_wake)/2));
                        all_pairs(pair_count).pre_wake_CC_2 = xcorr(pre_m1_unit(pre_wake_2),pre_dls_unit(pre_wake_2),100,'normalized');

                        post_nrem_1 = find(post_nrem==1);
                        post_nrem_1 = post_nrem_1(post_nrem_1<floor(length(post_nrem)/2));
                        all_pairs(pair_count).post_nrem_CC_1 = xcorr(post_m1_unit(post_nrem_1),post_dls_unit(post_nrem_1),300,'normalized');

                        post_nrem_2 = find(post_nrem==1);
                        post_nrem_2 = post_nrem_2(post_nrem_2>ceil(length(post_nrem)/2));
                        all_pairs(pair_count).post_nrem_CC_2 = xcorr(post_m1_unit(post_nrem_2),post_dls_unit(post_nrem_2),300,'normalized');

                        post_rem_1 = find(post_rem==1);
                        post_rem_1 = post_rem_1(post_rem_1<floor(length(post_rem)/2));
                        all_pairs(pair_count).post_rem_CC_1 = xcorr(post_m1_unit(post_rem_1),post_dls_unit(post_rem_1),100,'normalized');

                        post_rem_2 = find(post_rem==1);
                        post_rem_2 = post_rem_2(post_rem_2>ceil(length(post_rem)/2));
                        all_pairs(pair_count).post_rem_CC_2 = xcorr(post_m1_unit(post_rem_2),post_dls_unit(post_rem_2),100,'normalized');

                        post_wake_1 = find(post_wake==1);
                        post_wake_1 = post_wake_1(post_wake_1<floor(length(post_wake)/2));
                        all_pairs(pair_count).post_wake_CC_1 = xcorr(post_m1_unit(post_wake_1),post_dls_unit(post_wake_1),100,'normalized');

                        post_wake_2 = find(post_wake==1);
                        post_wake_2 = post_wake_2(post_wake_2>ceil(length(post_wake)/2));
                        all_pairs(pair_count).post_wake_CC_2 = xcorr(post_m1_unit(post_wake_2),post_dls_unit(post_wake_2),100,'normalized');
                        
                    %%% CHANGE IN CC BY BEHAVIORAL STATE
                    
%                         pre_cc_nrem = [];
%                         for pre_nrem_bouts = 1:length(pre_nrem_starts)
%                             tmp_m1_unit = pre_m1_unit(pre_nrem_starts(pre_nrem_bouts):pre_nrem_ends(pre_nrem_bouts));
%                             tmp_dls_unit = pre_dls_unit(pre_nrem_starts(pre_nrem_bouts):pre_nrem_ends(pre_nrem_bouts));
%                             tmp_pre_cc_nrem = [];
%                             tmp_mins = ceil(length(tmp_m1_unit)/60000);
%                             if tmp_mins == 1
%                                 tmp_cc = xcorr(tmp_m1_unit,tmp_dls_unit,200,'normalized');
%                                 tmp_pre_cc_nrem = [tmp_pre_cc_nrem mean(tmp_cc(190:200))-mean(tmp_cc(1:50))];
%                             else
%                                 for mins = 1:tmp_mins
%                                     if mins==tmp_mins
%                                         tmp_cc = xcorr(tmp_m1_unit(1+((mins-1)*60000):end),tmp_dls_unit(1+((mins-1)*60000):end),200,'normalized');
%                                         tmp_pre_cc_nrem = [tmp_pre_cc_nrem mean(tmp_cc(190:200))-mean(tmp_cc(1:50))];
%                                     else
%                                         tmp_cc = xcorr(tmp_m1_unit(1+((mins-1)*60000):((mins)*60000)),tmp_dls_unit(1+((mins-1)*60000):((mins)*60000)),200,'normalized');
%                                         tmp_pre_cc_nrem = [tmp_pre_cc_nrem mean(tmp_cc(190:200))-mean(tmp_cc(1:50))];
%                                     end
%                                 end
%                             end
%                             pre_cc_nrem = [pre_cc_nrem; tmp_pre_cc_nrem nan(1,60-length(tmp_pre_cc_nrem))];
%                         end
%                         all_pairs(pair_count).pre_CC_nrem_change = nanmean(pre_cc_nrem);
%                     
%                         pre_cc_rem = [];
%                         for pre_rem_bouts = 1:length(pre_rem_starts)
%                             tmp_m1_unit = pre_m1_unit(pre_rem_starts(pre_rem_bouts):pre_rem_ends(pre_rem_bouts));
%                             tmp_dls_unit = pre_dls_unit(pre_rem_starts(pre_rem_bouts):pre_rem_ends(pre_rem_bouts));
%                             tmp_pre_cc_rem = [];
%                             tmp_mins = ceil(length(tmp_m1_unit)/60000);
%                             if tmp_mins == 1
%                                 tmp_cc = xcorr(tmp_m1_unit,tmp_dls_unit,200,'normalized');
%                                 tmp_pre_cc_rem = [tmp_pre_cc_rem mean(tmp_cc(190:200))-mean(tmp_cc(1:50))];
%                             else
%                                 for mins = 1:tmp_mins
%                                     if mins==tmp_mins
%                                         tmp_cc = xcorr(tmp_m1_unit(1+((mins-1)*60000):end),tmp_dls_unit(1+((mins-1)*60000):end),200,'normalized');
%                                         tmp_pre_cc_rem = [tmp_pre_cc_rem mean(tmp_cc(190:200))-mean(tmp_cc(1:50))];
%                                     else
%                                         tmp_cc = xcorr(tmp_m1_unit(1+((mins-1)*60000):((mins)*60000)),tmp_dls_unit(1+((mins-1)*60000):((mins)*60000)),200,'normalized');
%                                         tmp_pre_cc_rem = [tmp_pre_cc_rem mean(tmp_cc(190:200))-mean(tmp_cc(1:50))];
%                                     end
%                                 end
%                             end
%                             pre_cc_rem = [pre_cc_rem; tmp_pre_cc_rem nan(1,60-length(tmp_pre_cc_rem))];
%                         end
%                         all_pairs(pair_count).pre_CC_rem_change = nanmean(pre_cc_rem);
%                
%                         pre_cc_wake = [];
%                         for pre_wake_bouts = 1:length(pre_wake_starts)
%                             tmp_m1_unit = pre_m1_unit(pre_wake_starts(pre_wake_bouts):pre_wake_ends(pre_wake_bouts));
%                             tmp_dls_unit = pre_dls_unit(pre_wake_starts(pre_wake_bouts):pre_wake_ends(pre_wake_bouts));
%                             tmp_pre_cc_wake = [];
%                             tmp_mins = ceil(length(tmp_m1_unit)/60000);
%                             if tmp_mins == 1
%                                 tmp_cc = xcorr(tmp_m1_unit,tmp_dls_unit,200,'normalized');
%                                 tmp_pre_cc_wake = [tmp_pre_cc_wake mean(tmp_cc(190:200))-mean(tmp_cc(1:50))];
%                             else
%                                 for mins = 1:tmp_mins
%                                     if mins==tmp_mins
%                                         tmp_cc = xcorr(tmp_m1_unit(1+((mins-1)*60000):end),tmp_dls_unit(1+((mins-1)*60000):end),200,'normalized');
%                                         tmp_pre_cc_wake = [tmp_pre_cc_wake mean(tmp_cc(190:200))-mean(tmp_cc(1:50))];
%                                     else
%                                         tmp_cc = xcorr(tmp_m1_unit(1+((mins-1)*60000):((mins)*60000)),tmp_dls_unit(1+((mins-1)*60000):((mins)*60000)),200,'normalized');
%                                         tmp_pre_cc_wake = [tmp_pre_cc_wake mean(tmp_cc(190:200))-mean(tmp_cc(1:50))];
%                                     end
%                                 end
%                             end
%                             pre_cc_wake = [pre_cc_wake; tmp_pre_cc_wake nan(1,60-length(tmp_pre_cc_wake))];
%                         end
%                         all_pairs(pair_count).pre_CC_wake_change = nanmean(pre_cc_wake);
%                     
%                         post_cc_nrem = [];
%                         for post_nrem_bouts = 1:length(post_nrem_starts)
%                             tmp_m1_unit = post_m1_unit(post_nrem_starts(post_nrem_bouts):post_nrem_ends(post_nrem_bouts));
%                             tmp_dls_unit = post_dls_unit(post_nrem_starts(post_nrem_bouts):post_nrem_ends(post_nrem_bouts));
%                             tmp_post_cc_nrem = [];
%                             tmp_mins = ceil(length(tmp_m1_unit)/60000);
%                             if tmp_mins == 1
%                                 tmp_cc = xcorr(tmp_m1_unit,tmp_dls_unit,200,'normalized');
%                                 tmp_post_cc_nrem = [tmp_post_cc_nrem mean(tmp_cc(190:200))-mean(tmp_cc(1:50))];
%                             else
%                                 for mins = 1:tmp_mins
%                                     if mins==tmp_mins
%                                         tmp_cc = xcorr(tmp_m1_unit(1+((mins-1)*60000):end),tmp_dls_unit(1+((mins-1)*60000):end),200,'normalized');
%                                         tmp_post_cc_nrem = [tmp_post_cc_nrem mean(tmp_cc(190:200))-mean(tmp_cc(1:50))];
%                                     else
%                                         tmp_cc = xcorr(tmp_m1_unit(1+((mins-1)*60000):((mins)*60000)),tmp_dls_unit(1+((mins-1)*60000):((mins)*60000)),200,'normalized');
%                                         tmp_post_cc_nrem = [tmp_post_cc_nrem mean(tmp_cc(190:200))-mean(tmp_cc(1:50))];
%                                     end
%                                 end
%                             end
%                             post_cc_nrem = [post_cc_nrem; tmp_post_cc_nrem nan(1,60-length(tmp_post_cc_nrem))];
%                         end
%                         all_pairs(pair_count).post_CC_nrem_change = nanmean(post_cc_nrem);
% 
%                         post_cc_rem = [];
%                         for post_rem_bouts = 1:length(post_rem_starts)
%                             tmp_m1_unit = post_m1_unit(post_rem_starts(post_rem_bouts):post_rem_ends(post_rem_bouts));
%                             tmp_dls_unit = post_dls_unit(post_rem_starts(post_rem_bouts):post_rem_ends(post_rem_bouts));
%                             tmp_post_cc_rem = [];
%                             tmp_mins = ceil(length(tmp_m1_unit)/60000);
%                             if tmp_mins == 1
%                                 tmp_cc = xcorr(tmp_m1_unit,tmp_dls_unit,200,'normalized');
%                                 tmp_post_cc_rem = [tmp_post_cc_rem mean(tmp_cc(190:200))-mean(tmp_cc(1:50))];
%                             else
%                                 for mins = 1:tmp_mins
%                                     if mins==tmp_mins
%                                         tmp_cc = xcorr(tmp_m1_unit(1+((mins-1)*60000):end),tmp_dls_unit(1+((mins-1)*60000):end),200,'normalized');
%                                         tmp_post_cc_rem = [tmp_post_cc_rem mean(tmp_cc(190:200))-mean(tmp_cc(1:50))];
%                                     else
%                                         tmp_cc = xcorr(tmp_m1_unit(1+((mins-1)*60000):((mins)*60000)),tmp_dls_unit(1+((mins-1)*60000):((mins)*60000)),200,'normalized');
%                                         tmp_post_cc_rem = [tmp_post_cc_rem mean(tmp_cc(190:200))-mean(tmp_cc(1:50))];
%                                     end
%                                 end
%                             end
%                             post_cc_rem = [post_cc_rem; tmp_post_cc_rem nan(1,60-length(tmp_post_cc_rem))];
%                         end
%                         all_pairs(pair_count).post_CC_rem_change = nanmean(post_cc_rem);
% 
%                         post_cc_wake = [];
%                         for post_wake_bouts = 1:length(post_wake_starts)
%                             tmp_m1_unit = post_m1_unit(post_wake_starts(post_wake_bouts):post_wake_ends(post_wake_bouts));
%                             tmp_dls_unit = post_dls_unit(post_wake_starts(post_wake_bouts):post_wake_ends(post_wake_bouts));
%                             tmp_post_cc_wake = [];
%                             tmp_mins = ceil(length(tmp_m1_unit)/60000);
%                             if tmp_mins == 1
%                                 tmp_cc = xcorr(tmp_m1_unit,tmp_dls_unit,200,'normalized');
%                                 tmp_post_cc_wake = [tmp_post_cc_wake mean(tmp_cc(190:200))-mean(tmp_cc(1:50))];
%                             else
%                                 for mins = 1:tmp_mins
%                                     if mins==tmp_mins
%                                         tmp_cc = xcorr(tmp_m1_unit(1+((mins-1)*60000):end),tmp_dls_unit(1+((mins-1)*60000):end),200,'normalized');
%                                         tmp_post_cc_wake = [tmp_post_cc_wake mean(tmp_cc(190:200))-mean(tmp_cc(1:50))];
%                                     else
%                                         tmp_cc = xcorr(tmp_m1_unit(1+((mins-1)*60000):((mins)*60000)),tmp_dls_unit(1+((mins-1)*60000):((mins)*60000)),200,'normalized');
%                                         tmp_post_cc_wake = [tmp_post_cc_wake mean(tmp_cc(190:200))-mean(tmp_cc(1:50))];
%                                     end
%                                 end
%                             end
%                             post_cc_wake = [post_cc_wake; tmp_post_cc_wake nan(1,60-length(tmp_post_cc_wake))];
%                         end
%                         all_pairs(pair_count).post_CC_wake_change = nanmean(post_cc_wake);

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

                    %%% SPIKE PROPERTIES (OVERALL FR, SPIKE WIDTH, PEAK-TO-VALLEY)
                    
                        all_pairs(pair_count).M1_width = mean(m1_spiking_width{m1_chan,m1_unit});
                        all_pairs(pair_count).M1_peak2valley = mean(m1_spiking_peak2valley{m1_chan,m1_unit});
                        all_pairs(pair_count).DLS_width = mean(dls_spiking_width{dls_chan,dls_unit});
                        all_pairs(pair_count).DLS_peak2valley = mean(dls_spiking_peak2valley{dls_chan,dls_unit}); 

                        pair_count = pair_count + 1;

                end
            end
            
            pause(0.1);
            
        end
    end
    
end
            
        