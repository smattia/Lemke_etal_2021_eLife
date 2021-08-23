%% SCRIPT TO GENERATE FIGURE 5 & FIGURE 5 FIGURE SUPPLEMENTS
% 
% to generate figures, download data from https://doi.org/10.7272/Q6KK9927

clc; clear all; close all;
data_path = 'D:\Lemke_etal_2021\Lemke_eLife_2021_data';

%% FIG 5 C | EXAMPLE CC COMPARISON OF WAVE VS. NREM 

    animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
    animal_days = {[1:8],[1:6],[2:10],[4:17],[4:13],[1:5 7:9]};

    for animal = 5%1:length(animal_id)

        CC_path = [data_path '\' animal_id{animal} '\cross_correlations'];

        nrem_cc = [];
        rem_cc = [];
        wake_cc = [];
        for day=animal_days{animal}
            load([CC_path '\day_' num2str(day) '\behavioral_state_cc.mat']);
            if isfield(all_pairs,'pre_CC_nrem_1')
                for n=1:length(all_pairs)
                    if ~isempty(all_pairs(n).pre_CC_nrem_1)
                        nrem_cc = [nrem_cc; mean([mean(repmat(all_pairs(n).pre_CC_nrem_1,100,1)-all_pairs(n).shuffle_pre_CC_nrem_1); ...
                            mean(repmat(all_pairs(n).pre_CC_nrem_2,100,1)-all_pairs(n).shuffle_pre_CC_nrem_2); ...
                            mean(repmat(all_pairs(n).post_CC_nrem_1,100,1)-all_pairs(n).shuffle_post_CC_nrem_1); ...
                            mean(repmat(all_pairs(n).post_CC_nrem_2,100,1)-all_pairs(n).shuffle_post_CC_nrem_2)])];
                        rem_cc = [rem_cc; mean([mean(repmat(all_pairs(n).pre_CC_rem_1,100,1)-all_pairs(n).shuffle_pre_CC_rem_1); ...
                            mean(repmat(all_pairs(n).pre_CC_rem_2,100,1)-all_pairs(n).shuffle_pre_CC_rem_2); ...
                            mean(repmat(all_pairs(n).post_CC_rem_1,100,1)-all_pairs(n).shuffle_post_CC_rem_1); ...
                            mean(repmat(all_pairs(n).post_CC_rem_2,100,1)-all_pairs(n).shuffle_post_CC_rem_2)])];
                        wake_cc = [wake_cc; mean([mean(repmat(all_pairs(n).pre_CC_wake_1,100,1)-all_pairs(n).shuffle_pre_CC_wake_1); ...
                            mean(repmat(all_pairs(n).pre_CC_wake_2,100,1)-all_pairs(n).shuffle_pre_CC_wake_2); ...
                            mean(repmat(all_pairs(n).post_CC_wake_1,100,1)-all_pairs(n).shuffle_post_CC_wake_1); ...
                            mean(repmat(all_pairs(n).post_CC_wake_2,100,1)-all_pairs(n).shuffle_post_CC_wake_2)])];
                    end
                end
            end
        end

        figure;
            subplot(2,1,1); hold on;
                plot(mean(nrem_cc))
                plot([101 101],[-3e-4 7e-4],'k')
                xlim([51 150])
                ylim([-3e-4 7e-4])
                title([animal_id{animal} ' day ' num2str(day) '  pair ' num2str(n)])
            subplot(2,1,2); hold on;
                plot(mean(wake_cc))
                plot([101 101],[-3e-4 7e-4],'k')
                xlim([51 150])
                ylim([-3e-4 7e-4])

    end

%% FIG 5 C | ALL CC COMPARISON OF WAKE VS. NREM

    animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
    animal_days = {[1:8],[1:6],[2:10],[4:17],[4:13],[1:5 7:9]};

    all_nrem_cc = [];
    all_rem_cc = [];
    all_wake_cc = [];

    for animal = 1:length(animal_id)

        CC_path = [data_path '\' animal_id{animal} '\cross_correlations'];

        nrem_cc = [];
        rem_cc = [];
        wake_cc = [];
        for day=animal_days{animal}
            load([CC_path '\day_' num2str(day) '\behavioral_state_cc.mat']);
            if isfield(all_pairs,'pre_CC_nrem_1')
                for n=1:length(all_pairs)
                    if ~isempty(all_pairs(n).pre_CC_nrem_1)
                        nrem_cc = [nrem_cc; mean([mean(repmat(all_pairs(n).pre_CC_nrem_1,100,1)-all_pairs(n).shuffle_pre_CC_nrem_1); ...
                            mean(repmat(all_pairs(n).pre_CC_nrem_2,100,1)-all_pairs(n).shuffle_pre_CC_nrem_2); ...
                            mean(repmat(all_pairs(n).post_CC_nrem_1,100,1)-all_pairs(n).shuffle_post_CC_nrem_1); ...
                            mean(repmat(all_pairs(n).post_CC_nrem_2,100,1)-all_pairs(n).shuffle_post_CC_nrem_2)])];
                        rem_cc = [rem_cc; mean([mean(repmat(all_pairs(n).pre_CC_rem_1,100,1)-all_pairs(n).shuffle_pre_CC_rem_1); ...
                            mean(repmat(all_pairs(n).pre_CC_rem_2,100,1)-all_pairs(n).shuffle_pre_CC_rem_2); ...
                            mean(repmat(all_pairs(n).post_CC_rem_1,100,1)-all_pairs(n).shuffle_post_CC_rem_1); ...
                            mean(repmat(all_pairs(n).post_CC_rem_2,100,1)-all_pairs(n).shuffle_post_CC_rem_2)])];
                        wake_cc = [wake_cc; mean([mean(repmat(all_pairs(n).pre_CC_wake_1,100,1)-all_pairs(n).shuffle_pre_CC_wake_1); ...
                            mean(repmat(all_pairs(n).pre_CC_wake_2,100,1)-all_pairs(n).shuffle_pre_CC_wake_2); ...
                            mean(repmat(all_pairs(n).post_CC_wake_1,100,1)-all_pairs(n).shuffle_post_CC_wake_1); ...
                            mean(repmat(all_pairs(n).post_CC_wake_2,100,1)-all_pairs(n).shuffle_post_CC_wake_2)])];
                    end
                end
            end
        end

        all_nrem_cc = [all_nrem_cc; nrem_cc];
        all_rem_cc = [all_rem_cc; rem_cc];
        all_wake_cc = [all_wake_cc; wake_cc];

        figure; hold on;
            plot(nanmean(nrem_cc),'b');
            plot(nanmean(rem_cc),'r');
            plot(nanmean(wake_cc),'k');
            line([101 101],[-5e-4 15e-4],'color','k')
            xlim([26 176])
            ylim([-6e-4 10e-4])
            title(num2str(size(nrem_cc,1)))

    end
    
%%% PLOT ALL

    figure;
        hold on;
        histogram(mean(all_nrem_cc(:,86:100),2)-mean(all_wake_cc(:,86:100),2),'DisplayStyle','Stairs')
        xlim([-1e-3 1e-3])
        title(num2str(sum(mean(all_nrem_cc(:,86:100),2)-mean(all_wake_cc(:,86:100),2)>0)/size(all_nrem_cc,1)))

    figure;
        subplot(1,2,1)
            hold on
            plot(nanmean(all_nrem_cc)+nanstd(all_nrem_cc)/sqrt(size(all_nrem_cc,1)),'b');
            plot(nanmean(all_nrem_cc)-nanstd(all_nrem_cc)/sqrt(size(all_nrem_cc,1)),'b');
            line([101 101],[-3e-4 10e-4],'color','k')
            xlim([26 176])
            ylim([-3e-4 8e-4])
        subplot(1,2,2)
            hold on
            plot(nanmean(all_wake_cc)+nanstd(all_wake_cc)/sqrt(size(all_wake_cc,1)),'k');
            plot(nanmean(all_wake_cc)-nanstd(all_wake_cc)/sqrt(size(all_wake_cc,1)),'k');
            line([101 101],[-3e-4 10e-4],'color','k')
            xlim([26 176])
            ylim([-3e-4 8e-4])

%% FIG 5 F | EXAMPLE CC COMPARISON OF NREM RHYTHMS

    animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
    animal_days = {[1:8],[1:6],[2:10],[4:17],[4:13],[1:5 7:9]};

    for animal = 3%:length(animal_id)

        CC_path = [data_path '\' animal_id{animal} '\cross_correlations'];

        spindle_cc = [];
        SO_cc = [];
        delta_cc = [];
        for day=animal_days{animal}
            load([CC_path '\day_' num2str(day) '\NREM_rhythms_cc.mat']);
            if isfield(all_pairs,'pre_spindle_cc')
                for n=1:length(all_pairs)
                    if ~isempty(all_pairs(n).pre_spindle_cc)
                        spindle_cc = [spindle_cc; mean(repmat(all_pairs(n).pre_spindle_cc',100,1)-all_pairs(n).pre_spindle_cc_shuffle); mean(repmat(all_pairs(n).post_spindle_cc',100,1)-all_pairs(n).post_spindle_cc_shuffle)];
                        SO_cc = [SO_cc; mean(repmat(all_pairs(n).pre_SO_cc',100,1)-all_pairs(n).pre_SO_cc_shuffle); mean(repmat(all_pairs(n).post_SO_cc',100,1)-all_pairs(n).post_SO_cc_shuffle)];
                        delta_cc = [delta_cc; mean(repmat(all_pairs(n).pre_delta_cc',100,1)-all_pairs(n).pre_delta_cc_shuffle); mean(repmat(all_pairs(n).post_delta_cc',100,1)-all_pairs(n).post_delta_cc_shuffle)];
                    end
                end
            end
        end

        figure;
            subplot(3,1,1); hold on;
                plot(mean(spindle_cc))
                plot([101 101],[ylim],'k')
                xlim([51 150]) 
                ylim([-1e-3 2e-3])
                title(animal_id{animal})                                   
            subplot(3,1,2); hold on;
                plot(mean(SO_cc))
                plot([101 101],[ylim],'k')
                xlim([51 150]) 
                ylim([-1e-3 2e-3])
            subplot(3,1,3); hold on;
                plot(mean(delta_cc))
                plot([101 101],[ylim],'k')
                xlim([51 150]) 
                ylim([-1e-3 2e-3])
                pause(0.1);

    end

%% FIG 5 F | ALL CC COMPARISON OF NREM RHYTHMS

    animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
    animal_days = {[1:8],[1:6],[2:10],[4:17],[4:13],[1:5 7:9]};

    all_spindle_cc = [];
    all_SO_cc = [];
    all_delta_cc = [];

    for animal = 1:length(animal_id)

        CC_path = [data_path '\' animal_id{animal} '\cross_correlations'];

        spindle_cc = [];
        SO_cc = [];
        delta_cc = [];
        for day=animal_days{animal}
            load([CC_path '\day_' num2str(day) '\NREM_rhythms_cc.mat']);
            if isfield(all_pairs,'pre_spindle_cc')
                for n=1:length(all_pairs)
                    if ~isempty(all_pairs(n).pre_spindle_cc)
                        spindle_cc = [spindle_cc; mean(repmat(all_pairs(n).pre_spindle_cc',100,1)-all_pairs(n).pre_spindle_cc_shuffle)];
                        SO_cc = [SO_cc; mean(repmat(all_pairs(n).pre_SO_cc',100,1)-all_pairs(n).pre_SO_cc_shuffle)];
                        delta_cc = [delta_cc; mean(repmat(all_pairs(n).pre_delta_cc',100,1)-all_pairs(n).pre_delta_cc_shuffle)];
                    end
                end
            end
        end

        all_spindle_cc = [all_spindle_cc; spindle_cc];
        all_SO_cc = [all_SO_cc; SO_cc];
        all_delta_cc = [all_delta_cc; delta_cc];

        figure;
        hold on
            plot(nanmean(spindle_cc),'b');
            plot(nanmean(SO_cc),'r');
            plot(nanmean(delta_cc),'k');
            line([101 101],[-3e-3 3e-3],'color','k')
            xlim([1 201])
            ylim([-3e-3 3e-3])
            title(num2str(size(spindle_cc,1)))

    end

    %%% PLOT ALL

    figure;
        hold on;
        histogram(mean(all_spindle_cc(:,86:100),2)-mean(all_SO_cc(:,86:100),2),[-6e-3:.5e-3:6e-3],'DisplayStyle','Stairs')
        title(num2str(sum(mean(all_spindle_cc(:,86:100),2)-mean(all_SO_cc(:,86:100),2)>0)/size(all_spindle_cc,1)))    
        xlim([-6e-3 6e-3])

    figure;
        hold on;
        histogram(mean(all_spindle_cc(:,86:100),2)-mean(all_delta_cc(:,86:100),2),[-6e-3:.5e-3:6e-3],'DisplayStyle','Stairs')
        title(num2str(sum(mean(all_spindle_cc(:,86:100),2)-mean(all_delta_cc(:,86:100),2)>0)/size(all_spindle_cc,1)))        
        xlim([-6e-3 6e-3])

    figure;
        subplot(1,3,1)
            hold on
            plot(nanmean(all_spindle_cc)+nanstd(all_spindle_cc)/sqrt(size(all_spindle_cc,1)),'b');
            plot(nanmean(all_spindle_cc)-nanstd(all_spindle_cc)/sqrt(size(all_spindle_cc,1)),'b');
            line([101 101],[-1e-3 2e-3],'color','k')
            xlim([26 176])
            ylim([-1e-3 2e-3])
        subplot(1,3,2)
            hold on
            plot(nanmean(all_SO_cc)+nanstd(all_SO_cc)/sqrt(size(all_SO_cc,1)),'r');
            plot(nanmean(all_SO_cc)-nanstd(all_SO_cc)/sqrt(size(all_SO_cc,1)),'r');
            line([101 101],[-1e-3 2e-3],'color','k')
            xlim([26 176])
            ylim([-1e-3 2e-3])
        subplot(1,3,3)
            hold on
            plot(nanmean(all_delta_cc)+nanstd(all_delta_cc)/sqrt(size(all_delta_cc,1)),'g');
            plot(nanmean(all_delta_cc)-nanstd(all_delta_cc)/sqrt(size(all_delta_cc,1)),'g');
            line([101 101],[-1e-3 2e-3],'color','k')
            xlim([26 176])
            ylim([-1e-3 2e-3])        

%% sFIG 1 | DLS RECORDING MSN VS. FSI

    animal_id = {'LC3','LC4','LC6','LC5','LC2'};
    animal_days = {[1:8],[1:6],[2:10],[4:17],[4:13]};

    %%% ALL ANIMAL PERC
    all_spike_width = [];
    all_spike_amp = [];
    perc_msn = [];

    for animal = 1:length(animal_id)

        spike_width = [];
        spike_FR = [];
        spike_amp = [];

        for day = animal_days{animal}
            load([data_path '\' animal_id{animal} '\cross_correlations\day_' num2str(day) '\CC_by_behavior.mat']);
            load([data_path '\' animal_id{animal} '\cross_correlations\day_' num2str(day) '\CC_by_sleep_rhythm.mat'],'all_dls_spindle');

            for unit = 1:length(all_dls_spindle)
                spike_width = [spike_width all_pairs(unit).DLS_width];
                spike_amp = [spike_amp all_pairs(unit).DLS_peak2valley];
                spike_FR = [spike_FR mean([all_pairs(unit).pre_DLS_FR_wake all_pairs(unit).pre_DLS_FR_nrem all_pairs(unit).pre_DLS_FR_rem all_pairs(unit).post_DLS_FR_wake all_pairs(unit).post_DLS_FR_nrem all_pairs(unit).post_DLS_FR_rem])]; 
            end

        end

        perc_msn = [perc_msn; length(spike_width(spike_width>10))/length(spike_width) 1-length(spike_width(spike_width>10))/length(spike_width)];
        all_spike_width = [all_spike_width spike_width];
        all_spike_amp = [all_spike_amp spike_amp];

    end

%%% NEED RAW DATA NOT ON DRYAD
% spike_path = {['E:\T398\neural_data\kilosort\2.26.19'], ...
%             ['E:\T398\neural_data\kilosort\2.27.19'], ...
%             ['E:\T398\neural_data\kilosort\2.28.19'], ...
%             ['E:\T398\neural_data\kilosort\3.1.19'], ...
%             ['E:\T398\neural_data\kilosort\3.2.19'], ...
%             ['E:\T398\neural_data\kilosort\3.3.19'], ...
%             ['E:\T398\neural_data\kilosort\3.4.19'], ...
%             ['E:\T398\neural_data\kilosort\3.5.19'], ...
%             ['E:\T398\neural_data\kilosort\3.6.19']};
% 
% spike_width = [];
% spike_amp = [];
%     
% for day = [1:5 7:9]
% 
%     opts = delimitedTextImportOptions("NumVariables", 2);
%     opts.DataLines = [2, Inf];
%     opts.Delimiter = "\t";
%     opts.VariableNames = ["cluster_id", "group"];
%     opts.VariableTypes = ["double", "categorical"];
%     opts = setvaropts(opts, 2, "EmptyFieldRule", "auto");
%     opts.ExtraColumnsRule = "ignore";
%     opts.EmptyLineRule = "read";
% 
%     cd([spike_path{day} '\DLS']);
%     load('rez2.mat')        
%     dls_ts = readNPY('spike_times.npy');
%     dls_id = readNPY('spike_clusters.npy');
%     dls_unit_class = readtable('cluster_group.tsv', opts);
%     dls_clust_id = zeros(1,max(dls_unit_class{:,1}));
%     dls_clust_id_idx = dls_unit_class{:,1}'+1;
%     idx_count = 1;
%     for clust_id = dls_clust_id_idx
%         tmp_qual = cellstr(dls_unit_class{idx_count,2});
%         if strcmp(tmp_qual{1,1},'noise')
%             dls_clust_id(clust_id) = 0;
%         else
%             dls_clust_id(clust_id) = 1;
%         end
%         idx_count = idx_count + 1;
%     end
%     dls_clust_id = find(dls_clust_id)-1;
%     idx = ismember(dls_id,dls_clust_id);       
%     
%     gwfparams.dataDir = [spike_path{day} '\DLS\'];     % KiloSort/Phy output folder
%     gwfparams.fileName = 'DLS.bin';                     % .dat file containing the raw 
%     gwfparams.dataType = 'int16';                       % Data type of .dat file (this should be BP filtered)
%     gwfparams.nCh = 38;                                 % Number of channels that were streamed to disk in .dat file
%     gwfparams.wfWin = [-40 41];                         % Number of samples before and after spiketime to include in waveform
%     gwfparams.nWf = 2000;                               % Number of waveforms per unit to pull out
%     gwfparams.spikeTimes =    dls_ts(idx);                   % Vector of cluster spike times (in samples) same length as .spikeClusters
%     gwfparams.spikeClusters = dls_id(idx);                   % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
%     
%     wf = getWaveForms(gwfparams)
% 
%     for unit = 1:length(wf.unitIDs)
%         
%         amp = [];
%         for n = 1:size(wf.waveFormsMean,2)
%             amp = [amp max(wf.waveFormsMean(unit,n,:))-min(wf.waveFormsMean(unit,n,:))];
%         end
%         [~, i] = max(amp);
%         tmp_waveform = squeeze(wf.waveFormsMean(unit,i,:));
%         [min_val, min_idx] = min(tmp_waveform);
%         [max_val, max_idx] = max(tmp_waveform);
% 
%         spike_width = [spike_width max_idx-min_idx];
%         spike_amp = [spike_amp max_val-min_val];
%         all_spike_width = [all_spike_width max_idx-min_idx];
%         all_spike_amp = [all_spike_amp max_val-min_val];
%         
%     end
%          
% end
% 
% perc_msn = [perc_msn; length(spike_width(spike_width>12))/length(spike_width) 1-length(spike_width(spike_width>12))/length(spike_width)];

perc_msn_overall = [length(all_spike_width(all_spike_width>12))/length(all_spike_width) 1-length(all_spike_width(all_spike_width>12))/length(all_spike_width)];
figure
bar([perc_msn; perc_msn_overall],'stacked')

%%% EXAMPLE CLUSTERING

animal = 5;

spike_width = [];
spike_FR = [];
spike_amp = [];

for day = animal_days{animal}
    load([data_path '\' animal_id{5} '\cross_correlations\day_' num2str(day) '\CC_by_behavior.mat']);
    load([data_path '\' animal_id{5} '\cross_correlations\day_' num2str(day) '\CC_by_sleep_rhythm.mat'],'all_dls_spindle');

    for unit = 1:length(all_dls_spindle)
        spike_width = [spike_width all_pairs(unit).DLS_width];
        spike_amp = [spike_amp all_pairs(unit).DLS_peak2valley];
        spike_FR = [spike_FR mean([all_pairs(unit).pre_DLS_FR_wake all_pairs(unit).pre_DLS_FR_nrem all_pairs(unit).pre_DLS_FR_rem all_pairs(unit).post_DLS_FR_wake all_pairs(unit).post_DLS_FR_nrem all_pairs(unit).post_DLS_FR_rem])]; 
    end

end

figure;
    subplot(2,2,2);
        hold on;
        scatter(spike_width(spike_width>9),spike_amp(spike_width>9),30,[0 0 0],'filled');
        scatter(spike_width(spike_width<9),spike_amp(spike_width<9),30,[1 0 0],'filled');
        xlim([0 16])
        ylim([20 120])
    subplot(2,2,4);
        hold on;
        histogram(spike_width(spike_width>9),[0:1:16],'DisplayStyle','Stairs','EdgeColor',[0 0 0])
        histogram(spike_width(spike_width<9),[0:1:16],'DisplayStyle','Stairs','EdgeColor',[1 0 0])
        xlim([0 16])
        ylim([0 50])
    subplot(2,2,1)
        hold on;        
        h1 = histogram(spike_amp(spike_width>9),[20:5:120],'DisplayStyle','Stairs','EdgeColor',[0 0 0])
        h2 = histogram(spike_amp(spike_width<9),[20:5:120],'DisplayStyle','Stairs','EdgeColor',[1 0 0])
        h1.Orientation = 'horizontal';
        h2.Orientation = 'horizontal';
        xlim([0 30])            
        ylim([20 120])

%% FIG 5 B & sFIG 2 | EXAMPLE JITTER CROSS CORRELATION

animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
animal_days = {[1:8],[1:6],[2:10],[4:17],[4:13],[1:5 7:9]};

for animal = 3%1:length(animal_id)
    
    CC_path = [data_path '\' animal_id{animal} '\cross_correlations'];
    
    for day=5%animal_days{animal}
        
        load([CC_path '\day_' num2str(day) '\behavioral_state_cc.mat']);
        
        for n=145%:length(all_pairs)

            figure;
                subplot(3,1,1); hold on;
                    plot(mean([all_pairs(n).pre_CC_nrem_1; all_pairs(n).pre_CC_nrem_2; all_pairs(n).post_CC_nrem_1; all_pairs(n).post_CC_nrem_2]))        
                    title([animal_id{animal} ' day ' num2str(day) '  pair ' num2str(n) ' - RAW CC'])
                    plot([101 101],[0 6e-3],'color','k')                            
                    xlim([51 150])
                    ylim([0 6e-3])
                subplot(3,1,2); hold on;
                    plot(prctile((all_pairs(n).shuffle_pre_CC_nrem_1+all_pairs(n).shuffle_pre_CC_nrem_2+all_pairs(n).shuffle_post_CC_nrem_1+all_pairs(n).shuffle_post_CC_nrem_2)./4,50))
                    title([animal_id{animal} ' day ' num2str(day) '  pair ' num2str(n) ' - MEAN JITTER CC'])
                    plot([101 101],[0 6e-3],'color','k')                            
                    xlim([51 150])
                    ylim([0 6e-3])
                subplot(3,1,3); hold on;
                    plot(mean([mean(repmat(all_pairs(n).pre_CC_nrem_1,100,1)-all_pairs(n).shuffle_pre_CC_nrem_1); ...
                                              mean(repmat(all_pairs(n).pre_CC_nrem_2,100,1)-all_pairs(n).shuffle_pre_CC_nrem_2); ...
                                              mean(repmat(all_pairs(n).post_CC_nrem_1,100,1)-all_pairs(n).shuffle_post_CC_nrem_1); ...
                                              mean(repmat(all_pairs(n).post_CC_nrem_2,100,1)-all_pairs(n).shuffle_post_CC_nrem_2)]));
                    title([animal_id{animal} ' day ' num2str(day) '  pair ' num2str(n) ' - RAW - MEAN JITTER CC'])
                    plot([101 101],[-2e-3 4e-3],'color','k')                            
                    xlim([51 150])
                    ylim([-2e-3 4e-3])

            figure;
                hold on;
                plot(mean([all_pairs(n).pre_CC_nrem_1; all_pairs(n).pre_CC_nrem_2; all_pairs(n).post_CC_nrem_1; all_pairs(n).post_CC_nrem_2]),'color','k')        
                plot(prctile((all_pairs(n).shuffle_pre_CC_nrem_1+all_pairs(n).shuffle_pre_CC_nrem_2+all_pairs(n).shuffle_post_CC_nrem_1+all_pairs(n).shuffle_post_CC_nrem_2)./4,99),'color',[1 0 0])
                plot(prctile((all_pairs(n).shuffle_pre_CC_nrem_1+all_pairs(n).shuffle_pre_CC_nrem_2+all_pairs(n).shuffle_post_CC_nrem_1+all_pairs(n).shuffle_post_CC_nrem_2)./4,95),'color',[5/6 0 1/6])
                plot(prctile((all_pairs(n).shuffle_pre_CC_nrem_1+all_pairs(n).shuffle_pre_CC_nrem_2+all_pairs(n).shuffle_post_CC_nrem_1+all_pairs(n).shuffle_post_CC_nrem_2)./4,90),'color',[4/6 0 2/6])
                plot(prctile((all_pairs(n).shuffle_pre_CC_nrem_1+all_pairs(n).shuffle_pre_CC_nrem_2+all_pairs(n).shuffle_post_CC_nrem_1+all_pairs(n).shuffle_post_CC_nrem_2)./4,80),'color',[3/6 0 3/6])
                plot(prctile((all_pairs(n).shuffle_pre_CC_nrem_1+all_pairs(n).shuffle_pre_CC_nrem_2+all_pairs(n).shuffle_post_CC_nrem_1+all_pairs(n).shuffle_post_CC_nrem_2)./4,70),'color',[2/6 0 4/6])
                plot(prctile((all_pairs(n).shuffle_pre_CC_nrem_1+all_pairs(n).shuffle_pre_CC_nrem_2+all_pairs(n).shuffle_post_CC_nrem_1+all_pairs(n).shuffle_post_CC_nrem_2)./4,60),'color',[1/6 0 5/6])
                plot(prctile((all_pairs(n).shuffle_pre_CC_nrem_1+all_pairs(n).shuffle_pre_CC_nrem_2+all_pairs(n).shuffle_post_CC_nrem_1+all_pairs(n).shuffle_post_CC_nrem_2)./4,50),'color',[0 0 1])
                xlim([51 150])
                
        end
    end
end
   
%% sFIG 2 | JITTER PERCENTILES

animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
animal_days = {[1:8],[1:6],[2:10],[4:17],[4:13],[1:5 7:9]};

p_val = [];
for animal = 1:length(animal_id)
    
    CC_path = [data_path '\' animal_id{animal} '\cross_correlations'];
    for day=animal_days{animal}
        load([CC_path '\day_' num2str(day) '\jitter_test.mat']);
        for pair = 1:length(all_pairs)
            p_val = [p_val all_pairs(pair).p_val];
        end
    end
end

figure;
    subplot(1,2,1)
        histogram(p_val,[0:0.01:1],'DisplayStyle','Stairs')
    subplot(1,2,2)
        hold on;
        histogram(p_val,[.9:0.001:1],'DisplayStyle','Stairs')
        plot([.99 .99],[0 100],'k')
        title(['Sig pairs: ' num2str(sum(p_val>.99)) ' - ' num2str(sum(p_val>.99)/length(p_val)) ' percent of pairs'])

%% sFIG 3 | FIRING RATE

animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
animal_days = {[1:8],[1:6],[2:10],[4:17],[4:13],[1:5 7:9]};

all_m1_pre_nrem_fr = [];
all_m1_pre_rem_fr = [];
all_m1_pre_wake_fr = [];
all_m1_post_nrem_fr = [];
all_m1_post_rem_fr = [];
all_m1_post_wake_fr = [];
all_m1_nrem_fr = [];
all_m1_rem_fr = [];
all_m1_wake_fr = [];

all_dls_pre_nrem_fr = [];
all_dls_pre_rem_fr = [];
all_dls_pre_wake_fr = [];
all_dls_post_nrem_fr = [];
all_dls_post_rem_fr = [];
all_dls_post_wake_fr = [];
all_dls_nrem_fr = [];
all_dls_rem_fr = [];
all_dls_wake_fr = [];

for animal = 1:length(animal_id)
    
    m1_pre_nrem_fr = [];
    m1_pre_rem_fr = [];
    m1_pre_wake_fr = [];
    m1_post_nrem_fr = [];
    m1_post_rem_fr = [];
    m1_post_wake_fr = [];
    
    dls_pre_nrem_fr = [];
    dls_pre_rem_fr = [];
    dls_pre_wake_fr = [];
    dls_post_nrem_fr = [];
    dls_post_rem_fr = [];
    dls_post_wake_fr = [];
    
    for day = animal_days{animal}
        
        load([data_path '\' animal_id{animal} '\cross_correlations\day_' num2str(day) '\CC_by_behavior.mat']);
        
        m1_chans = [];
        current_chan = [0 0];
        for n=1:length(all_pairs)
            if current_chan(1)==all_pairs(n).m1_chan_unit(1) && current_chan(2)==all_pairs(n).m1_chan_unit(2)
                continue
            else
                m1_pre_nrem_fr = [m1_pre_nrem_fr all_pairs(n).pre_M1_FR_nrem];
                m1_pre_rem_fr = [m1_pre_rem_fr all_pairs(n).pre_M1_FR_rem];
                m1_pre_wake_fr = [m1_pre_wake_fr all_pairs(n).pre_M1_FR_wake];
                m1_post_nrem_fr = [m1_post_nrem_fr all_pairs(n).post_M1_FR_nrem];
                m1_post_rem_fr = [m1_post_rem_fr all_pairs(n).post_M1_FR_rem];
                m1_post_wake_fr = [m1_post_wake_fr all_pairs(n).post_M1_FR_wake];
                current_chan = [all_pairs(n).m1_chan_unit(1) all_pairs(n).m1_chan_unit(2)];
                m1_chans = [m1_chans n];
            end
        end
        
        for n=1:(m1_chans(2)-1)
            dls_pre_nrem_fr = [dls_pre_nrem_fr all_pairs(n).pre_DLS_FR_nrem];
            dls_pre_rem_fr = [dls_pre_rem_fr all_pairs(n).pre_DLS_FR_rem];
            dls_pre_wake_fr = [dls_pre_wake_fr all_pairs(n).pre_DLS_FR_wake];
            dls_post_nrem_fr = [dls_post_nrem_fr all_pairs(n).post_DLS_FR_nrem];
            dls_post_rem_fr = [dls_post_rem_fr all_pairs(n).post_DLS_FR_rem];
            dls_post_wake_fr = [dls_post_wake_fr all_pairs(n).post_DLS_FR_wake];
        end
        
    end
    
    figure;
    hold on
    errorbar(1,nanmean([m1_pre_nrem_fr m1_post_nrem_fr]),nanstd([m1_pre_nrem_fr m1_post_nrem_fr])/sqrt(length([m1_pre_nrem_fr m1_post_nrem_fr])),'color','k');
    errorbar(2,nanmean([m1_pre_rem_fr m1_post_rem_fr]),nanstd([m1_pre_rem_fr m1_post_rem_fr])/sqrt(length([m1_pre_rem_fr m1_post_rem_fr])),'color','k');
    errorbar(3,nanmean([m1_pre_wake_fr m1_post_wake_fr]),nanstd([m1_pre_wake_fr m1_post_wake_fr])/sqrt(length([m1_pre_wake_fr m1_post_wake_fr])),'color','k');
    errorbar(5,nanmean([dls_pre_nrem_fr dls_post_nrem_fr]),nanstd([dls_pre_nrem_fr dls_post_nrem_fr])/sqrt(length([dls_pre_nrem_fr dls_post_nrem_fr])),'color','r');
    errorbar(6,nanmean([dls_pre_rem_fr dls_post_rem_fr]),nanstd([dls_pre_rem_fr dls_post_rem_fr])/sqrt(length([dls_pre_rem_fr dls_post_rem_fr])),'color','r');
    errorbar(7,nanmean([dls_pre_wake_fr dls_post_wake_fr]),nanstd([dls_pre_wake_fr dls_post_wake_fr])/sqrt(length([dls_pre_wake_fr dls_post_wake_fr])),'color','r');
    title([animal_id{animal} ' (NREM/REM/WAKE)'])
    xlim([0 8])
    
    all_m1_pre_nrem_fr = [all_m1_pre_nrem_fr m1_pre_nrem_fr];
    all_m1_pre_rem_fr = [all_m1_pre_rem_fr m1_pre_rem_fr];
    all_m1_pre_wake_fr = [all_m1_pre_wake_fr m1_pre_wake_fr];
    all_m1_post_nrem_fr = [all_m1_post_nrem_fr m1_post_nrem_fr];
    all_m1_post_rem_fr = [all_m1_post_rem_fr m1_post_rem_fr];
    all_m1_post_wake_fr = [all_m1_post_wake_fr m1_post_wake_fr];
    
    all_dls_pre_nrem_fr = [all_dls_pre_nrem_fr dls_pre_nrem_fr];
    all_dls_pre_rem_fr = [all_dls_pre_rem_fr dls_pre_rem_fr];
    all_dls_pre_wake_fr = [all_dls_pre_wake_fr dls_pre_wake_fr];
    all_dls_post_nrem_fr = [all_dls_post_nrem_fr dls_post_nrem_fr];
    all_dls_post_rem_fr = [all_dls_post_rem_fr dls_post_rem_fr];
    all_dls_post_wake_fr = [all_dls_post_wake_fr dls_post_wake_fr];
    
    all_m1_nrem_fr = [all_m1_nrem_fr mean([m1_pre_nrem_fr; m1_post_nrem_fr])];
    all_m1_rem_fr = [all_m1_rem_fr mean([m1_pre_rem_fr; m1_post_rem_fr])];
    all_m1_wake_fr = [all_m1_wake_fr mean([m1_pre_wake_fr; m1_post_wake_fr])];
    
    all_dls_nrem_fr = [all_dls_nrem_fr mean([dls_pre_nrem_fr; dls_post_nrem_fr])];
    all_dls_rem_fr = [all_dls_rem_fr mean([dls_pre_rem_fr; dls_post_rem_fr])];
    all_dls_wake_fr = [all_dls_wake_fr mean([dls_pre_wake_fr; dls_post_wake_fr])];
    
end

%%% PLOT 

    fig = figure;
    
        subplot(2,4,1)
            hold on
            histogram(all_m1_pre_nrem_fr,[1:1:30],'DisplayStyle','Stairs','Normalization','Probability')
            histogram(all_m1_post_nrem_fr,[1:1:30],'DisplayStyle','Stairs','Normalization','Probability')
            [h p] = kstest2(all_m1_pre_nrem_fr,all_m1_post_nrem_fr)
            title(['m1 nrem - ks test p: ' num2str(p)])
            ylim([0 .3])

        subplot(2,4,2)
            hold on
            histogram(all_m1_pre_rem_fr,[1:1:30],'DisplayStyle','Stairs','Normalization','Probability')
            histogram(all_m1_post_rem_fr,[1:1:30],'DisplayStyle','Stairs','Normalization','Probability')
            [h p] = kstest2(all_m1_pre_rem_fr,all_m1_post_rem_fr)
            title(['m1 rem - ks test p: ' num2str(p)])
            ylim([0 .3])

        subplot(2,4,3)
            hold on
            histogram(all_m1_pre_wake_fr,[1:1:30],'DisplayStyle','Stairs','Normalization','Probability')
            histogram(all_m1_post_wake_fr,[1:1:30],'DisplayStyle','Stairs','Normalization','Probability')
            [h p] = kstest2(all_m1_pre_wake_fr,all_m1_post_wake_fr)
            title(['m1 wake - ks test p: ' num2str(p)])
            ylim([0 .3])

        subplot(2,4,4) 
            hold on;
        
            errorbar(1,nanmean(all_m1_pre_nrem_fr),nanstd(all_m1_pre_nrem_fr)/sqrt(length(all_m1_pre_nrem_fr)),'color','r');
            errorbar(2,nanmean(all_m1_post_nrem_fr),nanstd(all_m1_post_nrem_fr)/sqrt(length(all_m1_post_nrem_fr)),'color','r');
            
            errorbar(4,nanmean(all_m1_pre_rem_fr),nanstd(all_m1_pre_rem_fr)/sqrt(length(all_m1_pre_rem_fr)),'color','r');
            errorbar(5,nanmean(all_m1_post_rem_fr),nanstd(all_m1_post_rem_fr)/sqrt(length(all_m1_post_rem_fr)),'color','r');
            
            errorbar(7,nanmean(all_m1_pre_wake_fr),nanstd(all_m1_pre_wake_fr)/sqrt(length(all_m1_pre_wake_fr)),'color','r');
            errorbar(8,nanmean(all_m1_post_wake_fr),nanstd(all_m1_post_wake_fr)/sqrt(length(all_m1_post_wake_fr)),'color','r');
            
            [h p_nrem ci stats_nrem] = ttest(all_m1_pre_nrem_fr,all_m1_post_nrem_fr);
            [h p_rem ci stats_rem] = ttest(all_m1_pre_rem_fr,all_m1_post_rem_fr);
            [h p_wake ci stats_wake] = ttest(all_m1_pre_wake_fr,all_m1_post_wake_fr);
            
            title(['M1 - NREM p=' num2str(p_nrem) ' - REM p=' num2str(p_rem) ' - WAKE p=' num2str(p_wake)]);
            xlim([0 9])
            ylim([3 5])
            
            disp(['M1 pre NREM units : ' num2str(sum(~isnan(all_m1_pre_nrem_fr)))])
            disp(['M1 post NREM units : ' num2str(sum(~isnan(all_m1_post_nrem_fr)))])     
            disp(['NREM: pre ' num2str(nanmean(all_m1_pre_nrem_fr)) '+-' num2str(nanstd(all_m1_pre_nrem_fr)/sqrt(length(all_m1_pre_nrem_fr)))]);
            disp(['NREM: post ' num2str(nanmean(all_m1_post_nrem_fr)) '+-' num2str(nanstd(all_m1_post_nrem_fr)/sqrt(length(all_m1_post_nrem_fr)))]);
            disp(['NREM t-test: t(' num2str(stats_nrem.df) ') = ' num2str(stats_nrem.tstat) ', P=' num2str(num2str(nanmean(p_nrem)))]);
            
            disp(['M1 pre REM units : ' num2str(sum(~isnan(all_m1_pre_rem_fr)))])
            disp(['M1 post REM units : ' num2str(sum(~isnan(all_m1_post_rem_fr)))])    
            disp(['REM: pre ' num2str(nanmean(all_m1_pre_rem_fr)) '+-' num2str(nanstd(all_m1_pre_rem_fr)/sqrt(length(all_m1_pre_rem_fr)))])
            disp(['REM: post ' num2str(nanmean(all_m1_post_rem_fr)) '+-' num2str(nanstd(all_m1_post_rem_fr)/sqrt(length(all_m1_post_rem_fr)))])
            disp(['REM t-test: t(' num2str(stats_rem.df) ') = ' num2str(stats_rem.tstat) ', P=' num2str(num2str(nanmean(p_rem)))]);
            
            disp(['M1 pre wake units : ' num2str(sum(~isnan(all_m1_pre_wake_fr)))])
            disp(['M1 post wake units : ' num2str(sum(~isnan(all_m1_post_wake_fr)))])    
            disp(['wake: pre ' num2str(nanmean(all_m1_pre_wake_fr)) '+-' num2str(nanstd(all_m1_pre_wake_fr)/sqrt(length(all_m1_pre_wake_fr)))])
            disp(['wake: post ' num2str(nanmean(all_m1_post_wake_fr)) '+-' num2str(nanstd(all_m1_post_wake_fr)/sqrt(length(all_m1_post_wake_fr)))])
            disp(['wake t-test: t(' num2str(stats_wake.df) ') = ' num2str(stats_wake.tstat) ', P=' num2str(num2str(nanmean(p_wake)))]);
            
        subplot(2,4,5)
            hold on
            histogram(all_dls_pre_nrem_fr,[1:1:30],'DisplayStyle','Stairs','Normalization','Probability')
            histogram(all_dls_post_nrem_fr,[1:1:30],'DisplayStyle','Stairs','Normalization','Probability')
            [h p] = kstest2(all_dls_pre_nrem_fr,all_dls_post_nrem_fr)
            title(['dls nrem - ks test p: ' num2str(p)])
            ylim([0 .3])

        subplot(2,4,6)
            hold on
            histogram(all_dls_pre_rem_fr,[1:1:30],'DisplayStyle','Stairs','Normalization','Probability')
            histogram(all_dls_post_rem_fr,[1:1:30],'DisplayStyle','Stairs','Normalization','Probability')
            [h p] = kstest2(all_dls_pre_rem_fr,all_dls_post_rem_fr)
            title(['dls rem - ks test p: ' num2str(p)])
            ylim([0 .3])

        subplot(2,4,7)
            hold on
            histogram(all_dls_pre_wake_fr,[1:1:30],'DisplayStyle','Stairs','Normalization','Probability')
            histogram(all_dls_post_wake_fr,[1:1:30],'DisplayStyle','Stairs','Normalization','Probability')
            [h p] = kstest2(all_dls_pre_wake_fr,all_dls_post_wake_fr)
            title(['dls wake - ks test p: ' num2str(p)])
            ylim([0 .3])

        subplot(2,4,8) 
            hold on;
        
            errorbar(1,nanmean(all_dls_pre_nrem_fr),nanstd(all_dls_pre_nrem_fr)/sqrt(length(all_dls_pre_nrem_fr)),'color','b');
            errorbar(2,nanmean(all_dls_post_nrem_fr),nanstd(all_dls_post_nrem_fr)/sqrt(length(all_dls_post_nrem_fr)),'color','b');
            errorbar(4,nanmean(all_dls_pre_rem_fr),nanstd(all_dls_pre_rem_fr)/sqrt(length(all_dls_pre_rem_fr)),'color','b');
            errorbar(5,nanmean(all_dls_post_rem_fr),nanstd(all_dls_post_rem_fr)/sqrt(length(all_dls_post_rem_fr)),'color','b');
            errorbar(7,nanmean(all_dls_pre_wake_fr),nanstd(all_dls_pre_wake_fr)/sqrt(length(all_dls_pre_wake_fr)),'color','b');
            errorbar(8,nanmean(all_dls_post_wake_fr),nanstd(all_dls_post_wake_fr)/sqrt(length(all_dls_post_wake_fr)),'color','b');
            
            [h p_nrem ci stats_nrem] = ttest(all_dls_pre_nrem_fr,all_dls_post_nrem_fr);
            [h p_rem ci stats_rem] = ttest(all_dls_pre_rem_fr,all_dls_post_rem_fr);
            [h p_wake ci stats_wake] = ttest(all_dls_pre_wake_fr,all_dls_post_wake_fr);
            
            title(['DLS - NREM p=' num2str(p_nrem) ' - REM p=' num2str(p_rem) ' - WAKE p=' num2str(p_wake)])
            xlim([0 9])
            ylim([1.5 3])
            
            disp(['DLS pre NREM units : ' num2str(sum(~isnan(all_dls_pre_nrem_fr)))])
            disp(['DLS post NREM units : ' num2str(sum(~isnan(all_dls_post_nrem_fr)))])            
            disp(['NREM: pre ' num2str(nanmean(all_dls_pre_nrem_fr)) '+-' num2str(nanstd(all_dls_pre_nrem_fr)/sqrt(length(all_dls_pre_nrem_fr)))]);
            disp(['NREM: post ' num2str(nanmean(all_dls_post_nrem_fr)) '+-' num2str(nanstd(all_dls_post_nrem_fr)/sqrt(length(all_dls_post_nrem_fr)))]);
            disp(['NREM t-test: t(' num2str(stats_nrem.df) ') = ' num2str(stats_nrem.tstat) ', P=' num2str(num2str(nanmean(p_nrem)))]);
            
            disp(['DLS pre REM units : ' num2str(sum(~isnan(all_dls_pre_rem_fr)))])
            disp(['DLS post REM units : ' num2str(sum(~isnan(all_dls_post_rem_fr)))])    
            disp(['REM: pre ' num2str(nanmean(all_dls_pre_rem_fr)) '+-' num2str(nanstd(all_dls_pre_rem_fr)/sqrt(length(all_dls_pre_rem_fr)))])
            disp(['REM: post ' num2str(nanmean(all_dls_post_rem_fr)) '+-' num2str(nanstd(all_dls_post_rem_fr)/sqrt(length(all_dls_post_rem_fr)))])
            disp(['REM t-test: t(' num2str(stats_rem.df) ') = ' num2str(stats_rem.tstat) ', P=' num2str(num2str(nanmean(p_rem)))]);
            
            disp(['DLS pre wake units : ' num2str(sum(~isnan(all_dls_pre_wake_fr)))])
            disp(['DLS post wake units : ' num2str(sum(~isnan(all_dls_post_wake_fr)))])    
            disp(['wake: pre ' num2str(nanmean(all_dls_pre_wake_fr)) '+-' num2str(nanstd(all_dls_pre_wake_fr)/sqrt(length(all_dls_pre_wake_fr)))])
            disp(['wake: post ' num2str(nanmean(all_dls_post_wake_fr)) '+-' num2str(nanstd(all_dls_post_wake_fr)/sqrt(length(all_dls_post_wake_fr)))])
            disp(['wake t-test: t(' num2str(stats_wake.df) ') = ' num2str(stats_wake.tstat) ', P=' num2str(num2str(nanmean(p_wake)))]);
            
%% sFIG 4 | COMPARE SO AND DELTA

animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
animal_days = {[1:8],[1:6],[2:10],[4:17],[4:13],[1:5 7:9]};

delta_peak = [];
delta_valley = [];
delta_time = [];
mean_delta_peak = [];
mean_delta_valley = [];
mean_delta_time = [];

SO_peak = [];
SO_valley = [];
SO_time = [];
mean_SO_peak = [];
mean_SO_valley = [];
mean_SO_time = [];

for animal = 1:length(animal_id)
    for day = animal_days{animal}
        
        block_name = dir([data_path '\' animal_id{animal} '\sleep_activity\day_' num2str(day) '*.mat']);
        
        for block = 1:length(block_name)
        
            clearvars sleep_rhythms
            load([data_path '\' animal_id{animal} '\sleep_activity\' block_name(block).name])
            
            delta_peak = [delta_peak sleep_rhythms.so_delta.delta_peaks'];
            delta_valley = [delta_valley sleep_rhythms.so_delta.delta_troughs'];
            delta_time = [delta_time sleep_rhythms.so_delta.delta_dur'];
            mean_delta_peak = [mean_delta_peak mean(sleep_rhythms.so_delta.delta_peaks)];
            mean_delta_valley = [mean_delta_valley mean(sleep_rhythms.so_delta.delta_troughs)];
            mean_delta_time = [mean_delta_time median(sleep_rhythms.so_delta.delta_dur)];

            SO_peak = [SO_peak sleep_rhythms.so_delta.SO_peaks'];
            SO_valley = [SO_valley sleep_rhythms.so_delta.SO_troughs'];
            SO_time = [SO_time sleep_rhythms.so_delta.SO_dur'];
            mean_SO_peak = [mean_SO_peak mean(sleep_rhythms.so_delta.SO_peaks)];
            mean_SO_valley = [mean_SO_valley mean(sleep_rhythms.so_delta.SO_troughs)];
            mean_SO_time = [mean_SO_time median(sleep_rhythms.so_delta.SO_dur)];

        end
    end
end

figure;

    subplot(1,3,1); hold on;
        histogram(mean_delta_peak,[0:.2:4],'DisplayStyle','Stairs');
        histogram(mean_SO_peak,[0:.2:4],'DisplayStyle','Stairs');
        title(['peak ' num2str(length(mean_delta_peak)) ' sessions'])
        
    subplot(1,3,2); hold on;
        histogram(mean_delta_valley,[-3:.2:0],'DisplayStyle','Stairs');
        histogram(mean_SO_valley,[-3:.2:0],'DisplayStyle','Stairs');
        title('valley')

    subplot(1,3,3); hold on;
        histogram(mean_delta_time,[0:0.02:0.4],'DisplayStyle','Stairs');
        histogram(mean_SO_time,[0:0.02:0.4],'DisplayStyle','Stairs');
        title('duration')
               
%% sFIG 5 | NREM RHYTHM MODULATION

animal_id = {'LC3','LC4','LC6','LC5','LC2','LC1'};
animal_days = {[1:8],[1:6],[2:10],[4:17],[4:13],[1:5 7:9]};

spindle_m1_raster = [];
so_m1_raster = [];
delta_m1_raster = [];
control_m1_raster = [];

spindle_dls_raster = [];
so_dls_raster = [];
delta_dls_raster = [];
control_dls_raster = [];

for animal = 1:length(animal_id)

    animal_path = [data_path '\' animal_id{animal} '\cross_correlations'];

    for day = animal_days{animal}

        load([animal_path '\day_' num2str(day) '\CC_by_sleep_rhythm.mat']); 

        %%% FIND MINIMUM NUMBER OF RHYTHMS
        min_pre = min([size(all_m1_spindle(1).pre_m1_spindle_raster,1) size(all_m1_delta(1).pre_m1_delta_raster,1) size(all_m1_SO(1).pre_m1_SO_raster,1)]);
        min_post = min([size(all_m1_spindle(1).post_m1_spindle_raster,1) size(all_m1_delta(1).post_m1_delta_raster,1) size(all_m1_SO(1).post_m1_SO_raster,1)]);

        for unit = 1:length(all_m1_spindle)
            spindle_m1_raster = [spindle_m1_raster; mean([all_m1_spindle(unit).pre_m1_spindle_raster(randperm(size(all_m1_spindle(unit).pre_m1_spindle_raster,1),min_pre),:); all_m1_spindle(unit).post_m1_spindle_raster(randperm(size(all_m1_spindle(unit).post_m1_spindle_raster,1),min_post),:)])];
            delta_m1_raster = [delta_m1_raster; mean([all_m1_delta(unit).pre_m1_delta_raster(randperm(size(all_m1_delta(unit).pre_m1_delta_raster,1),min_pre),:); all_m1_delta(unit).post_m1_delta_raster(randperm(size(all_m1_delta(unit).post_m1_delta_raster,1),min_post),:)])];
            so_m1_raster = [so_m1_raster; mean([all_m1_SO(unit).pre_m1_SO_raster(randperm(size(all_m1_SO(unit).pre_m1_SO_raster,1),min_pre),:); all_m1_SO(unit).post_m1_SO_raster(randperm(size(all_m1_SO(unit).post_m1_SO_raster,1),min_post),:)])];
            control_m1_raster = [control_m1_raster; mean([all_m1_spindle(unit).pre_m1_spindle_control_raster(randperm(size(all_m1_spindle(unit).pre_m1_spindle_control_raster,1),min_pre),:); all_m1_spindle(unit).post_m1_spindle_control_raster(randperm(size(all_m1_spindle(unit).post_m1_spindle_control_raster,1),min_post),:)])];
        end

        for unit = 1:length(all_dls_spindle)
            spindle_dls_raster = [spindle_dls_raster; mean([all_dls_spindle(unit).pre_dls_spindle_raster(randperm(size(all_dls_spindle(unit).pre_dls_spindle_raster,1),min_pre),:); all_dls_spindle(unit).post_dls_spindle_raster(randperm(size(all_dls_spindle(unit).post_dls_spindle_raster,1),min_post),:)])];
            delta_dls_raster = [delta_dls_raster; mean([all_dls_delta(unit).pre_dls_delta_raster(randperm(size(all_dls_delta(unit).pre_dls_delta_raster,1),min_pre),:); all_dls_delta(unit).post_dls_delta_raster(randperm(size(all_dls_delta(unit).post_dls_delta_raster,1),min_post),:)])];
            so_dls_raster = [so_dls_raster; mean([all_dls_SO(unit).pre_dls_SO_raster(randperm(size(all_dls_SO(unit).pre_dls_SO_raster,1),min_pre),:); all_dls_SO(unit).post_dls_SO_raster(randperm(size(all_dls_SO(unit).post_dls_SO_raster,1),min_post),:)])];
            control_dls_raster = [control_dls_raster; mean([all_dls_spindle(unit).pre_dls_spindle_control_raster(randperm(size(all_dls_spindle(unit).pre_dls_spindle_control_raster,1),min_pre),:); all_dls_spindle(unit).post_dls_spindle_control_raster(randperm(size(all_dls_spindle(unit).post_dls_spindle_control_raster,1),min_post),:)])];
        end

    end
end

fig = figure;
hold on;
    cdfplot(100*(max(control_m1_raster(:,51:150)')-min(control_m1_raster(:,51:150)')));
    cdfplot(100*(max(spindle_m1_raster(:,51:150)')-min(spindle_m1_raster(:,51:150)')));
    cdfplot(100*(max(so_m1_raster(:,51:150)')-min(so_m1_raster(:,51:150)')));
    cdfplot(100*(max(delta_m1_raster(:,51:150)')-min(delta_m1_raster(:,51:150)')));
    xlim([0 15])

    [h p] = kstest2(100*(max(control_m1_raster(:,51:150)')-min(control_m1_raster(:,51:150)')),100*(max(spindle_m1_raster(:,51:150)')-min(spindle_m1_raster(:,51:150)')))
    [h p] = kstest2(100*(max(control_m1_raster(:,51:150)')-min(control_m1_raster(:,51:150)')),100*(max(so_m1_raster(:,51:150)')-min(so_m1_raster(:,51:150)')))
    [h p] = kstest2(100*(max(control_m1_raster(:,51:150)')-min(control_m1_raster(:,51:150)')),100*(max(delta_m1_raster(:,51:150)')-min(delta_m1_raster(:,51:150)')))

fig = figure;

    subplot(1,3,1)
        hold on;
        [xd, yd, delta, deltaCI, pval, cpval, sig] = shiftdhd_pbci(100*(max(spindle_m1_raster(:,51:150)')-min(spindle_m1_raster(:,51:150)')),100*(max(control_m1_raster(:,51:150)')-min(control_m1_raster(:,51:150)')),[.2 .4 .6 .8],2000,0.5,0)
        plot(delta,'color','b')    
        for n=1:length(delta)
            errorbar(n,delta(n),abs(delta(n)-deltaCI(n,1)),abs(delta(n)-deltaCI(n,2)),'color','b')    
        end
        xlim([0 5])
        plot([0 5],[0 0],'color','k')
        ylim([-.5 2.5])

    subplot(1,3,2)
        hold on;
        [xd, yd, delta, deltaCI, pval, cpval, sig] = shiftdhd_pbci(100*(max(so_m1_raster(:,51:150)')-min(so_m1_raster(:,51:150)')),100*(max(control_m1_raster(:,51:150)')-min(control_m1_raster(:,51:150)')),[.2 .4 .6 .8],2000,0.5,0)
        plot(delta,'color','r')    
        for n=1:length(delta)
            errorbar(n,delta(n),abs(delta(n)-deltaCI(n,1)),abs(delta(n)-deltaCI(n,2)),'color','r')    
        end
        xlim([0 5])
        plot([0 5],[0 0],'color','k')
        ylim([-.5 2.5])

    subplot(1,3,3)
        hold on;
        [xd, yd, delta, deltaCI, pval, cpval, sig] = shiftdhd_pbci(100*(max(delta_m1_raster(:,51:150)')-min(delta_m1_raster(:,51:150)')),100*(max(control_m1_raster(:,51:150)')-min(control_m1_raster(:,51:150)')),[.2 .4 .6 .8],2000,0.5,0)
        plot(delta,'color','g')    
        for n=1:length(delta)
            errorbar(n,delta(n),abs(delta(n)-deltaCI(n,1)),abs(delta(n)-deltaCI(n,2)),'color','g')    
        end
        xlim([0 5])
        plot([0 5],[0 0],'color','k')
        ylim([-.5 2.5])

fig = figure;
hold on;
    cdfplot(100*(max(control_dls_raster(:,51:150)')-min(control_dls_raster(:,51:150)')));
    cdfplot(100*(max(spindle_dls_raster(:,51:150)')-min(spindle_dls_raster(:,51:150)')));
    cdfplot(100*(max(so_dls_raster(:,51:150)')-min(so_dls_raster(:,51:150)')));
    cdfplot(100*(max(delta_dls_raster(:,51:150)')-min(delta_dls_raster(:,51:150)')));
    xlim([0 15])

    [h p] = kstest2(100*(max(control_dls_raster(:,51:150)')-min(control_dls_raster(:,51:150)')),100*(max(spindle_dls_raster(:,51:150)')-min(spindle_dls_raster(:,51:150)')))
    [h p] = kstest2(100*(max(control_dls_raster(:,51:150)')-min(control_dls_raster(:,51:150)')),100*(max(so_dls_raster(:,51:150)')-min(so_dls_raster(:,51:150)')))
    [h p] = kstest2(100*(max(control_dls_raster(:,51:150)')-min(control_dls_raster(:,51:150)')),100*(max(delta_dls_raster(:,51:150)')-min(delta_dls_raster(:,51:150)')))

fig = figure;

    subplot(1,3,1)
        hold on;
        [xd, yd, delta, deltaCI, pval, cpval, sig] = shiftdhd_pbci(100*(max(spindle_dls_raster(:,51:150)')-min(spindle_dls_raster(:,51:150)')),100*(max(control_dls_raster(:,51:150)')-min(control_dls_raster(:,51:150)')),[.2 .4 .6 .8],2000,0.5,0)
        plot(delta,'color','b')    
        for n=1:length(delta)
            errorbar(n,delta(n),abs(delta(n)-deltaCI(n,1)),abs(delta(n)-deltaCI(n,2)),'color','b')    
        end
        xlim([0 5])
        plot([0 5],[0 0],'color','k')
        ylim([-.5 2.5])

    subplot(1,3,2)
        hold on;
        [xd, yd, delta, deltaCI, pval, cpval, sig] = shiftdhd_pbci(100*(max(so_dls_raster(:,51:150)')-min(so_dls_raster(:,51:150)')),100*(max(control_dls_raster(:,51:150)')-min(control_dls_raster(:,51:150)')),[.2 .4 .6 .8],2000,0.5,0)
        plot(delta,'color','r')    
        for n=1:length(delta)
            errorbar(n,delta(n),abs(delta(n)-deltaCI(n,1)),abs(delta(n)-deltaCI(n,2)),'color','r')    
        end
        xlim([0 5])
        plot([0 5],[0 0],'color','k')
        ylim([-.5 2.5])

    subplot(1,3,3)
        hold on;
        [xd, yd, delta, deltaCI, pval, cpval, sig] = shiftdhd_pbci(100*(max(delta_dls_raster(:,51:150)')-min(delta_dls_raster(:,51:150)')),100*(max(control_dls_raster(:,51:150)')-min(control_dls_raster(:,51:150)')),[.2 .4 .6 .8],2000,0.5,0)
        plot(delta,'color','g')    
        for n=1:length(delta)
            errorbar(n,delta(n),abs(delta(n)-deltaCI(n,1)),abs(delta(n)-deltaCI(n,2)),'color','g')    
        end
        xlim([0 5])
        plot([0 5],[0 0],'color','k')
        ylim([-.5 2.5])
