function [m1_reach, dls_reach] = reach_modulation(m1_spiking,dls_spiking,beh_markers)
%
% function to compute reach modulation
%
    
%% M1

    unit_count = 1;
    for chan = 1:size(m1_spiking,2)
        for unit = 1:size(m1_spiking,3)
            
            if isempty(m1_spiking{1,chan,unit})
                continue
            end
                
            raster = [];
            raster_10ms = [];
            raster_40ms = [];
            raster_50ms = [];
            raster_100ms = [];
            raster_200ms = [];
                
            baseline = [];
            baseline_10ms = [];
            baseline_40ms = [];
            baseline_50ms = [];
            baseline_100ms = [];
            baseline_200ms = [];
            
            for block = 1:length(beh_markers)
                for trial = 1:length(beh_markers{block})
                    if isnan(beh_markers{block}(trial))
                        continue
                    end
                    
                    raster = [raster; histcounts(m1_spiking{block,chan,unit}(m1_spiking{block,chan,unit}>beh_markers{block}(trial)-5 & m1_spiking{block,chan,unit}<beh_markers{block}(trial)+5)-beh_markers{block}(trial),[-5:0.025:5])];
                    raster_10ms = [raster_10ms; histcounts(m1_spiking{block,chan,unit}(m1_spiking{block,chan,unit}>beh_markers{block}(trial)-5 & m1_spiking{block,chan,unit}<beh_markers{block}(trial)+5)-beh_markers{block}(trial),[-5:0.010:5])];
                    raster_40ms = [raster_40ms; histcounts(m1_spiking{block,chan,unit}(m1_spiking{block,chan,unit}>beh_markers{block}(trial)-5 & m1_spiking{block,chan,unit}<beh_markers{block}(trial)+5)-beh_markers{block}(trial),[-5:0.040:5])];
                    raster_50ms = [raster_50ms; histcounts(m1_spiking{block,chan,unit}(m1_spiking{block,chan,unit}>beh_markers{block}(trial)-5 & m1_spiking{block,chan,unit}<beh_markers{block}(trial)+5)-beh_markers{block}(trial),[-5:0.050:5])];
                    raster_100ms = [raster_100ms; histcounts(m1_spiking{block,chan,unit}(m1_spiking{block,chan,unit}>beh_markers{block}(trial)-5 & m1_spiking{block,chan,unit}<beh_markers{block}(trial)+5)-beh_markers{block}(trial),[-5:0.100:5])];
                    raster_200ms = [raster_200ms; histcounts(m1_spiking{block,chan,unit}(m1_spiking{block,chan,unit}>beh_markers{block}(trial)-5 & m1_spiking{block,chan,unit}<beh_markers{block}(trial)+5)-beh_markers{block}(trial),[-5:0.200:5])];
                    
                    baseline = [baseline; histcounts(m1_spiking{block,chan,unit}(m1_spiking{block,chan,unit}>beh_markers{block}(trial)-15 & m1_spiking{block,chan,unit}<beh_markers{block}(trial)-5)-beh_markers{block}(trial)-10,[-5:0.025:5])];
                    baseline_10ms = [baseline_10ms; histcounts(m1_spiking{block,chan,unit}(m1_spiking{block,chan,unit}>beh_markers{block}(trial)-15 & m1_spiking{block,chan,unit}<beh_markers{block}(trial)-5)-beh_markers{block}(trial)-10,[-5:0.010:5])];
                    baseline_40ms = [baseline_40ms; histcounts(m1_spiking{block,chan,unit}(m1_spiking{block,chan,unit}>beh_markers{block}(trial)-15 & m1_spiking{block,chan,unit}<beh_markers{block}(trial)-5)-beh_markers{block}(trial)-10,[-5:0.040:5])];
                    baseline_50ms = [baseline_50ms; histcounts(m1_spiking{block,chan,unit}(m1_spiking{block,chan,unit}>beh_markers{block}(trial)-15 & m1_spiking{block,chan,unit}<beh_markers{block}(trial)-5)-beh_markers{block}(trial)-10,[-5:0.050:5])];
                    baseline_100ms = [baseline_100ms; histcounts(m1_spiking{block,chan,unit}(m1_spiking{block,chan,unit}>beh_markers{block}(trial)-15 & m1_spiking{block,chan,unit}<beh_markers{block}(trial)-5)-beh_markers{block}(trial)-10,[-5:0.100:5])];
                    baseline_200ms = [baseline_200ms; histcounts(m1_spiking{block,chan,unit}(m1_spiking{block,chan,unit}>beh_markers{block}(trial)-15 & m1_spiking{block,chan,unit}<beh_markers{block}(trial)-5)-beh_markers{block}(trial)-10,[-5:0.200:5])];
                
                end
            end
            
            m1_reach(unit_count).m1_chan_unit = [chan unit];
            
            m1_reach(unit_count).raster = raster;
            m1_reach(unit_count).raster_10ms = raster_10ms;
            m1_reach(unit_count).raster_40ms = raster_40ms;
            m1_reach(unit_count).raster_50ms = raster_50ms;
            m1_reach(unit_count).raster_100ms = raster_100ms;
            m1_reach(unit_count).raster_200ms = raster_200ms;
            
            m1_reach(unit_count).baseline = baseline;
            m1_reach(unit_count).baseline_10ms = baseline_10ms;
            m1_reach(unit_count).baseline_40ms = baseline_40ms;
            m1_reach(unit_count).baseline_50ms = baseline_50ms;
            m1_reach(unit_count).baseline_100ms = baseline_100ms;
            m1_reach(unit_count).baseline_200ms = baseline_200ms;
            
            unit_count = unit_count + 1;

        end
    end
    
%% DLS

    unit_count = 1;
    for chan = 1:size(dls_spiking,2)
        for unit = 1:size(dls_spiking,3)
            
            if isempty(dls_spiking{1,chan,unit})
                continue
            end
                
            raster = [];
            raster_10ms = [];
            raster_40ms = [];
            raster_50ms = [];
            raster_100ms = [];
            raster_200ms = [];
                
            baseline = [];
            baseline_10ms = [];
            baseline_40ms = [];
            baseline_50ms = [];
            baseline_100ms = [];
            baseline_200ms = [];
            
            for block = 1:length(beh_markers)
                for trial = 1:length(beh_markers{block})
                    if isnan(beh_markers{block}(trial))
                        continue
                    end
                    
                    raster = [raster; histcounts(dls_spiking{block,chan,unit}(dls_spiking{block,chan,unit}>beh_markers{block}(trial)-5 & dls_spiking{block,chan,unit}<beh_markers{block}(trial)+5)-beh_markers{block}(trial),[-5:0.025:5])];
                    raster_10ms = [raster_10ms; histcounts(dls_spiking{block,chan,unit}(dls_spiking{block,chan,unit}>beh_markers{block}(trial)-5 & dls_spiking{block,chan,unit}<beh_markers{block}(trial)+5)-beh_markers{block}(trial),[-5:0.010:5])];
                    raster_40ms = [raster_40ms; histcounts(dls_spiking{block,chan,unit}(dls_spiking{block,chan,unit}>beh_markers{block}(trial)-5 & dls_spiking{block,chan,unit}<beh_markers{block}(trial)+5)-beh_markers{block}(trial),[-5:0.040:5])];
                    raster_50ms = [raster_50ms; histcounts(dls_spiking{block,chan,unit}(dls_spiking{block,chan,unit}>beh_markers{block}(trial)-5 & dls_spiking{block,chan,unit}<beh_markers{block}(trial)+5)-beh_markers{block}(trial),[-5:0.050:5])];
                    raster_100ms = [raster_100ms; histcounts(dls_spiking{block,chan,unit}(dls_spiking{block,chan,unit}>beh_markers{block}(trial)-5 & dls_spiking{block,chan,unit}<beh_markers{block}(trial)+5)-beh_markers{block}(trial),[-5:0.100:5])];
                    raster_200ms = [raster_200ms; histcounts(dls_spiking{block,chan,unit}(dls_spiking{block,chan,unit}>beh_markers{block}(trial)-5 & dls_spiking{block,chan,unit}<beh_markers{block}(trial)+5)-beh_markers{block}(trial),[-5:0.200:5])];
                    
                    baseline = [baseline; histcounts(dls_spiking{block,chan,unit}(dls_spiking{block,chan,unit}>beh_markers{block}(trial)-15 & dls_spiking{block,chan,unit}<beh_markers{block}(trial)-5)-beh_markers{block}(trial)-10,[-5:0.025:5])];
                    baseline_10ms = [baseline_10ms; histcounts(dls_spiking{block,chan,unit}(dls_spiking{block,chan,unit}>beh_markers{block}(trial)-15 & dls_spiking{block,chan,unit}<beh_markers{block}(trial)-5)-beh_markers{block}(trial)-10,[-5:0.010:5])];
                    baseline_40ms = [baseline_40ms; histcounts(dls_spiking{block,chan,unit}(dls_spiking{block,chan,unit}>beh_markers{block}(trial)-15 & dls_spiking{block,chan,unit}<beh_markers{block}(trial)-5)-beh_markers{block}(trial)-10,[-5:0.040:5])];
                    baseline_50ms = [baseline_50ms; histcounts(dls_spiking{block,chan,unit}(dls_spiking{block,chan,unit}>beh_markers{block}(trial)-15 & dls_spiking{block,chan,unit}<beh_markers{block}(trial)-5)-beh_markers{block}(trial)-10,[-5:0.050:5])];
                    baseline_100ms = [baseline_100ms; histcounts(dls_spiking{block,chan,unit}(dls_spiking{block,chan,unit}>beh_markers{block}(trial)-15 & dls_spiking{block,chan,unit}<beh_markers{block}(trial)-5)-beh_markers{block}(trial)-10,[-5:0.100:5])];
                    baseline_200ms = [baseline_200ms; histcounts(dls_spiking{block,chan,unit}(dls_spiking{block,chan,unit}>beh_markers{block}(trial)-15 & dls_spiking{block,chan,unit}<beh_markers{block}(trial)-5)-beh_markers{block}(trial)-10,[-5:0.200:5])];
                
                end
            end
            
            dls_reach(unit_count).dls_chan_unit = [chan unit];
            
            dls_reach(unit_count).raster = raster;
            dls_reach(unit_count).raster_10ms = raster_10ms;
            dls_reach(unit_count).raster_40ms = raster_40ms;
            dls_reach(unit_count).raster_50ms = raster_50ms;
            dls_reach(unit_count).raster_100ms = raster_100ms;
            dls_reach(unit_count).raster_200ms = raster_200ms;
            
            dls_reach(unit_count).baseline = baseline;
            dls_reach(unit_count).baseline_10ms = baseline_10ms;
            dls_reach(unit_count).baseline_40ms = baseline_40ms;
            dls_reach(unit_count).baseline_50ms = baseline_50ms;
            dls_reach(unit_count).baseline_100ms = baseline_100ms;
            dls_reach(unit_count).baseline_200ms = baseline_200ms;
            
            unit_count = unit_count + 1;

        end
    end
    
%% PLOT
% 
%     tmp_m1 = [];
%     for unit = 1:length(m1_reach)
%         tmp_m1 = [tmp_m1; mean(m1_reach(unit).raster)];
%         
%         figure;
%         subplot(2,1,1)
%         imagesc(m1_reach(unit).raster)
%         xlim([1 400])
%         subplot(2,1,2);
%         plot(mean(m1_reach(unit).raster))
%         xlim([1 400])
% 
%     end
% 
%     tmp_dls = [];
%     for unit = 1:length(dls_reach)
%         tmp_dls = [tmp_dls; mean(dls_reach(unit).raster)];
%         
%         figure;
%         subplot(2,1,1)
%         imagesc(dls_reach(unit).raster)
%         xlim([1 400])
%         subplot(2,1,2);
%         plot(mean(dls_reach(unit).raster))
%         xlim([1 400])
%         
%     end
% 
%     figure;
%     subplot(1,2,1)
%         hold on;
%         plot(mean(tmp_m1),'r');
%         plot(mean(tmp_dls),'b');
%     subplot(1,2,2)
%         hold on;
%         plot(zscore(mean(tmp_m1)),'r');
%         plot(zscore(mean(tmp_dls)),'b');
%         
end
