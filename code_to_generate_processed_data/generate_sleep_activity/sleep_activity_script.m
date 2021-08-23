%% T102

    clc; clear all; close all;
    sleep_blocks = {[3 9 13],[14 19 21],[22 27],[30 33 37],[39 45],[48 49 54],[57 60],[62 65]};
    sleep_vids = {{'T102-(2015_2_11)-(10h50m)-SPONT1.avi','T102-(2015_2_11)-(13h57m)-SPONT1.avi','T102-(2015_2_11)-(16h46m)-SPONT1.avi'}, ...
                  {'T102-(2015_2_12)-(11h25m)-SPONT1.avi','T102-(2015_2_12)-(14h35m)-SPONT1.avi','T102-(2015_2_12)-(16h51m)-SPONT1.avi'}, ...
                    {'T102-(2015_2_13)-(10h5m)-SPONT1.avi','T102-(2015_2_13)-(13h32m)-SPONT1.avi'}, ...
                    {'T102-(2015_2_14)-(9h49m)-SPONT1.avi','T102-(2015_2_14)-(13h32m)-SPONT1.avi','T102-(2015_2_14)-(15h58m)-SPONT1.avi'}, ...
                    {'T102-(2015_2_15)-(10h7m)-SPONT1.avi','T102-(2015_2_15)-(13h25m)-SPONT1.avi'}, ...
                    {'T102-(2015_2_18)-(11h5m)-SPONT1.avi','T102-(2015_2_18)-(13h8m)-SPONT1.avi','T102-(2015_2_18)-(14h30m)-SPONT1.avi'}, ...
                    {'T102-(2015_2_19)-(9h53m)-SPONT1.avi','T102-(2015_2_19)-(12h42m)-SPONT1.avi'}, ...
                    {'T102-(2015_2_20)-(9h32m)-SPONT1.avi','T102-(2015_2_20)-(12h19m)-SPONT1.avi'}};

    for day = 1:length(sleep_blocks)
        for block = 1:length(sleep_blocks{day})

            wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\Block-' num2str(sleep_blocks{day}(block))],'STORE','Wave'); 
            video_start_stop = [find(wave.streams.Wave.data(1,:)>1.4,1,'first')/wave.streams.Wave.fs find(wave.streams.Wave.data(1,:)>1.4,1,'last')/wave.streams.Wave.fs];
            lfp = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\Block-' num2str(sleep_blocks{day}(block))],'Type',4);    
            save_path = 'C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T102\sleep_activity\';
            [beh_state, sleep_rhythms] = detect_sleep(lfp.streams.LFPs.data(1:16,:),lfp.streams.LFPs.data(17:32,:),lfp.streams.LFPs.fs, ['\\MyCloudPR4100\Public\SL_Backups\T102\sleep_data\day_' num2str(day) '\' sleep_vids{day}{block}], video_start_stop,.25,day,block,[],[],save_path,[]);

        end
    end

%% T107

    clc; clear all; close all;
    sleep_blocks = {[1 2],[3 6],[9 14],[18 24],[26 36],[41 47]};
    sleep_vids = {{'T107-(2015_4_1)-(10h19m)-SPONT1.avi','T107-(2015_4_1)-(13h18m)-SPONT1.avi'}, ...
                  {'T107-(2015_4_2)-(10h43m)-SPONT1.avi','T107-(2015_4_2)-(13h15m)-SPONT1.avi'}, ...
                    {'T107-(2015_4_3)-(11h3m)-SPONT1.avi','T107-(2015_4_3)-(15h30m)-SPONT1.avi'}, ...
                    {'T107-(2015_4_4)-(10h3m)-SPONT1.avi','T107-(2015_4_4)-(13h20m)-SPONT1.avi'}, ...
                    {'T107-(2015_4_5)-(12h16m)-SPONT1.avi','T107-(2015_4_5)-(15h3m)-SPONT1.avi'}, ...
                    {'T107-(2015_4_6)-(11h27m)-SPONT1.avi','T107-(2015_4_6)-(14h7m)-SPONT1.avi'}};

    for day = 1:length(sleep_blocks)
        for block = 1:length(sleep_blocks{day})

            wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T107\neural_data\Block-' num2str(sleep_blocks{day}(block))],'STORE','Wave'); 
            video_start_stop = [find(wave.streams.Wave.data(1,:)>1.4,1,'first')/wave.streams.Wave.fs find(wave.streams.Wave.data(1,:)>1.4,1,'last')/wave.streams.Wave.fs];
            lfp = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T107\neural_data\Block-' num2str(sleep_blocks{day}(block))],'Type',4);    
            save_path = 'C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T107\sleep_activity\';
            [beh_state, sleep_rhythms] = detect_sleep(lfp.streams.LFPs.data(1:16,:),lfp.streams.LFPs.data(17:32,:),lfp.streams.LFPs.fs, ['\\MyCloudPR4100\Public\SL_Backups\T107\sleep_data\day_' num2str(day) '\' sleep_vids{day}{block}], video_start_stop,.5,day,block,[],[],save_path,[]);

        end
    end

%% T200

    clc; clear all; close all;
    sleep_blocks = {[1],[2 8],[9 11 14],[15 16 19],[20 23 24],[25 27],[28 30],[32 36],[37 41],[42 46 47]};
    sleep_vids = {{'T200-(2016_7_9)-(16h6m)-SPONT.avi'}, ...
                  {'T200-(2016_7_11)-(9h41m)-SPONT.avi','T200-(2016_7_11)-(16h22m)-SPONT.avi'}, ...
                    {'T200-(2016_7_12)-(9h20m)-SPONT.avi','T200-(2016_7_12)-(11h29m)-SPONT.avi','T200-(2016_7_12)-(14h30m)-SPONT.avi'}, ...
                    {'T200-(2016_7_13)-(8h36m)-SPONT.avi','T200-(2016_7_13)-(10h57m)-SPONT.avi','T200-(2016_7_13)-(12h58m)-SPONT.avi'}, ...
                    {'T200-(2016_7_14)-(10h59m)-SPONT.avi','T200-(2016_7_14)-(14h53m)-SPONT.avi','T200-(2016_7_14)-(17h35m)-SPONT.avi'}, ...
                    {[],'T200-(2016_7_15)-(11h57m)-SPONT.avi'}, ...
                    {'T200-(2016_7_16)-(11h49m)-SPONT.avi','T200-(2016_7_16)-(13h56m)-SPONT.avi'}, ...
                    {'T200-(2016_7_18)-(8h59m)-SPONT.avi','T200-(2016_7_18)-(13h30m)-SPONT.avi'}, ...
                    {'T200-(2016_7_19)-(8h55m)-SPONT.avi','T200-(2016_7_19)-(13h42m)-SPONT.avi'}, ...
                    {'T200-(2016_7_20)-(8h49m)-SPONT.avi','T200-(2016_7_20)-(13h36m)-SPONT.avi','T200-(2016_7_20)-(14h46m)-SPONT.avi'}};

    for day = 1:length(sleep_blocks)
        for block = 1:length(sleep_blocks{day})

            wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\Block-' num2str(sleep_blocks{day}(block))],'STORE','Wave'); 
            video_start_stop = [find(wave.streams.Wave.data(1,:)>1.4,1,'first')/wave.streams.Wave.fs find(wave.streams.Wave.data(1,:)>1.4,1,'last')/wave.streams.Wave.fs];
            lfp = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\Block-' num2str(sleep_blocks{day}(block))],'Type',4);    
            save_path = 'C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T200\sleep_activity\';
            [beh_state, sleep_rhythms] = detect_sleep(lfp.streams.LFPs.data(17:32,:),lfp.streams.LFPs.data(1:16,:),lfp.streams.LFPs.fs, ['\\MyCloudPR4100\Public\SL_Backups\T200\sleep_data\day_' num2str(day-1) '\' sleep_vids{day}{block}], video_start_stop,.5,day,block,[],[],save_path,[]);

        end
    end

%% T201

    clc; clear all; close all;
    sleep_blocks = {[1],[2 3],[4 5 7],[8 10],[11 13],[15 17],[18 20],[21 26],[27 29],[30 32],[33 34 37],[38 40],[41 44],[45 48],[49 52],[53 56],[57 61]};
    sleep_vids = {{'T201-(2016_7_28)-(16h32m)-SPONT.avi'}, ...
                    {'T201-(2016_7_29)-(11h5m)-SPONT.avi','T201-(2016_7_29)-(13h31m)-SPONT.avi'}, ...
                    {'T201-(2016_8_1)-(8h42m)-SPONT.avi','T201-(2016_8_1)-(10h4m)-SPONT.avi','T201-(2016_8_1)-(14h51m)-SPONT.avi'}, ...
                    {'T201-(2016_8_2)-(9h52m)-SPONT.avi','T201-(2016_8_2)-(13h52m)-SPONT.avi'}, ...
                    {'T201-(2016_8_3)-(10h1m)-SPONT.avi','T201-(2016_8_3)-(15h1m)-SPONT.avi'}, ...
                    {'T201-(2016_8_4)-(11h57m)-SPONT.avi','T201-(2016_8_4)-(14h37m)-SPONT.avi'}, ...
                    {'T201-(2016_8_5)-(13h19m)-SPONT.avi','T201-(2016_8_5)-(16h22m)-SPONT.avi'}, ...
                    {'T201-(2016_8_6)-(11h7m)-SPONT.avi','T201-(2016_8_6)-(15h20m)-SPONT.avi'}, ...
                    {'T201-(2016_8_7)-(11h32m)-SPONT.avi','T201-(2016_8_7)-(14h42m)-SPONT.avi'}, ...
                    {'T201-(2016_8_8)-(9h27m)-SPONT.avi','T201-(2016_8_8)-(14h0m)-SPONT.avi'}, ...
                    {'T201-(2016_8_9)-(10h23m)-SPONT.avi','T201-(2016_8_9)-(12h56m)-SPONT.avi','T201-(2016_8_9)-(16h7m)-SPONT.avi'}, ...
                    {[],'T201-(2016_8_10)-(13h2m)-SPONT.avi'}, ...
                    {[],'T201-(2016_8_11)-(15h40m)-SPONT.avi'}, ...
                    {'T201-(2016_8_12)-(10h24m)-SPONT.avi','T201-(2016_8_12)-(15h9m)-SPONT.avi'}, ...
                    {'T201-(2016_8_13)-(11h49m)-SPONT.avi','T201-(2016_8_13)-(15h54m)-SPONT.avi'}, ...
                    {'T201-(2016_8_14)-(11h38m)-SPONT.avi','T201-(2016_8_14)-(15h54m)-SPONT.avi'}, ...
                    {'T201-(2016_8_15)-(9h20m)-SPONT.avi','T201-(2016_8_15)-(13h58m)-SPONT.avi'}};

    for day = 1:length(sleep_blocks)
        for block = 1:length(sleep_blocks{day})

            if day == 13 && block == 1
                lfp = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(sleep_blocks{day}(block))],'Type',4,'T2',7200);
                video_start_stop = [];
            else
                wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(sleep_blocks{day}(block))],'STORE','Wave'); 
                video_start_stop = [find(wave.streams.Wave.data(1,:)>1.4,1,'first')/wave.streams.Wave.fs find(wave.streams.Wave.data(1,:)>1.4,1,'last')/wave.streams.Wave.fs];
                lfp = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(sleep_blocks{day}(block))],'Type',4);    
            end
            save_path = 'C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T201\sleep_activity\';
            [beh_state, sleep_rhythms] = detect_sleep(lfp.streams.LFPs.data(1:2:63,:),lfp.streams.LFPs.data(34:2:64,:),lfp.streams.LFPs.fs, ['\\MyCloudPR4100\Public\SL_Backups\T201\sleep_data\day_' num2str(day-2) '\' sleep_vids{day}{block}], video_start_stop,.5,day,block,[],[],save_path,[]);
       
        end
    end

%% SLDD2

clc; clear all; close all;

sleep_blocks = {{'SLDD2-191113-143329'}, ...
                {'SLDD2-191114-121834'}, ...
                {'SLDD2-191115-154323'}, ...
                {'SLDD2-191116-145119','SLDD2-191116-193647'}, ...
                {'SLDD2-191117-150115','SLDD2-191117-181931'}, ...
                {'SLDD2-191118-154100','SLDD2-191118-205741'}, ...
                {'SLDD2-191119-160043','SLDD2-191119-193444'}, ...
                {'SLDD2-191120-140404','SLDD2-191120-170507'}, ...
                {'SLDD2-191121-170233','SLDD2-191121-192200'}, ...
                {'SLDD2-191122-152109','SLDD2-191122-184340','SLDD2-191122-204428'}, ...
                {'SLDD2-191123-173250','SLDD2-191123-205347'}, ...
                {'SLDD2-191124-152841','SLDD2-191124-183235'}, ...
                {'SLDD2-191125-141910','SLDD2-191125-171539'}, ...
                {'SLDD2-191126-173517','SLDD2-191126-195858'}};

sleep_vids = {{'SLDD2-(2019_11_13)-(14h34m)-SPONT.avi'}, ...
              {'SLDD2-(2019_11_14)-(12h18m)-SPONT.avi'}, ...
              {'SLDD2-(2019_11_15)-(15h43m)-SPONT.avi'}, ...
              {'SLDD2-(2019_11_16)-(14h56m)-SPONT.avi','SLDD2-(2019_11_16)-(19h37m)-SPONT.avi'}, ...
              {'SLDD2-(2019_11_17)-(15h2m)-SPONT.avi','SLDD2-(2019_11_17)-(18h19m)-SPONT.avi'}, ...
              {'SLDD2-(2019_11_18)-(15h41m)-SPONT.avi','SLDD2-(2019_11_18)-(20h57m)-SPONT.avi'}, ...
              {'SLDD2-(2019_11_19)-(16h0m)-SPONT.avi','SLDD2-(2019_11_19)-(19h35m)-SPONT.avi'}, ...
              {'SLDD2-(2019_11_20)-(14h4m)-SPONT.avi','SLDD2-(2019_11_20)-(17h7m)-SPONT.avi'}, ...
              {'SLDD2-(2019_11_21)-(17h3m)-SPONT.avi','SLDD2-(2019_11_21)-(19h22m)-SPONT.avi'}, ...
              {'SLDD2-(2019_11_22)-(15h21m)-SPONT.avi','SLDD2-(2019_11_22)-(18h44m)-SPONT.avi','SLDD2-(2019_11_22)-(20h44m)-SPONT.avi'}, ...
              {'SLDD2-(2019_11_23)-(17h33m)-SPONT.avi','SLDD2-(2019_11_23)-(20h54m)-SPONT.avi'}, ...
              {'SLDD2-(2019_11_24)-(15h29m)-SPONT.avi','SLDD2-(2019_11_24)-(18h33m)-SPONT.avi'}, ...
              {'SLDD2-(2019_11_25)-(14h19m)-SPONT.avi','SLDD2-(2019_11_25)-(17h16m)-SPONT.avi'}, ...
              {'SLDD2-(2019_11_26)-(17h35m)-SPONT.avi','SLDD2-(2019_11_26)-(19h59m)-SPONT.avi'}};
    
for day = 1:length(sleep_blocks)
    for block = 1:length(sleep_blocks{day})
        
        m1_lfp = [];
        for chan = [8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33]
            data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\RS4_data\' sleep_blocks{day}{block}],'CHANNEL',chan);
            m1_lfp = [m1_lfp; downsample(data.RSn1.data,20)];
        end
        
        dls_lfp = [];
        for chan = [2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63]
            data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\RS4_data\' sleep_blocks{day}{block}],'CHANNEL',chan);
            dls_lfp = [dls_lfp; downsample(data.RSn1.data,20)];
        end
        
        lfp_samp_rate = data.RSn1.fs/20;

        wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\TDT_tanks\' sleep_blocks{day}{block}],'TYPE',[2]);
        if isempty(wave.epocs)
            video_start_stop = [];
        else
            video_start_stop = [min(wave.epocs.PtC0.onset) max(wave.epocs.PtC0.onset)];
        end        
        
        save_path = 'C:\Users\Stefan\Desktop\Lemke_etal_2021\data\SLDD2\sleep_activity\';
        if day==12 && block==2
            [beh_state, sleep_rhythms] = detect_sleep(double(m1_lfp),double(dls_lfp),lfp_samp_rate, [], video_start_stop,.5,day,block,[],[],save_path,[]);
        else
            [beh_state, sleep_rhythms] = detect_sleep(double(m1_lfp),double(dls_lfp),lfp_samp_rate, ['\\MyCloudPR4100\Public\SL_Backups\SLDD2\sleep_data\day_' num2str(day-3) '\' sleep_vids{day}{block}], video_start_stop,.5,day,block,[],[],save_path,[]);
        end
    end
end

%% T398

clc; clear all; close all;

    sleep_blocks = {{'T398-190226-112203','T398-190226-160611'}, ...
                      {'T398-190227-113133','T398-190227-154618'}, ...
                      {'T398-190228-105336','T398-190228-141352'}, ...
                      {'T398-190301-112846','T398-190301-151942'}, ...
                      {'T398-190302-131102','T398-190302-162611'}, ...
                      {'T398-190303-140804','T398-190303-171207'}, ...
                      {'T398-190304-114807','T398-190304-145843'}, ...
                      {'T398-190305-130625','T398-190305-161031'}, ...
                      {'T398-190306-091252','T398-190306-124648'}};

    sleep_vids = {{'T398-(2019_2_26)-(11h22m)-SPONT.avi','T398-(2019_2_26)-(16h7m)-SPONT.avi'}, ...
                  {'T398-(2019_2_27)-(11h32m)-SPONT.avi','T398-(2019_2_27)-(15h46m)-SPONT.avi'}, ...
                  {'T398-(2019_2_28)-(10h54m)-SPONT.avi','T398-(2019_2_28)-(14h14m)-SPONT.avi'}, ...
                  {'T398-(2019_3_1)-(11h29m)-SPONT.avi','T398-(2019_3_1)-(15h28m)-SPONT.avi'}, ...
                  {'T398-(2019_3_2)-(13h11m)-SPONT.avi','T398-(2019_3_2)-(16h26m)-SPONT.avi'}, ...
                  {'T398-(2019_3_3)-(14h8m)-SPONT','T398-(2019_3_3)-(17h12m)-SPONT'}, ...
                  {'T398-(2019_3_4)-(11h48m)-SPONT.avi','T398-(2019_3_4)-(14h59m)-SPONT.avi'}, ...
                  {'T398-(2019_3_5)-(13h7m)-SPONT.avi','T398-(2019_3_5)-(16h10m)-SPONT.avi'}, ...
                  {'T398-(2019_3_6)-(9h13m)-SPONT.avi','T398-(2019_3_6)-(12h47m)-SPONT.avi'}};

    % BANKS
    bank_1_2 = {[5 6 7 8 9 21 39 40 41 42 43 46 51 52 53 54 56 57 58 59 60 61 100 102 105 107 108 110 111 113 114 115 116 119 120 122 123], ...
        [129 130 133 140 141 153 159 160 166 168 170 173 176 179 180 181 187 188 194 195 196 197 198 199 200 201 202 203 205 0206 207 211 212 214 223 224 231]};
    bank_2_1 = {[1 2 3 5 13 19 25 31 32 38 40 42 44 45 48 51 52 53 60 67 68 69 70 71 72 73 74 75 77 78 79 83 84 86 95 96 103], ...
        [135 136 137 139 140 149 167 168 169 170 171 179 180 181 182 184 185 186 187 188 189 228 233 235 236 238 239 241 242 243 244 247 248 249 250 251 255]};

    % BANK ORDER (M1 then DLS)
    bank_order = {[1 2], [1 2], [1 2], [2 1], [1 2], [2 1], [1 2], [1 2], [2 1]};
    
for day = 1:length(sleep_blocks)
    for block = 1:length(sleep_blocks{day})
        
        if bank_order{day}(1)==1
            M1_bank = bank_1_2{1};
            M1_offset = 1;
            DLS_bank = bank_1_2{2};
            DLS_offset = 129;
        elseif bank_order{day}(1)==2
            M1_bank = bank_2_1{2};
            M1_offset = 129;
            DLS_bank = bank_2_1{1};
            DLS_offset = 1;
        end
    
        M1_mapping = {(M1_offset+[15 16 14 17 13 18 12 19]), ...         %01
                     (M1_offset+[11 20 10 21 9 22 8 23]), ...            %02
                     [], ...                                             %03
                     [], ...                                             %04
                     (M1_offset+[63 32 62 33 61 34 60 35]), ...          %05
                     (M1_offset+[59 36 58 37 57 38 56 39]), ...          %06
                     [], ...                                             %07
                     [], ...                                             %08
                     (M1_offset+[79 80 78 81 77 82 76 83]), ...          %09
                     (M1_offset+[75 84 74 85 73 86 72 87]), ...          %10
                     [], ...                                             %11
                     [], ...                                             %12
                     (M1_offset+[127 96 126 97 125 98 124 99]), ...      %13
                     (M1_offset+[123 100 122 101 121 102 120 103]), ...  %14
                     [], ...                                             %15
                     [], ...                                             %16
                     (M1_offset+[111 112 110 113 109 114 108 115]), ...  %17
                     (M1_offset+[107 116 106 117 105 118 104 119]), ...  %18
                     [], ...                                             %19
                     [], ...                                             %20
                     (M1_offset+[95 64 94 65 93 66 92 67]), ...          %21
                     (M1_offset+[91 68 90 69 89 70 88 71]), ...          %22
                     [], ...                                             %23
                     [], ...                                             %24
                     (M1_offset+[47 48 46 49 45 50 44 51]), ...          %25
                     (M1_offset+[43 52 42 53 41 54 40 55]), ...          %26
                     [], ...                                             %27
                     [], ...                                             %28
                     (M1_offset+[31 0 30 1 29 2 28 3]), ...              %29
                     (M1_offset+[27 4 26 5 25 6 24 7]), ...              %30
                     [], ...                                             %31
                     []};                                                %32

        DLS_mapping = {[], ...                                           %01
                     [], ...                                             %02
                     [], ...                                             %03
                     [], ...                                             %04
                     [], ...                                             %05
                     [], ...                                             %06
                     [], ...                                             %07
                     [], ...                                             %08
                     [], ...                                             %09
                     [], ...                                             %10
                     [], ...                                             %11
                     [], ...                                             %12
                     [], ...                                             %13
                     [], ...                                             %14
                     [], ...                                             %15
                     [], ...                                             %16
                     (DLS_offset+[15 16 14 17 13 18 12 19]), ...         %17
                     (DLS_offset+[11 20 10 21 9 22 8 23]), ...           %18
                     (DLS_offset+[7 24 6 25 5 26 4 27]), ...             %19
                     (DLS_offset+[3 28 2 29 1 30 0 31]), ...             %20
                     (DLS_offset+[63 32 62 33 61 34 60 35]), ...         %21
                     (DLS_offset+[59 36 58 37 57 38 56 39]), ...         %22
                     (DLS_offset+[55 40 54 41 53 42 52 43]), ...         %2 
                     (DLS_offset+[51 44 50 45 49 46 48 47]), ...         %24
                     (DLS_offset+[79 80 78 81 77 82 76 83]), ...         %25
                     (DLS_offset+[75 84 74 85 73 86 72 87]), ...         %26
                     (DLS_offset+[71 88 70 89 69 90 68 91]), ...         %27
                     (DLS_offset+[67 92 66 93 65 94 64 95]), ...         %28
                     (DLS_offset+[127 96 126 97 125 96 124 95]), ...     %29
                     (DLS_offset+[123 100 122 101 121 102 120 103]), ... %30
                     (DLS_offset+[119 104 118 105 117 106 116 107]), ... %31
                     (DLS_offset+[115 108 114 109 113 110 112 111])};    %32

        if day==7 && block==2
            video_start_stop = [10 7210];
        else
            wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\TDT_tanks\M1-DLS_256Probe-190222-144047\' sleep_blocks{day}{block}]);
            video_start_stop = [min(wave.epocs.PC1_.onset) max(wave.epocs.PC1_.onset)];
            if day == 2 && block == 1
                video_start_stop(2) = video_start_stop(2)+7200;
            end
        end

        m1_lfp = [];
        for shank = [1 2 5 6 9 10 13 14 17 18 21 22 25 26 29 30]
         disp(shank);
         [tmp_chans, ~] = intersect(M1_mapping{shank},M1_bank);
         for electrode = 1:length(tmp_chans)
             data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\RS4_data\M1-DLS_256Probe-190222-144047\' sleep_blocks{day}{block}],'CHANNEL',tmp_chans(electrode));
             m1_lfp = [m1_lfp; downsample(data.RSn1.data,20)];
         end
        end
        dls_lfp = [];
        for shank = 17:32
         disp(shank);
         [tmp_chans, ~] = intersect(DLS_mapping{shank},DLS_bank);
         for electrode = 1:length(tmp_chans)
             data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\RS4_data\M1-DLS_256Probe-190222-144047\' sleep_blocks{day}{block}],'CHANNEL',tmp_chans(electrode));
             dls_lfp = [dls_lfp; downsample(data.RSn1.data,20)];
         end
        end
        lfp_samp_rate = data.RSn1.fs/20;
        
        save_path = 'C:\Users\Stefan\Desktop\Lemke_etal_2021\data\T398\sleep_activity\';
        
        [beh_state, sleep_rhythms] = detect_sleep(double(m1_lfp),double(dls_lfp),lfp_samp_rate, ['\\MyCloudPR4100\Public\SL_Backups\T398\sleep_data\day_' num2str(day) '\' sleep_vids{day}{block}], video_start_stop,0.01, day, block,[],[],save_path,1);
            
    end
end
