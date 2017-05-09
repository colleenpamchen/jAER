function Parse_DVS128_Recording( filename, streams, save_information )
% Parses recorded binary eDVS128 data into files.
%
% This funtion reads binary eDVS128 files recorded with the java recording program and saves the events as a matrix in new
% files. Additionally, it generates a continuous timestamp (microseconds) for the events. If multiple sensors have been used 
% to record the data, then the program will produce a single file for every sensor.
%
% At the end of the code segment "split all events into several streams" you can choose which dataformat should be used for 
% the new files which contain the parsed events. By default they are saved as mat files. By default the new files will be 
% named like the input files, but with an appending "-parsed". This can be changed in the code section "parameters". 
% Additionally a number "-x" is added to the filename end to indicate which sensor generated these events.
% 
% !!! This function works only with the data fromat which uses 6 byte per event !!!
%
% ###########################################################################################################################
% 
% Input Parameters:
%
% "filename" :
% 
% String containing the file name which should be parsed. The file should be a txt file (the file exension can be changed in
% the code section "parameters") and the file name has to exclude the file suffix. For example for the file test.txt the 
% parameter "filename" should be 'test'.
%
% "streams" (default value 1) :
% 
% Optional parameter. If multiple sensor have been recorded at the same time, this parameters has to be set to the number of
% sensors.
%
% "save_information" (default value 0) :
%
% Optional parameter. Decides whether the plot is saved as matlab figure and as png file. Aditionally and the results which 
% are plotted into the command window are saved as an excel file.
% 0 = no ; 1 = yes 
%
% ###########################################################################################################################
%
% Output:
%
% The function does not return parameters, but it saves a single file for every sensor containing the output. These files 
% contain a matrix with the event data. Each line in the matrix contains one event and they are sorted by their timestamp. 
% Data format:
% [ x-position ; y-position ; polarity ; continuous timestamp (miscroseconds) ]
%
% Additionally there are some plots visualizing the original and continuous event timestamps.
%
% Furthermore 4 values are printed for every sensor:
%
% RealEventCount    : number of events captured with one sensor
% ElapsedTimeS      : time of last event in the record of one sensor (seconds)
% ActiveTime        : time between the first and the last event of the record of one sensor (seconds)
% EventRate         : event rate averaged over the whole recording (events/second)
%
% ###########################################################################################################################
%
% original author           : Jörg Conradt
% further development by    : Martin Medler 
% version                   : 1.1
% date                      : 14.10.2012


%% parameters

if nargin < 3
    save_information = 0;                   % decide if data characteristics are saved
end
if nargin < 2
    streams = 1;                            % number of cameras used to capture the eDVS128 data
end

fileSuffix      = '.txt';                   % file suffix for the input file
save_file_name  = [filename '-parsed'];     % name to save the new files

% various plot styles to distinguish the different sensors
PlotStyle{1} = 'b.';
PlotStyle{2} = 'r.';
PlotStyle{3} = 'g.';
PlotStyle{4} = 'y.';
PlotStyle{5} = 'c.';
PlotStyle{6} = 'k.';
PlotStyle{7} = 'm.';


%% create excel file (if user has chosen to save the plots)

if save_information == 1
    excelname = [save_file_name '.xls'];
    xlswrite(excelname, {'Data Block Errors'}, 1, 'A1');
    xlswrite(excelname, {'Event Count'}, 1, 'A2');
    xlswrite(excelname, {'Position n'}, 1, 'B2');
    xlswrite(excelname, {'Sync Errors'}, 1, 'D1');
    data_err_count = 4;
    sync_err_count = 4;
end


%% open file

f = fopen([filename fileSuffix], 'rb');        % open file
d = fread(f, inf, '*uint8');                   % read all data
fclose(f);
l = size(d,1);                                 % determine exact file length


%% translate binary data into "events"

% upper bound estimate of number of events
estUpperEventCount = ceil(l/7)    %#ok<*NOPRT>

% pre-allocate memory
eventS = zeros(estUpperEventCount,1)-1;          
eventX = zeros(estUpperEventCount,1);
eventY = zeros(estUpperEventCount,1);
eventP = zeros(estUpperEventCount,1);
eventT = zeros(estUpperEventCount,1);

n = 1;                          % index into binary data
c = 0;                          % counter of events extracted
currentSensor = 1;
eventCount=0;

while n<l                       % iterate over binary data
    if (d(n) < 10)              % special char, indicating data block from new sensor
        currentSensor = d(n);
        n = n+1;
        if (eventCount ~= 1001)
            disp(['Error! EventCount: ' num2str(eventCount)]);
            disp(['Position n = ' num2str(n)]);
            if save_information == 1
                xlswrite(excelname, eventCount, 1, ['A' num2str(data_err_count)]);
                xlswrite(excelname, n, 1, ['B' num2str(data_err_count)]);
                data_err_count = data_err_count + 1;
            end
        end
        eventCount = 0;        
    end
    
    c = c+1;                    % new event, remember data
    eventCount = eventCount + 1;    

    eventS(c) = currentSensor;
    eventX(c) = d(n)-32;          % x-position of event
    eventY(c) = d(n+1)-32;        % y-position of event
    eventP(c) = (bitand((d(n+2)-32), 64) == 64); % polarity of event
    eventT(c) = 128*128*128*double(bitand((d(n+2)-32), 15)) + 128*128*double(d(n+3)-32) + 128*double(d(n+4)-32) + double(d(n+5)-32);    % timestamp of event
    
    n = n+6;                      % advance counter by 5 (bytes read)

    if (n<l)
        if (d(n) ~= 10)
            disp(['sync error at byte position ' num2str(n)]);            
            % discard previous event
            eventS(c) = -1;
            eventX(c) = 0;
            eventY(c) = 0;
            eventP(c) = 0;
            eventT(c) = 0;
            c = c-1;
            while (n<l) && (d(n) ~= 10)
                n = n+1;
            end
            if save_information == 1
                xlswrite(excelname, n, 1, ['D' num2str(sync_err_count)]);
                sync_err_count = sync_err_count + 1;
            end
        end
    end

    n = n+1;                    % advance to next block start
end

%% prepare the figure for the plots

vis_fig_h = figure(1);
clf(vis_fig_h);
set(vis_fig_h, 'Position', [50, 50, 1400, 825], 'PaperPositionMode', 'auto');

plot_a = subplot(4,1,1);
hold on;
xlabel 'time (microseconds)'
ylabel 'time difference (microseconds)'
title 'time difference between successive events'
plot_b = subplot(4,1,2);
hold on;
xlabel 'event'
ylabel 'time difference (microseconds)'
title 'time difference between successive events'
plot_c = subplot(4,1,3);
hold on;
xlabel 'event'
ylabel 'timestamp (microseconds)'
title 'original timestamps of events'
plot_d = subplot(4,1,4);
hold on;
xlabel 'event'
ylabel 'timestamp (microseconds)'
title 'continuous timestamp of events'


%% split all events into several streams

RealEventCount = zeros(1, streams);
ElapsedTimeS = zeros(1, streams);
ActiveTime = zeros(1, streams);

for s = 1:(streams)
    eventXS = eventX(eventS==(s-1));
    eventYS = eventY(eventS==(s-1));
    eventPS = eventP(eventS==(s-1));
    eventTS = eventT(eventS==(s-1));    
    
    % fix timestamp
    timeOffset = 0;
    eventTSC = eventTS;
    for n=2:size(eventTS,1)
        if (eventTS(n) < eventTS(n-1))
            timeOffset = timeOffset + 32*1024*1024;
        end
        eventTSC(n) = eventTS(n)+timeOffset;
    end
    TimeDiff = eventTSC(2:end) - eventTSC(1:end-1);    
    
    % plot data
    plot(plot_a, eventTSC(1:end-1), TimeDiff, PlotStyle{s})
    plot(plot_b, TimeDiff, PlotStyle{s})
    plot(plot_c, eventTS, PlotStyle{s})
    plot(plot_d, eventTSC, PlotStyle{s})
    
    % get characteristics of data stream
    RealEventCount(s) = size(eventTS, 1);
    ElapsedTimeS(s) = max(eventTSC)/1000000;
    ActiveTime(s) = (max(eventTSC) - min(eventTSC))/1000000;
    
    % save data
    e = [eventXS, eventYS, eventPS, eventTSC];     %#ok<NASGU>
    save([save_file_name '-' num2str(s) '.mat'], 'e', '-mat');
%     save([save_file_name '-' num2str(s) '.txt'], 'e', '-ascii');
    
end

linkaxes([plot_b plot_c plot_d], 'x');


%% save results of data analasis

if save_information == 1    
    xlswrite(excelname, {'Event Count'}, 1, 'F1');
    xlswrite(excelname, RealEventCount, 1, 'F3');
    xlswrite(excelname, {'Elapsed Time'}, 1, 'F5');
    xlswrite(excelname, round(ElapsedTimeS*100)/100, 1, 'F7');
    xlswrite(excelname, {'Active Time'}, 1, 'F10');
    xlswrite(excelname, round(ActiveTime*100)/100, 1, 'F12');
    xlswrite(excelname, {'Event Rate'}, 1, 'F15');
    xlswrite(excelname, round(RealEventCount./ActiveTime), 1, 'F17');    
    
    saveas(vis_fig_h, save_file_name, 'fig');
    saveas(vis_fig_h, save_file_name, 'png');
end

%% show characteristics of data

fprintf('\ndata characteristics for all streams \n')

RealEventCount  %#ok<NOPRT>
ElapsedTimeS    %#ok<NOPRT>
ActiveTime      %#ok<NOPRT>
EventRate = RealEventCount/ActiveTime  %#ok<NASGU>


