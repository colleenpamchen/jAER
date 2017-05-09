function Analyze_Event_Timediffs(filename, sensors, save)
% Analyze the time difference between successive events to search for recording errors.
% 
% This function analyzes the time differences between successive events in a set of eDVS128 data. It generates a logarithmic
% plot visualizing these time differences and counts the occurences of big time differences. Additionally it prints several
% data characeristics into the command window (more infos in the section Output).
%
% It counts 4 kinds of errors. Time differences which are:
% - greater than 10^3 and smaller than 10^4 microseconds
% - greater than 10^4 and smaller than 10^5 microseconds
% - greater than 10^5 and smaller than 10^6 microseconds
% - greater than 10^6 microseconds
%
% ###########################################################################################################################
%
% Input Parameters:
%
% "filename" :
% 
% String containing the file name of the data files. The strng has to exclude the file identifyer and the number identifying
% the sensor. The files have to be mat files. If the recording has been done with several sensors, the function expects a 
% single file for every sensor and each file with a number at the end of the filename which identifies the sensor which 
% recorded these events. 
% For example if you recorded data with 2 sensors and generated the files "test-parsed-1.mat" and  "test-parsed-2.mat", then
% you should choose 'test-parsed' as "filename" and a value of 2 for the input variable "sensors".
%
% "sensors" (default value 1) :
% 
% Optional parameter. If multiple sensor have been recorded at the same time, this parameters has to be set to the number of
% sensors.
%
% "save" (default value 0) :
%
% Optional parameter. Decides whether the plot is saved as matlab figure and as png file. Aditionally and the results which 
% are plotted into the command window are saved as an excel file.
% 0 = no ; 1 = yes 
%
%% ###########################################################################################################################
%
% Output:
%
% The function generates a plot which visualizes the time differences between the events.
%
% Furthermore several values are printed for every sensor:
%
% RealEventCount    : Number of events captured with one sensor
% ActiveTime        : Time between the first and the last event of the record of one sensor (seconds)
% EventRate         : Event rate averaged over the whole recording (events/second)
% Error Count       : Total count for the 4 error categories for all sensors
% Error             : Matrix containing the detailed error counts. Every column represents a sensor and every line an error
%                     category
%
% ###########################################################################################################################
%
% author        : Martin Medler
% version       : 1.2
% date          : 20.10.2012


%% check input

if nargin < 3
    save = 0;
end
if nargin < 2
    sensors = 1;
end


%% parameters

% error thresholds to analyze data quality
err_thresh_1 = 10^3;
err_thresh_2 = 10^4;
err_thresh_3 = 10^5;
err_thresh_4 = 10^6;

% upper limit for error plot
y_thresh = 10^7;

% various plot styles to distinguish different sensors
PlotStyle{1} = 'b.';
PlotStyle{2} = 'r.';
PlotStyle{3} = 'g.';
PlotStyle{4} = 'm.';
PlotStyle{5} = 'c.';
PlotStyle{6} = 'k.';
PlotStyle{7} = 'y.';

% filename to save the figures and the exsel file
save_ana_name  = [filename '-ana'];     


%% create figure

vis_fig_h = figure(1);
clf(vis_fig_h);
set(vis_fig_h, 'Position', [50, 50, 1500, 825], 'PaperPositionMode', 'auto');

plot_a = subplot(2,1,1);
hold on;
xlabel 'time (microseconds)'
ylabel 'time difference (microseconds)'
title 'time difference between successive events'
plot_b = subplot(2,1,2);
hold on;
xlabel 'events'
ylabel 'time difference (microseconds)'
title 'time difference between successive events'


%% analyze data

% allocate memory
Event_Rate  = zeros(1, sensors);
Numb_Events = zeros(1, sensors);
Act_Time    = zeros(1, sensors);
Elap_Time   = zeros(1, sensors);
Errors      = zeros(4, sensors);
Err_Count   = zeros(1, 4);

for k = 1:sensors
    % read mat file containing events
    Input = load([filename '-' num2str(k) '.mat'], '-mat');
    names = fieldnames(Input);
    if length(names) == 1
        Time = Input.(names{1})(:,4);
        clearvars Input;
    else
        fprintf(2, 'error: more then one variable within mat file\n\n');
        return
    end
    
    % basic parameters
    Numb_Events(k)  = length(Time);
    Elap_Time(k)    = Time(end);
    Act_Time(k)     = ( Time(end)-Time(1) )/1000000;
    Event_Rate(k)   = Numb_Events(k)/Act_Time(k);
    
    % time differences
    Time_Diff = Time(2:end) - Time(1:end-1);
    plot(plot_a, Time(1:end-1), Time_Diff, PlotStyle{k});
    plot(plot_b, Time_Diff, PlotStyle{k});

    % find errors (big time difference between events)
    Errors(1, k) = length(find( (Time_Diff > err_thresh_1) & (Time_Diff < err_thresh_2) ));
    Errors(2, k) = length(find( (Time_Diff > err_thresh_2) & (Time_Diff < err_thresh_3) ));
    Errors(3, k) = length(find( (Time_Diff > err_thresh_3) & (Time_Diff < err_thresh_4) ));
    Errors(4, k) = length(find( (Time_Diff > err_thresh_4) ));
    
    % count errors
    Err_Count(1) = Err_Count(1) + length(find( (Time_Diff > err_thresh_1) & (Time_Diff < err_thresh_2) ));
    Err_Count(2) = Err_Count(2) + length(find( (Time_Diff > err_thresh_2) & (Time_Diff < err_thresh_3) ));
    Err_Count(3) = Err_Count(3) + length(find( (Time_Diff > err_thresh_3) & (Time_Diff < err_thresh_4) ));
    Err_Count(4) = Err_Count(4) + length(find( (Time_Diff > err_thresh_4) ));    
end


%% format plot

% plot legend
Legend = cell(1, sensors);
for k = 1:sensors
    Legend{k} = ['sensor ' num2str(k)];
end

% create YTicks for logarithmic y-axis 
Log_Ticks = ones(1, round(log(y_thresh)+5));
for k = 2:length(Log_Ticks)
    Log_Ticks(k) = Log_Ticks(k-1) * 10;
end

% scale plot
set(plot_a, 'YScale', 'log', 'YTick', Log_Ticks);
axis(plot_a, [0, max(Elap_Time), 10^0, y_thresh]);
legend(plot_a, Legend, 'Location', 'NorthEastOutside');
set(plot_b, 'YScale', 'log', 'YTick', Log_Ticks);
axis(plot_b, [0, max(Numb_Events), 10^0, y_thresh]);
legend(plot_b, Legend, 'Location', 'NorthEastOutside');

drawnow


%% save plots and excel file (if user has chosen to)

% save results
if save == 1
    excelname = [save_ana_name '.xls'];
    
    xlswrite(excelname, {'Events'}, 1, 'A1');
    xlswrite(excelname, Numb_Events, 1, 'B2');    
    xlswrite(excelname, {'Act. Time'}, 1, 'A4');
    xlswrite(excelname, round(Act_Time*100)/100, 1, 'B5');    
    xlswrite(excelname, {'Event Rate'}, 1, 'A7');
    xlswrite(excelname, round(Event_Rate), 1, 'B8');    
    xlswrite(excelname, {'Errors'}, 1, 'A10');
    xlswrite(excelname, {'10^3 - 10^4'}, 1, 'B11');
    xlswrite(excelname, {'10^4 - 10^5'}, 1, 'B12');
    xlswrite(excelname, {'10^5 - 10^6'}, 1, 'B13');
    xlswrite(excelname, {' > 10^6'}, 1, 'B14');
    xlswrite(excelname, Errors, 1, 'D11');
    xlswrite(excelname, Err_Count', 1, 'K11');  
    
    saveas(vis_fig_h, save_ana_name, 'fig');
    saveas(vis_fig_h, save_ana_name, 'png');
end

%% check if plot contains all errors

if max(Time_Diff) > y_thresh
    fprintf(2,'\n\n!!! diff value out of border !!!\n\n');
    fprintf('max time diff between successive events : %d\n', max(Time_Diff));
end


%% print results into command window
fprintf('\nNumber Events:\n\n');
fprintf('%-10d ', Numb_Events);
fprintf('\n\nActive Time:\n\n');
fprintf('%-10.1f ', Act_Time);
fprintf('\n\nEvent Rate:\n\n');
fprintf('%-10.0f ', Event_Rate);
fprintf('\n\n\nError Count:\n\n');
fprintf('%-6d ', Err_Count);
fprintf('\n')
Errors %#ok<NOPRT>

end