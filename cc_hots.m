%% PRE-RECORDED DVS DATASET from Jorg Conradt 
% data already in matrix format: 
    % Column 1: x coordinate (from 0 to 127) 
    % Column 2: y coordinate (from 0 to 127)
    % Column 3: event polarity [0 off | 1 on]
    % Column 4: timestamps with 1us time tick (microseconds?miliseconds? how to convert to seconds?)
    % Events = [ x-coordinate , y-coordinate , polarity , timestamp (microseconds) ]
Events = load('pen_vertical.dvs');
% Events = load('pen_horizontal.dvs');
% Events = load('spinner.dvs');
% Events = load('hand_vertical.dvs');
% Events = load('hand_horizontal.dvs');
    
%% video visualizing the events 
visualize_events(Events); 

%% Batch incoming events: 
num_events = size(Events,1) ; 

batchsize = 1000;
   idx=[1:batchsize];
for i = 1:100
    data = Events(idx,:); 
    idx = idx+batchsize
end
