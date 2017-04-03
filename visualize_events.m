function Visualize_Events( Events )
% Visualizes eDVS128 event data in a video.
%
% This function visualizes eDVS128 event data as video within a matlab figure. It does not matter if the event data 
% coordinates are integers ranging from 0-127 (e.g. unchanged data format of sensor) or if the values are double values in a 
% larger range (e.g. after distortion correction). The center of the event coordinate system should be 63.5, but this 
% parameter can be adjusted.
% On events (polarity 1) are visualized as green dots and off events ( polarity 0) are visualized as red dots.
% It is possible to control the function online with the keyboard, which is explained below.
% The code section "parameters" contains several parameters, which can be used to adjust the function behavior.
%
% ###########################################################################################################################
%
% Input Parameters:
%
% "Events" :
%
% Matrix containing the event data. Each line should represent one event and the events should be ordered by time. The 
% timestamp needs to be a continuous timestamp. Data format for one line:
% 
% [ x-position ; y-position ; polarity ; continuous timestamp (microseconds) ]
%
% This parameter can also be a string containing the path to a file containing the event data which will be loaded.
%
% ###########################################################################################################################
%
% Keyboard Control:
%
% To use the keyboard control, the video figure has to be in focus.
% Pressing "q" ends the program. 
% Pressing "p" switches the program into step by step mode. Pressing "p" agains switches back to the normal live mode. In the
% step by step mode the function is paused after a video frame has been plotted and the next frame will not be plotted, until
% a key has been pressed.
%
% ###########################################################################################################################
%
% author        : Martin Medler
% version       : 1.1
% date          : 14.10.2012

%% load input data

if ischar(Events)
    Events = importdata(Events);
end


%% parameters

% center of the event coordinate system (127/2)
coord_mid = 63.5;

% time interval which is visulalized as one image
vis_delta = 0.05 * 1000000;                                                % 20 frames/second

% obtain resolution and corresponding offset for event data
max_ev_coord = zeros(1,4);
max_ev_coord(1:2) = abs(min(Events(:,1:2)-coord_mid));
max_ev_coord(3:4) = max(Events(:,1:2)-coord_mid);
resolution = ceil(max(max_ev_coord))*2 + 6;
offset = floor( (resolution-2*coord_mid)/2 ) + 1;

% colormap to show discrete colours:
% [ red ; yellow ; black ; blue ; green ]
% [ -5  ;   -2   ;   0   ;   2  ;    5  ]
myCMap = [ 1 , 0 , 0 ;
    1 , 1 , 0 ;
    0 , 0 , 0 ;
    0 , 0 , 1 ;
    0 , 1 , 0 ];


%% variables for program control

global quitProg
global stepByStep
global pressedButton
quitProg = 0;
stepByStep = 0;
pressedButton = 0;


%% figure for visualization

vis_fig_h = figure(1);
clf(vis_fig_h);
vis_axes_h = axes();
set(vis_fig_h, 'Colormap', myCMap, 'WindowKeyPressFcn', @KeyContrProg);


%% iterate over events and visualize them

Image               = zeros(resolution);
next_vis            = Events(1,4) + vis_delta;
start_time          = Events(1,4);
last_shown_event    = 0;
tic

for k = 1:1:size(Events, 1)    
    % check if user aborted the program
    if quitProg == 1
        disp('program terminated by user')
        return
    end
    
    % limit event processing rate to fit event timestamp
    while toc * 1000000 < Events(k, 4) - start_time
    end
    
    % visualize events in a regular timer interval
    if Events(k, 4) > next_vis
        next_vis = next_vis + vis_delta;
                
        % write events into image matrix
        Image(:) = 0;
        EventPointer = (last_shown_event+1:k);        
        X = ceil( Events(EventPointer,1) ) + offset;
        Y = ceil( Events(EventPointer,2) ) + offset;
        PrintPointer = ((X-1) * resolution) + Y;        
        Image(PrintPointer) = Events(EventPointer,3)*10 - 5;
        
        % plot image
        imagesc(Image, 'Parent', vis_axes_h);
        title(vis_axes_h, [num2str( Events(EventPointer(end), 4)/1000000 ) 's (Event-Time)']);
        axis(vis_axes_h, 'image');
        caxis(vis_axes_h, [-5, 5]);
        drawnow;
        last_shown_event = k;

        % wait for key press (program control)
        if stepByStep == 1
            pressedButton = 0;
            while pressedButton == 0
                pause(0.01);
            end
        end
    end
end

end


%% key control function

function KeyContrProg(~,eventData)
global quitProg
global stepByStep
global pressedButton

pressedButton = 1;

if eventData.Character == 'p'
    if stepByStep == 0
        stepByStep = 1;
    else
        stepByStep = 0;
    end
end

if eventData.Character == 'q'
    quitProg = 1;
end

end
