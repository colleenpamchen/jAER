function Calib_Distortion_Paras( Events_Orig )
% Function to manually find barrel distortion parameters of eDVS128 data. 
%
% This function is used to find the distortion parameters for eDVS128 event data with a trial and error approach. 
% The function displays a figure with 2 videos. On the left side the unchanged event data and on the right side as comparison
% the event data which has been changed accordingly to the currently chosen distortion parameters. These videos are played 
% endlessly. 
% A second figure contains a small GUI which allows to change the distortion parameters live on the fly. The right video
% containig the corrected data is adjusted online to each distortion parameter change.
%
% On events (polarity 1) are visualized as green dots and off events ( polarity 0) are visualized as red dots.
% The currently chosen distortion parameters are printed in the command window if you end the function with the key "q". Do 
% not forget to copy them down before closing matlab.
% It is possible to control the function online with the keyboard, which is
% explained below.
% The code section "parameters" contains several variables which can be used to customize the function to personal needs.
% 
% ###########################################################################################################################
%
% Input Parameters:
%
% "Events_Orig" :
%
% Matrix containing the event data. Each line should represent one event and  the events should be ordered by time. The 
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
% If the distortion parameters are changed during a paused frame, the frame will be updated accrodingly. 
%
% ###########################################################################################################################
%
% GUI:
%
% The GUI contains 3 slider bars, which allow to change the distortion parameters. The text box above the slider shows the 
% currently chosen value. 
% Below the slider are two text boxes "med impact" and "max impact". These tell you how much this parameter will change the 
% event coordinates. The "max impact" shows you how much the event coordinate will be increased, if the event has the maximum
% distance from the coordinate center. "Med impact" shows you the same thing for half the maximum distance. There is no
% distinction between the change for x and y because we assume for these values, that the event coordinates lie on the 
% line x = y.
% Pay attention to release the mouse button over the slider, while manipulating the slider bar. Clicking the slider, but 
% releasing the mouse button somewhere else fails to update the slider value.
%
% ###########################################################################################################################
%
% General Remarks:
%
% Use recordings with clear straight edges which should be evenly distributed through the whole image (either by recording 
% a big pattern, or recording a single edge and moving the sensor, so the edge travels through the whole image). This allows 
% you to clearly see the distortion and if the distortion correction works.
% Try to find a solution with only k1 and k2. Use k3 only if is really neccessary. More variables means more work to find a 
% solution.
% k1 has a linear, k2 a quadratic and k3 a cubic impact on the event coordinates depending on the distance between the event 
% and the coordiante center. Therefore, if the middle of the image looks good, but the edges are still distorted, try to 
% change a higher k, like k2.
%
% ###########################################################################################################################
%
% author            : Martin Medler 
% version           : 1.1
% date              : 14.10.2012


%% load input data

if ischar(Events_Orig)
    Events_Orig = importdata(Events_Orig);
end


%% parameters

% min and max values for the distortion parameters available via the slider
min_k1 = 0;
max_k1 = 0.02;
min_k2 = 0;
max_k2 = 0.0002;
min_k3 = 0;
max_k3 = 0.000005;

% center of the original coordinate system (127/2)
global coord_mid
coord_mid = 63.5;

% time interval which is visulalized as one image
vis_delta = 0.05 * 1000000;                                                 % 20 frames/second

global k1;
global k2;
global k3;
% inital values for distortion parameters
k1 = 0.0013;
k2 = 0.000016;
k3 = 0.0;

% process events with inital distortion parameter values
Events_Undist = CorrectDistortion(Events_Orig, coord_mid, k1, k2, k3);
% obtain required image resolution and corresponding offset (used for plotting)
[resol_undist, offset_undist] = GetResolution(Events_Undist, coord_mid);
[resol_orig, offset_orig] = GetResolution(Events_Orig, coord_mid);

% colormap to show discrete colours:
% [ red ; yellow ; black ; blue ; green ]
% [ -5  ;   -2   ;   0   ;   2  ;    5  ]
myCMap = [ 1 , 0 , 0 ;
    1 , 1 , 0 ;
    0 , 0 , 0 ;
    0 , 0 , 1 ;
    0 , 1 , 0 ];


%% figure for visualization

vis_fig_h           = figure(1);
clf(vis_fig_h);
vis_axes_orig_h     = subplot(1,2,1);
vis_axes_undist_h   = subplot(1,2,2);
set(vis_fig_h, 'Colormap', myCMap, 'WindowKeyPressFcn',@KeyContrProg);


%% variables for program control

global quitProg
global stepByStep
global pressedButton
global recalcEvents;
quitProg = 0;
stepByStep = 0;
pressedButton = 0;
recalcEvents = 0;


%% GUI for parameters

contr_fig_h = figure(2);
clf(contr_fig_h);

hSliderK1 = uicontrol(contr_fig_h, 'style', 'text', 'Position', [190,385,190,20], ... 
    'FontSize', 10, 'String', ['K1  :  ' num2str(k1, '%.8f')]);
hMedImpactK1 = uicontrol(contr_fig_h, 'style', 'text', 'Position', [70,320,190,20], ... 
    'FontSize', 10, 'String', ['med impact K1  :  ' num2str( sqrt( 2*(coord_mid/2)^2 )*k1*coord_mid/2, '%.3g' )]);
hMaxImpactK1 = uicontrol(contr_fig_h, 'style', 'text', 'Position', [320,320,190,20], ... 
    'FontSize', 10, 'String', ['max impact K1  :  ' num2str( sqrt( 2*(coord_mid)^2 )*k1*coord_mid, '%.3g' )]);
uicontrol(contr_fig_h, 'style', 'slider', 'Min', min_k1, 'Max', max_k1, 'SliderStep', [0.005,0.1], ...
    'Value', k1, 'Position', [135,355,300,17], 'Callback', {@SliderK1, hSliderK1, hMedImpactK1, hMaxImpactK1, coord_mid});

hSliderK2 = uicontrol(contr_fig_h, 'style', 'text', 'Position', [190,240,190,20], ...
    'FontSize', 10, 'String', ['K2  :  ' num2str(k2, '%.8f')]);
hMedImpactK2 = uicontrol(contr_fig_h, 'style', 'text', 'Position', [70,175,190,20], ...
    'FontSize', 10, 'String', ['med impact K2  :  ' num2str( sqrt( 2*(coord_mid/2)^2 )^2*k2*coord_mid/2, '%.3g' )]);
hMaxImpactK2 = uicontrol(contr_fig_h, 'style', 'text', 'Position', [320,175,190,20], ... 
    'FontSize', 10, 'String', ['max impact K2  :  ' num2str( sqrt( 2*(coord_mid)^2 )^2*k2*coord_mid, '%.3g' )]);
uicontrol(contr_fig_h, 'style', 'slider', 'Min', min_k2, 'Max', max_k2, 'SliderStep', [0.005,0.1], ...
    'Value', k2, 'Position', [135,210,300,17], 'Callback', {@SliderK2, hSliderK2, hMedImpactK2, hMaxImpactK2, coord_mid});

hSliderK3 = uicontrol(contr_fig_h, 'style', 'text', 'Position', [190,95,190,20], ...
    'FontSize', 10, 'String', ['K3  :  ' num2str(k3, '%.8f')]);
hMedImpactK3 = uicontrol(contr_fig_h, 'style', 'text', 'Position', [70,30,190,20], ...
    'FontSize', 10, 'String', ['med impact K3  :  ' num2str( sqrt( 2*(coord_mid/2)^2 )^3*k3*coord_mid/2, '%.3g' )]);
hMaxImpactK3 = uicontrol(contr_fig_h, 'style', 'text', 'Position', [320,30,190,20], ...
    'FontSize', 10, 'String', ['max impact K3  :  ' num2str( sqrt( 2*(coord_mid)^2 )^3*k3*coord_mid, '%.3g' )]);
uicontrol(contr_fig_h, 'style', 'slider', 'Min', min_k3, 'Max', max_k3, 'SliderStep', [0.005,0.1], ...
    'Value', k3, 'Position', [135,65,300,17], 'Callback', {@SliderK3, hSliderK3, hMedImpactK3, hMaxImpactK3, coord_mid});


%% iterate endlessly over events and visualize them

start_time = Events_Orig(1,4);

while 1    
    % initilaize parameters for next run
    last_shown_event = 0;
    next_vis = Events_Orig(1,4) + vis_delta;
    tic
    
    for k = 1:1:size(Events_Undist, 1)        
        % check if user aborted the program
        if quitProg == 1
            disp('program terminated by user')
            fprintf('\nk1 : %.8f \nk2 : %.8f \nk3 : %.8f \n\n' , k1, k2, k3);
            return
        end
        
        % check if events have to be recalculated with new distortion parameters
        if recalcEvents == 1
            recalcEvents = 0;
            Events_Undist = CorrectDistortion(Events_Orig, coord_mid, k1, k2, k3);
            [resol_undist, offset_undist] = GetResolution(Events_Undist, coord_mid);
        end
        
        % limit event processing rate to fit event timestamp
        while toc * 1000000 < Events_Undist(k, 4) - start_time
        end
        
        % visualize events in a regular timer interval
        if Events_Orig(k,4) > next_vis
            next_vis = next_vis + vis_delta;
            
            % plot events
            EventPointer = (last_shown_event+1:k);           
            PlotImage(Events_Orig, EventPointer, offset_orig, resol_orig, vis_axes_orig_h);            
            PlotImage(Events_Undist, EventPointer, offset_undist, resol_undist, vis_axes_undist_h);            
            drawnow;
            last_shown_event = k;
            
            % in pause mode -> continue porogram after a key has been pressed
            if stepByStep == 1                
                pressedButton = 0;                
                while pressedButton == 0
                    % check if events have too be recalculated and
                    % visualized during pause
                    if recalcEvents == 1
                        recalcEvents = 0;
                        Events_Undist = CorrectDistortion(Events_Orig, coord_mid, k1, k2, k3);
                        [resol_undist, offset_undist] = GetResolution(Events_Undist, coord_mid);                        
                        PlotImage(Events_Undist, EventPointer, offset_undist, resol_undist, vis_axes_undist_h);                        
                        drawnow;
                    end                    
                    pause(0.01);
                end
            end            
        end        
    end
end

end


%% distortion correction function

function Events = CorrectDistortion(Events, coord_mid, k1, k2, k3)

R = sqrt( (Events(:,1)-coord_mid).^2 + (Events(:,2)-coord_mid).^2 );
Events(:,1) = ( Events(:,1)-coord_mid ) .* ( 1 + k1*R + k2*R.^2 + k3*R.^3 ) + coord_mid;
Events(:,2) = ( Events(:,2)-coord_mid ) .* ( 1 + k1*R + k2*R.^2 + k3*R.^3 ) + coord_mid;

end


%% function to plot images

function PlotImage(Events, EventPointer, offset, resolution, axes_handle)
Image = zeros(resolution);

X = ceil( Events(EventPointer,1) ) + offset;
Y = ceil( Events(EventPointer,2) ) + offset;
PrintPointer = ((X-1) * resolution) + Y;

Image(PrintPointer) = Events(EventPointer,3)*10 - 5;

imagesc(Image, 'Parent', axes_handle);
title(axes_handle, [num2str( Events(EventPointer(end), 4)/1000000 ) 's (Event-Time)']);
axis(axes_handle, 'image');
caxis(axes_handle, [-5, 5]);

end


%% function to obtain image resolution fitting the evant data

function [resolution, offset] = GetResolution(Events, coord_mid)
max_ev_coord = zeros(1,4);

max_ev_coord(1:2) = abs(min(Events(:,1:2)-coord_mid));
max_ev_coord(3:4) = max(Events(:,1:2)-coord_mid);
resolution = ceil(max(max_ev_coord))*2 + 6;
offset = floor( (resolution-2*coord_mid)/2 ) + 1;

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


%% GUI callback functions

function SliderK1(hObj,~, k1_text_h, med_imp_h, max_imp_h, coord_mid)
global recalcEvents;
global k1;

recalcEvents = 1;
k1 = get(hObj,'Value');

set(k1_text_h, 'String', ['K1: ' num2str(k1, '%.8f')]);
set(med_imp_h, 'String', ['max impact K1  :  ' num2str( sqrt( 2*(coord_mid/2)^2 )*k1*coord_mid/2, '%.3g' )]);
set(max_imp_h, 'String', ['max impact K1  :  ' num2str( sqrt( 2*(coord_mid)^2 )*k1*coord_mid, '%.3g' )]);

end

function SliderK2(hObj,~, k2_text_h, med_imp_h, max_imp_h, coord_mid)
global recalcEvents;
global k2;

recalcEvents = 1;
k2 = get(hObj,'Value');

set(k2_text_h, 'String', ['K2: ' num2str(k2, '%.8f')]);
set(med_imp_h, 'String', ['med impact K2  :  ' num2str( sqrt( 2*(coord_mid/2)^2 )^2*k2*coord_mid/2, '%.3g' )]);
set(max_imp_h, 'String', ['max impact K2  :  ' num2str( sqrt( 2*(coord_mid)^2 )^2*k2*coord_mid, '%.3g' )]);

end

function SliderK3(hObj,~, k3_text_h, med_imp_h, max_imp_h, coord_mid)
global recalcEvents;
global k3;

recalcEvents = 1;
k3 = get(hObj,'Value');

set(k3_text_h, 'String', ['K3: ' num2str(k3, '%.8f')]);
set(med_imp_h, 'String', ['med impact K3  :  ' num2str( sqrt( 2*(coord_mid/2)^2 )^3*k3*coord_mid/2, '%.3g' )]);
set(max_imp_h, 'String', ['max impact K3  :  ' num2str( sqrt( 2*(coord_mid)^2 )^3*k3*coord_mid, '%.3g' )]);

end
