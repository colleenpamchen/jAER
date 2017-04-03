%% colleen's working progress 
% [CIN]=loadaerdat_matlab('DVS128-2017-03-21_penOrientation.aedat');
% [CIN]=loadaerdat_matlab('DVS128-2017-03-21T16-00-20-0700-0.aedat');
[CIN]=loadaerdat_matlab('DVS128-2017-03-21-short.aedat');

% Column 1: timestamps with 1us time tick
% Columns 2-3: ignore them (they are meant for simulator AERST (see Perez-Carrasco et al, IEEE
% TPAMI, Nov. 2013)).
% Column 4: x coordinate (from 0 to 127)
% Column 5: y coordinate (from 0surf(Son) to 127)
% Column 6: event polarity
% Events = [ x-coordinate , y-coordinate , polarity , timestamp (microseconds) ]

n = size(CIN,1);
ts = CIN(:,1); 
x0 = CIN(:,4);
y0 = CIN(:,5);
xy0 = CIN(:,4:5); 
P = CIN(:,6);
    
%%
for k=1:size(P,1)
    if P(k) == -1
       P(k) = 0;
    end
end
logical_pol = logical(P);

x_incr = x0(logical_pol);
y_incr = y0(logical_pol);
ts_incr = ts(logical_pol);

x_decr = x0(~logical_pol);
y_decr = y0(~logical_pol);
ts_decr = ts(~logical_pol);

res_x = 128;
res_y = 128;

% epg (events per graph) = it is associated with the scatterd mode. It
% represents the number of events displayed in each image. Default: 20000 
epg = 20000;


for i=ts(1): epg :ts(end)
    figure; 
    %it creates a mask to identify the events with polarity = 1, in the time interval we are going to represent
    incr_mask = ts_incr >= i & ts_incr < i + epg;
    %It plots events with polarity = 1
    scatter3(x_incr(incr_mask), ts_incr(incr_mask), y_incr(incr_mask), 1, 'g');
    hold on
    %it creates a mask to identify the events with polarity = 0, in the time interval we are going to represent
    decr_mask = ts_decr >= i & ts_decr < i+epg;
    %It plots events with polarity = 0
    scatter3(x_decr(decr_mask), ts_decr(decr_mask), y_decr(decr_mask), 1, 'r');
    %axes formatting
        axis([0 res_x i i+epg-1 0 res_y]);
        strx = sprintf('X Pixel [0 %f]',res_x);
        xlabel(strx);
        ylabel('Timestamps [\mus]');
        strz = sprintf('Y Pixel [0 %f]',res_y);
        zlabel(strz);
    
end

    

%%

% Events = zeros(n,4);
Events = [x0,y0,P,ts];
visualize_events(Events); 

% center of the original coordinate system (128/2) = 64
coord_mid = 64;
R = 2 ; 

[64-2 , 64+2] 

xmask = ;
ymask = ;



Son = nan(128);
Soff = nan(128);
deltaT = zeros(128);
ts_prev=0; 
tau = 20000 ;
% tau = 50000 ; % from paper: 50 milliseconds, convert this into microseconds. 
% this loops through the entire recording. 
% TODO: change it so that it takes incoming data one by one 
for i = 1:n
    
   x=xx(i)+1; % fix MATLAB's indexing starting with 1 
   y=yy(i)+1; 
   ts_current = ts(i); 
   
   if P(i)==1 % for ON polarity events only 
       deltaT = deltaT + (ts_current - ts_prev);
       ts_prev = ts_current; 
       deltaT(x,y)=0; 
       Son = exp( -(deltaT)/tau )  ; 
      % for a location (x,y), take the current time at THAT location minus previous time at that location  
   end  
   
%    if P(i)==-1 % for OFF polarity events only 
%        deltaT = deltaT + (ts_current - ts_prev);
%        ts_prev=ts_current; 
%        deltaT(x,y)=0; 
%        Soff= exp( -(deltaT)/tau )  ; 
%       % for a location (x,y), take the current time at THAT location minus previous time at that location  
%    end 

%    imagesc(Son)
   surf(Son)
   pause(1/2)

   
end

surf(Son) 


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

























% This function should take the following arguments: Xt,Yt,tau
% 
% Xt is a 2-D array of spatial locations from a single polarity (ON/OFF) channel in the DAVIS camera.  
% Yt is the output of the function from the previous timestep.  
% tau is a system parameter that will determine the how the filter shapes its output over time. 
% 
% It should provide the following outputs: Y(t+1)
% 
% Y(t+1) is the result of combining Xt with Yt.  
% I believe if tau = 1, then Y(t+1) = (Xt + Yt)/2.... i.e. the time-average of inputs.  
% To produce Y(t+1) for arbitrary values of tau, use the function Benosman et al. refer to as S(x,y,t)