%% colleen's working progress 
% [CIN]=loadaerdat_matlab('DVS128-2017-03-21_penOrientation.aedat');
[CIN]=loadaerdat_matlab('DVS128-2017-03-21T16-00-20-0700-0.aedat');

% Column 1: timestamps with 1us time tick
% Columns 2-3: ignore them (they are meant for simulator AERST (see Perez-Carrasco et al, IEEE
% TPAMI, Nov. 2013)).
% Column 4: x coordinate (from 0 to 127)
% Column 5: y coordinate (from 0 to 127)
% Column 6: event polarity
% [ x-coordinate , y-coordinate , polarity , timestamp (microseconds) ]

n=size(CIN,1);
nn = ceil(n / 100) ; 

ts= CIN(:,1); 
xx=CIN(:,4);
yy=CIN(:,5);
xy=CIN(:,4:5); 
P=CIN(:,6);

Son = nan(128);
Soff = nan(128);

deltaT = zeros(128);
ts_prev=0; 
tau = 50000 ; % from paper: 50 milliseconds, conver this into microseconds. 
% this loops through the entire recording. 
% TODO: change it so that it takes incoming data one by one 

for i = 1: n
    
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
%    imagesc(S)
%    F(i)= getframe; 
   
end

surf(Son) 


% movie(F)

% imagesc(S)




























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