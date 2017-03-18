%% colleen's working progress 
[CIN]=loadaerdat_matlab('DVS.aedat');

% Column 1: timestamps with 1us time tick
% Columns 2-3: ignore them (they are meant for simulator AERST (see Perez-Carrasco et al, IEEE
% TPAMI, Nov. 2013)).
% Column 4: x coordinate (from 0 to 127)
% Column 5: y coordinate (from 0 to 127)
% Column 6: event polarity

n=size(CIN,1);

ts= CIN(:,1); 
xx=CIN(:,4);
yy=CIN(:,5);
xy=CIN(:,4:5); 
P=CIN(:,6);

S = NaN(128);
deltaT = zeros(128);
ts_prev=0; 
tau = 1 ; 

% this loops through the entire recording. 
% TODO: change it so that it takes incoming data one by one 

% x=xx(1)+1; % fix MATLAB's indexing starting with 1 
% y=yy(1)+1; 
% deltaT(x,y) = ts(1);
% S(x,y) = exp(deltaT(x,y)) ;
% S(x,y) = exp( (-ts(i)+ deltaT(x,y))/tau )  ;
   
for i = 1:n
    
   x=xx(i)+1; % fix MATLAB's indexing starting with 1 
   y=yy(i)+1; 
   ts_current = ts(i); 
   deltaT = deltaT + (ts_current - ts_prev);
   ts_prev=ts_current; 
   deltaT(x,y)=0; 
   S= exp( -(deltaT)/tau )  ; 
 
    % for a location (x,y), take the current time at THAT location minus previous time at that location  
       
end

imagesc(S)



























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