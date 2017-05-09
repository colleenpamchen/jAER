function [Events] = Correct_Barrel_Distortion(Events, k1, k2, k3, coord_mid)
% Correct barrel distortion of eDVS128 data.
%
% This function corrects the barrel distortion in eDVS128 data and returns an array with the corrected data. 
% Only the parameter "Events" is mandatory. However, without at least a value for "k1", the function will not change the 
% data.
%
% This function assumes that the barrel distortion is symetric and centered at the middle of the lens. A barrel distortion is
% stronger near the image border and weak in the image center and therefore the function computes the distance of every event
% to the distortion center (center of the image) and then uses this value  to adjust the x and y coordinates with the 
% following formula (exemplary for x):
%
% r = sqrt( x^2 + y^2 ) 
% x_new = (x_new - x_center) * (1 + k_1*r + k_2*r^2 + k_3*r^3 + ...) + x_center
%
% where "r" is the distance from the event coordinates to the image center, "x_center" is the image center coordiante and
% "k_1, k_2, ..." are the distortion parameters which are unique to every lens. This method allows to use an arbitrary amount
% of distortion parameters, where a high amount of parameters generally is used for a very precise distortion correction. 
% However for the lenses used with the eDVS128 sensor is is sufficient to use 2 or 3 parameters. These parameters can be 
% found by trial and error with the matlab script "CalibDistortionParas".
% 
% ###########################################################################################################################
%
% Input Parameters:
%
% "Events" :
%
% Matrix containing the eDVS128 event data. Every line should represent one event:   [ x-position ; y-positiotn ; ... ]
%
% "k1, k2, k3" (default value 0) :
%
% Distortion parameters. 
%
% "coord_mid" (default value 63.5) :
%
% The coordinate which represents the middle for the coordinate system of the events. For example, the eDVS128 sensor 
% produces event coordinates ranging from 0 to 127. Therefore, the middle of the event coordiante system is 127/2 which 
% yields 63.5.
%
% ###########################################################################################################################
%
% Output Parameters:
%
% The output matrix "Events" contains the corrected event coordianates and has the same data format as the input.
%
% ###########################################################################################################################
%
% author        : Martin Medler
% version       : 1.1
% date          : 14.10.2012


%% parameters

if nargin < 5
    coord_mid = 63.5;
end
if nargin < 4
    k3 = 0.0;
end
if nargin < 3
    k2 = 0.0;
end
if nargin < 2
    k1 = 0.0;
end


%% correct the distortion

R = sqrt( (Events(:,1)-coord_mid).^2 + (Events(:,2)-coord_mid).^2 );
Events(:,1) = ( Events(:,1)-coord_mid ) .* ( 1 + k1*R + k2*R.^2 + k3*R.^3 ) + coord_mid;
Events(:,2) = ( Events(:,2)-coord_mid ) .* ( 1 + k1*R + k2*R.^2 + k3*R.^3 ) + coord_mid;

end
 