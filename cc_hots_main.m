%% PRE-RECORDED DVS DATASET from Jorg Conradt 
% data already in matrix format: 
    % Column 1: x coordinate (from 0 to 127) 
    % Column 2: y coordinate (from 0 to 127)
    % Column 3: event polarity [0 off | 1 on]
    % Column 4: timestamps with 1us time tick (microseconds?miliseconds? how to convert to seconds?)
    % Events = [ x-coordinate , y-coordinate , polarity , timestamp (microseconds) ]
Events1_v = load('sample_eDVS_data/pen_vertical.dvs');
Events1_h = load('sample_eDVS_data/pen_horizontal.dvs');
Events2 = load('sample_eDVS_data/spinner.dvs');
Events3_v = load('sample_eDVS_data/hand_vertical.dvs');
Events3_h = load('sample_eDVS_data/hand_horizontal.dvs');

Events= Events2; 

%% video visualizing the events 
visualize_events(Events); 

%% Batch incoming events: 
num_events = size(Events,1) ; 
batchsize = 1000;
end_video = floor(num_events/batchsize) ;

%% INITIAL prototypes 
Cidx=[1:10]; % Cn=10
pixels= 128; 
npixels= pixels*pixels;  % 16384
Cn_on = zeros(length(Cidx),npixels);
Cn_off = zeros(length(Cidx),npixels);

%%

tau = 20000 ; % PARAMETER taken from Lagorce et. al. 

eidx=[1:batchsize];
for i = 1:end_video                 %end_video
    data = Events(eidx,:); 
%     visualize_events(data); 
    [CON, COFF]= hots(data, Cn_on, Cn_off, tau, Cidx) ; 
            Cn_on = CON;
            Cn_off = COFF;

    eidx = eidx+batchsize; 
%     pause 

end


%% View the accum. prototype : ON 
pCON=reshape(CON,[length(Cidx),128,128]); 
for p = 1:length(Cidx)
    figure; contour(squeeze(pCON(p,:,:))) ; 
    title(sprintf('On prototype %d',p)),colorbar 
end
pCOFF=reshape(COFF,[length(Cidx),128,128]); 
for p = 1:length(Cidx)
    figure; contour(squeeze(pCOFF(p,:,:))) ; 
    title(sprintf('Off prototype %d',p)), colorbar 
end

%%
pCON=reshape(CON,[length(Cidx),128,128]); %for plotting/displaying purpose 

figure
subplot(ceil(length(Cidx)/2),2,1)
contour( squeeze(pCON(1,:,:) ) )  , colorbar
title( sprintf('ON PROTOTYPE 1' ))
subplot(ceil(length(Cidx)/2),2,2)
contour( squeeze(pCON(2,:,:) ) ) , colorbar
title( sprintf('ON PROTOTYPE 2' ))
subplot(ceil(length(Cidx)/2),2,3)
contour( squeeze(pCON(3,:,:) ) ) , colorbar
title( sprintf('ON PROTOTYPE 3'))
subplot(ceil(length(Cidx)/2),2,4)
contour( squeeze(pCON(4,:,:) ) ) , colorbar
title( sprintf('ON PROTOTYPE 4' ))
subplot(ceil(length(Cidx)/2),2,5)
contour( squeeze(pCON(5,:,:) ) ) , colorbar
title( sprintf('ON PROTOTYPE 5' ))
subplot(ceil(length(Cidx)/2),2,6)
contour( squeeze(pCON(6,:,:) ) ) , colorbar
title( sprintf('ON PROTOTYPE 6' ))
subplot(ceil(length(Cidx)/2),2,7)
contour( squeeze(pCON(7,:,:) ) ) , colorbar
title( sprintf('ON PROTOTYPE 7' ))
subplot(ceil(length(Cidx)/2),2,8)
contour( squeeze(pCON(8,:,:) ) ) , colorbar
title( sprintf('ON PROTOTYPE 8' ))
subplot(ceil(length(Cidx)/2),2,9)
contour( squeeze(pCON(9,:,:) ) ) , colorbar
title( sprintf('ON PROTOTYPE 9' ))
subplot(ceil(length(Cidx)/2),2,10)
contour( squeeze(pCON(10,:,:) ) ) , colorbar
title( sprintf('ON PROTOTYPE 10' ))
suptitle('ON polarity Prototypes')


%% View the prototype : OFF
pCOFF=reshape(COFF,[length(Cidx),128,128]);

figure
subplot(ceil(length(Cidx)/2),2,1)
contour( squeeze(pCOFF(1,:,:) ) ) , colorbar
title( sprintf('OFF PROTOTYPE 1' ))
subplot(ceil(length(Cidx)/2),2,2)
contour( squeeze(pCOFF(2,:,:) ) ) , colorbar
title( sprintf('OFF PROTOTYPE 2' ))
subplot(ceil(length(Cidx)/2),2,3)
contour( squeeze(pCOFF(3,:,:) ) ) , colorbar
title( sprintf('OFF PROTOTYPE 3'))
subplot(ceil(length(Cidx)/2),2,4)
contour( squeeze(pCOFF(4,:,:) ) ) , colorbar
title( sprintf('OFF PROTOTYPE 4' ))
subplot(ceil(length(Cidx)/2),2,5)
contour( squeeze(pCOFF(5,:,:) ) ) , colorbar
title( sprintf('OFF PROTOTYPE 5' ))
subplot(ceil(length(Cidx)/2),2,6)
contour( squeeze(pCOFF(6,:,:) ) ) , colorbar
title( sprintf('OFF PROTOTYPE 6' ))
subplot(ceil(length(Cidx)/2),2,7)
contour( squeeze(pCOFF(7,:,:) ) ) , colorbar
title( sprintf('OFF PROTOTYPE 7' ))
subplot(ceil(length(Cidx)/2),2,8)
contour( squeeze(pCOFF(8,:,:) ) ) , colorbar
title( sprintf('OFF PROTOTYPE 8' ))
subplot(ceil(length(Cidx)/2),2,9)
contour( squeeze(pCOFF(9,:,:) ) ) , colorbar
title( sprintf('OFF PROTOTYPE 9' ))
subplot(ceil(length(Cidx)/2),2,10)
contour( squeeze(pCOFF(10,:,:) ) ) , colorbar
title( sprintf('OFF PROTOTYPE 10' ))
suptitle('OFF polarity Prototypes')


%%  NMF 

bf=2;
ck = squeeze(COFF(10,:,:));
[W,h]=nnmf(ck,bf);
imagesc(W)


opt = statset('MaxIter',5,'Display','iter');
[W0,H0] = nnmf(X,5,'replicates',10,...
                   'options',opt,...
                   'algorithm','mult');

%% image J. fiji. volumetric viewer .....   
% isosurface()  
% 
help isosurface 







