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

Events= Events3_h; 

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

%% MAIN loop that takes incoming batched events (1000) and 
% 

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


%%  NMF 

bf=3;
% ck = squeeze(COFF(10,:,:));
[W,h]=nnmf(COFF(:,:,:),bf);
figure; imagesc(W)


h_disp = reshape(h,[bf,128,128]);

for i = 1:bf
figure; contour(squeeze(h_disp(i,:,:)))
end


opt = statset('MaxIter',5,'Display','iter');
[W0,H0] = nnmf(X,5,'replicates',10,...
                   'options',opt,...
                   'algorithm','mult');

%% image J. fiji. volumetric viewer .....   
% isosurface()  
% 
help isosurface 







