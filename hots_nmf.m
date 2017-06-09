%% NMF of time-surfaces 

% load .dvs files
Events1_v = load('sample_eDVS_data/pen_vertical.dvs');
Events1_h = load('sample_eDVS_data/pen_horizontal.dvs');
Events2 = load('sample_eDVS_data/spinner.dvs');
Events3_v = load('sample_eDVS_data/hand_vertical.dvs');
Events3_h = load('sample_eDVS_data/hand_horizontal.dvs');
Events= Events3_h; 

% load .mat files 
%  load('matfiles_prototypes/hand_h_CON.mat')
%  load('matfiles_prototypes/hand_h_COFF.mat')


%% video visualizing the events 
% visualize_events(Events); 

%% Batch incoming events: 
num_events = size(Events,1) ; 
batchsize = 1000;
end_video = floor(num_events/batchsize) ;

%% INITIAL prototypes 
NUM_prototypes= 10;
Cidx=[1:NUM_prototypes]; % Cn=10
pixels= 128; 
npixels= pixels*pixels;  % 16384
Cn_on = zeros(length(Cidx),npixels);
Cn_off = zeros(length(Cidx),npixels);

%% MAIN loop that takes incoming batched events (1000) and 
tau = 20000 ; % PARAMETER taken from Lagorce et. al. 
eidx=[1:batchsize];
for i = 1:end_video              
    data = Events(eidx,:); 
    [CON ,COFF, ONS, OFFS ]= hots(data, Cn_on, Cn_off, tau, Cidx) ; 
            Cn_on = CON;
            Cn_off = COFF;
    eidx = eidx+batchsize; 
end

%% SAVE the prototypes as a video .avi file 
savePrototype(CON,COFF,Cidx)

%% quick View the accum. prototype : ON 
pCON=reshape(CON,[length(Cidx),128,128]); 
pCOFF=reshape(COFF,[length(Cidx),128,128]); 
for p = 1:length(Cidx)
    figure; 
    subplot(2,1,1)
    contour(squeeze(pCON(p,:,:))) ; 
    title(sprintf('On prototype %d',p)),colorbar 
    subplot(2,1,2)
    contour(squeeze(pCOFF(p,:,:))) ; 
    title(sprintf('Off prototype %d',p)), colorbar 
end


%%  NMF 

bf=5;
% ck = squeeze(COFF(10,:,:));
[W,H,D]=nnmf(COFF(:,:),bf);
figure; imagesc(W)
h_disp = reshape(H,[bf,128,128]);
for i = 1:bf
figure; contour(squeeze(h_disp(i,:,:)))
end


[Won,Hon,Don]=nnmf(CON(:,:),bf);
figure; imagesc(Won)
hon_disp = reshape(Hon,[bf,128,128]);
for i = 1:bf
figure; contour(squeeze(hon_disp(i,:,:)))
end

% now RECONSTRUCT input .......... how???? 


re = W*H ; 




% opt = statset('MaxIter',1000,'Display','iter');
% [W0,H0] = nnmf(X,5,'replicates',10,...
%                    'options',opt,...
%                    'algorithm','mult');









