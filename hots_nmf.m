%% NMF of time-surfaces 
% [ ] TIME STEPS ... it's ok if it's sparse. loop through time steps
% instead of loop through events, that way you can reconstruct the SEQUENCE
% of input events
% [ ] SPATIAL RECEPTIVE FIELDS. implement the spatial neighborhood. 


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
batchsize = 100;
end_video = floor(num_events/batchsize) ;

%% INITIAL prototypes 
NUM_prototypes= 20;
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
%     [CON ,COFF, ONS, OFFS ]= hots(data, Cn_on, Cn_off, tau, Cidx) ; 
[CON ,COFF, ONS, OFFS ]= test(data, Cn_on, Cn_off, tau, Cidx) ; 
            Cn_on = CON;
            Cn_off = COFF;
    eidx = eidx+batchsize; 
end
save CON.mat
save COFF.mat
save ONS.mat
save OFFS.mat

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

bf=20;
% ck = squeeze(COFF(10,:,:));
[Woff,Hoff,Doff]=nnmf(COFF(:,:),bf);
figure; imagesc(Woff)
hoff_disp = reshape(Hoff,[bf,128,128]);
for i = 1:bf
figure; contour(squeeze(hoff_disp(i,:,:)))
title('NMF basis vectors. Parts based representation')
end
% now RECONSTRUCT input .......... how???? 
reoff = Woff*Hoff ; % COFF
reconoff=reshape(reoff,[length(Cidx),128,128]); 
for i = 1:NUM_prototypes
figure; contour(squeeze(reconoff(i,:,:)))
end
for ii =21:50
    figure;
    subplot(1,2,1)
    contour(squeeze(pCOFF(ii,:,:)))
    title('Original OFF_prototype')
    subplot(1,2,2)
    contour(squeeze(reconoff(ii,:,:)))
    title('Reconstructed')
end


% ON 
[Won,Hon,Don]=nnmf(CON(:,:),bf);
figure; imagesc(Won)
hon_disp = reshape(Hon,[bf,128,128]);
for i = 1:bf
figure; contour(squeeze(hon_disp(i,:,:)))
title('basis vectors. parts based representation')
end
reon = Won*Hon ; % COFF
reconon=reshape(reon,[length(Cidx),128,128]); 
for i = 1:NUM_prototypes
figure; contour(squeeze(reconon(i,:,:)))
end
for ii = 131:150
    figure;
    subplot(1,2,1)
    contour(squeeze(pCON(ii,:,:)))
    title('Original ON_prototype')
    subplot(1,2,2)
    contour(squeeze(reconon(ii,:,:)))
    title('Reconstructed')
end




% opt = statset('MaxIter',1000,'Display','iter');
% [W0,H0] = nnmf(X,5,'replicates',10,...
%                    'options',opt,...
%                    'algorithm','mult');


