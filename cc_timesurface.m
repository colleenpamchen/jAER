%% colleen's working progress 
[CIN]=loadaerdat_matlab('DVS128-2017-03-21_penOrientation.aedat');
% [CIN]=loadaerdat_matlab('DVS128-2017-03-21T16-00-20-0700-0.aedat');
% [CIN]=loadaerdat_matlab('DVS128-2017-03-21-short.aedat');

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
    
%% video visualizing the events 
Events = [x0,y0,P,ts];
visualize_events(Events); 

%% TRUNCATED video 
Events = [x0(idxstart:idxend),y0(idxstart:idxend),P(idxstart:idxend),ts(idxstart:idxend)];
visualize_events(Events); 

%% Truncate video 
idxstart = 1000 ;
idxend = 127000 ;

for k=1:size(P,1)
    if P(k) == -1
       P(k) = 0;
    end
end
logical_pol = logical(P);
% ON polarity events 
x_incr = x0(logical_pol);
y_incr = y0(logical_pol);
ts_incr = ts(logical_pol);
% OFF polarity events 
x_decr = x0(~logical_pol);
y_decr = y0(~logical_pol);
ts_decr = ts(~logical_pol);

res_x = 128;
res_y = 128;

%% display the data using 3D scatterplot
% epg (events per graph) = it is associated with the scatterd mode. It
% represents the number of events displayed in each image. Default: 20000 
% epg = 20000;
epg = 50000;
% for i=ts(1): epg :ts(end)
for i=ts(5000): epg :ts(10000)
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

%% define local neighborhood U [ux, uy] 
% center of the original coordinate system (128/2) = 64
coord_mid = 64;
R = 2 ; 
ux = coord_mid-R : coord_mid+R ; 
uy = coord_mid-R : coord_mid+R ; 
% [ {62:66}  {62:66} ] 
x1=x0+1;
y1=y0+1;
xmask = x1 >= 62 & x1 <= 66 ;
ymask = y1 >= 62 & y1 <= 66 ;

% TODO: apply the spatial mask to sensor data

%% Build time surfaces for ON / OFF events 

% ON polarity events 
    Son = nan(128); 
    ONS=[]; 
    deltaTon = zeros(128);
    ts_prev=0; 
    tau = 20000 ;
    x = x_incr+1;
    y = y_incr+1;
    
    for i = idxstart: length(ts_incr)/1.25
    ts_current = ts_incr(i) ; 
    deltaTon = deltaTon + (ts_current - ts_prev) ; 
    ts_prev = ts_current;
    deltaTon( x(i), y(i) ) = 0;
    Son = exp( -(deltaTon)/tau )  ;  
    ONS(:,:,i) = Son ;
    end

% OFF polarity events 
    Soff = nan(128);
    % SOFF = zeros(128); 
    OFFS=[]; 
    deltaToff = zeros(128);
    ts_prev=0; 
    tau = 20000 ;
    x = x_decr+1;
    y = y_decr+1;
    
    for i = idxstart:length(ts_decr)/1.25
    ts_current = ts_decr(i) ; 
    deltaToff = deltaToff + (ts_current - ts_prev) ; 
    ts_prev = ts_current;
    deltaToff( x(i), y(i) ) = 0;
    Soff = exp( -(deltaToff)/tau )  ;  
    OFFS(:,:,i) = Soff ;
    end

% Display the surface and contounrs for ON and OFF events 
% e=2050; % event snapshot 
on = size(ONS,3);
off= size(OFFS,3);

figure
subplot(2,2,1)
surf(ONS(:,:,on) )
subplot(2,2,2)  
contour(ONS(:,:,on)) 
title('ON events')
subplot(2,2,3)
surf(OFFS(:,:,off) )
subplot(2,2,4)  
contour(OFFS(:,:,off)) 
title('OFF events')


% figure
% subplot(2,2,1)
% surf(ONS(:,:,on) )
% subplot(2,2,2)  
% contour(ONS(:,:,on)) 
% title('ON events')
% subplot(2,2,3)
% surf(OFFS(:,:,999) )
% subplot(2,2,4)  
% contour(OFFS(:,:,999)) 
% title('OFF events')

% TODO: visualization. play these in a movie. convert events into time.


%% PROTOTYPE construction using nearest neighbor and distance metrics 
% theta = acos(dot(x,y)/(norm(x,2)*norm(y,2)));
% this is calculating the similarity between the two surfaces in 
% terms of their degree orientation. so if they are two parallel lines, 
% they'll have a value of 1, if they are orthogonal/perpendicular, they 
% have value of 0. this is then multiplied to the prototype and subtracted 
% from the surface as the "additive difference" between the incoming event 
% surface and the prototype 'class' it belongs to. this term should be a scalar
% beta = dot(Ck,Si) / (norm(Ck)*norm(Si)) ;  
% beta = 1 - pdist([Ck;Si],'cosine') ;

Cidx=[idxstart:10:idxstart+100]; % N=10 
npixels= 128*128; 
Cn = reshape( permute( ONS(:,:,Cidx(:)), [3 1 2] ), [length(Cidx), npixels]);
CON = zeros(size(Cn));
pk = 0; 
beta = -1;

% ON polarity events loop
for i=idxstart: on
% for i=999:10:1999

    test_on_event = reshape(ONS(:,:,i),[1,npixels]); % Si 
    D_on = pdist2(Cn, test_on_event); 
        [val,idx]=sort(D_on);
    %  VECTORIZEd data 
    Ck = Cn(idx(1),:); % Ck: the closest CLUSTER center   
    Si = test_on_event; % incoming event Surface 
   
    % UPDATE RULE :
    alpha = 0.01 / (1+ (pk/20000) );
    beta = dot(Ck,Si) / (norm(Ck)*norm(Si)) ;  % beta = 1 - pdist([Ck;Si],'cosine') ;
    
    Ck = Ck + alpha.*(Si - beta.*Ck) ; 
    pk=pk+1; 

    CON(idx(1),:) = Ck; 
    CON=reshape(CON,[length(Cidx),128,128]);
     
   
%  idx(1)
% figure
% subplot(1,2,1)
% contour(ONS(:,:,i) )
% title('test')
% subplot(1,2,2)  
% contour( squeeze(C(idx(1),:,:) ) ) 
% title( sprintf('PROTOTYPE|train %d', idx(1)) )

end

figure
subplot(1,2,1)
contour(ONS(:,:,i) )
title('test')
subplot(1,2,2)  
contour( squeeze(CON(idx(1),:,:) ) ) 
title( sprintf('PROTOTYPE|train %d', idx(1)) )

%% OFF polarity events loop

Cidx=[idxstart:10:idxstart+100]; % N=10
% Cidx=[idxstart:idxstart+10]; % this still captured the line ... 
npixels = 128*128; % 16384
aa=permute( OFFS(:,:,Cidx(:)), [3 1 2] );
Cn = reshape( aa, [length(Cidx), npixels]);

COFF = zeros(size(Cn));
pk = 0; 
beta = -1;  

for i= idxstart : off
 
   
    test_off_event = reshape(OFFS(:,:,i),[1,npixels]);
    D_off = pdist2(Cn,test_off_event); 
        [val,idx]=sort(D_off);
    Ck = Cn(idx(1),:); % Ck: the closest CLUSTER center   
    Si = test_off_event; % incoming event Surface 
       
    % UPDATE RULE :
    alpha = 0.01 / (1+ (pk/20000) );
    beta = dot(Ck,Si) / (norm(Ck)*norm(Si)) ;  % beta = 1 - pdist([Ck;Si],'cosine') ;
    
    Ck = Ck + alpha.*(Si - beta.*Ck) ; 
    pk=pk+1; 

    COFF(idx(1),:) = Ck; 
    COFF=reshape(COFF,[length(Cidx),128,128]);
     
% idx(1)
% figure
% subplot(1,2,1)
% contour(OFFS(:,:,i) )
% title('test')
% subplot(1,2,2)  
% contour( squeeze(COFF(idx(1),:,:) ) ) 
% title('PROTOTYPE|train') 

end

figure
subplot(1,2,1)
contour(OFFS(:,:,i) )
title('test')
subplot(1,2,2)  
contour( squeeze(COFF(idx(1),:,:) ) ) 
title( sprintf('PROTOTYPE|train %d', idx(1)) )

%% View the prototype : ON 
figure
subplot(ceil(length(Cidx)/2),2,1)
contour( squeeze(CON(1,:,:) ) ) 
title( sprintf('ON PROTOTYPE 1' ))
subplot(ceil(length(Cidx)/2),2,2)
contour( squeeze(CON(2,:,:) ) ) 
title( sprintf('ON PROTOTYPE 2' ))
subplot(ceil(length(Cidx)/2),2,3)
contour( squeeze(CON(3,:,:) ) ) 
title( sprintf('ON PROTOTYPE 3'))
subplot(ceil(length(Cidx)/2),2,4)
contour( squeeze(CON(4,:,:) ) ) 
title( sprintf('ON PROTOTYPE 4' ))
subplot(ceil(length(Cidx)/2),2,5)
contour( squeeze(CON(5,:,:) ) ) 
title( sprintf('ON PROTOTYPE 5' ))
subplot(ceil(length(Cidx)/2),2,6)
contour( squeeze(CON(6,:,:) ) ) 
title( sprintf('ON PROTOTYPE 6' ))
subplot(ceil(length(Cidx)/2),2,7)
contour( squeeze(CON(7,:,:) ) ) 
title( sprintf('ON PROTOTYPE 7' ))
subplot(ceil(length(Cidx)/2),2,8)
contour( squeeze(CON(8,:,:) ) ) 
title( sprintf('ON PROTOTYPE 8' ))
subplot(ceil(length(Cidx)/2),2,9)
contour( squeeze(CON(9,:,:) ) ) 
title( sprintf('ON PROTOTYPE 9' ))
subplot(ceil(length(Cidx)/2),2,10)
contour( squeeze(CON(10,:,:) ) ) 
title( sprintf('ON PROTOTYPE 10' ))
% subplot(6,2,11)
% contour( squeeze(CON(11,:,:) ) ) 
% title( sprintf('PROTOTYPE|train 11' ))

%% View the prototype : OFF
figure
subplot(ceil(length(Cidx)/2),2,1)
contour( squeeze(COFF(1,:,:) ) ) 
title( sprintf('OFF PROTOTYPE 1' ))
subplot(ceil(length(Cidx)/2),2,2)
contour( squeeze(COFF(2,:,:) ) ) 
title( sprintf('OFF PROTOTYPE 2' ))
subplot(ceil(length(Cidx)/2),2,3)
contour( squeeze(COFF(3,:,:) ) ) 
title( sprintf('OFF PROTOTYPE 3'))
subplot(ceil(length(Cidx)/2),2,4)
contour( squeeze(COFF(4,:,:) ) ) 
title( sprintf('OFF PROTOTYPE 4' ))
subplot(ceil(length(Cidx)/2),2,5)
contour( squeeze(COFF(5,:,:) ) ) 
title( sprintf('OFF PROTOTYPE 5' ))
subplot(ceil(length(Cidx)/2),2,6)
contour( squeeze(COFF(6,:,:) ) ) 
title( sprintf('OFF PROTOTYPE 6' ))
subplot(ceil(length(Cidx)/2),2,7)
contour( squeeze(COFF(7,:,:) ) ) 
title( sprintf('OFF PROTOTYPE 7' ))
subplot(ceil(length(Cidx)/2),2,8)
contour( squeeze(COFF(8,:,:) ) ) 
title( sprintf('OFF PROTOTYPE 8' ))
subplot(ceil(length(Cidx)/2),2,9)
contour( squeeze(COFF(9,:,:) ) ) 
title( sprintf('OFF PROTOTYPE 9' ))
subplot(ceil(length(Cidx)/2),2,10)
contour( squeeze(COFF(10,:,:) ) ) 
title( sprintf('OFF PROTOTYPE 10' ))
% subplot(6,2,11)
% contour( squeeze(C(11,:,:) ) ) 
% title( sprintf('PROTOTYPE|train 11' ))

%% INITIAL prototypes : 
figure; subplot(ceil(length(Cidx)/2),2,1); contour( ONS(:,:,Cidx(1)))
subplot(ceil(length(Cidx)/2),2,2); contour( ONS(:,:,Cidx(2)))
subplot(ceil(length(Cidx)/2),2,3); contour( ONS(:,:,Cidx(3)))
subplot(ceil(length(Cidx)/2),2,4); contour( ONS(:,:,Cidx(4)))
subplot(ceil(length(Cidx)/2),2,5); contour( ONS(:,:,Cidx(5)))
subplot(ceil(length(Cidx)/2),2,6); contour( ONS(:,:,Cidx(6)))
subplot(ceil(length(Cidx)/2),2,7); contour( ONS(:,:,Cidx(7)))
subplot(ceil(length(Cidx)/2),2,8); contour( ONS(:,:,Cidx(8)))
subplot(ceil(length(Cidx)/2),2,9); contour( ONS(:,:,Cidx(9)))
subplot(ceil(length(Cidx)/2),2,10); contour( ONS(:,:,Cidx(10)))

 
figure; subplot(ceil(length(Cidx)/2),2,1); contour( OFFS(:,:,Cidx(1)))
subplot(ceil(length(Cidx)/2),2,2); contour( OFFS(:,:,Cidx(2)))
subplot(ceil(length(Cidx)/2),2,3); contour( OFFS(:,:,Cidx(3)))
subplot(ceil(length(Cidx)/2),2,4); contour( OFFS(:,:,Cidx(4)))
subplot(ceil(length(Cidx)/2),2,5); contour( OFFS(:,:,Cidx(5)))
subplot(ceil(length(Cidx)/2),2,6); contour( OFFS(:,:,Cidx(6)))
subplot(ceil(length(Cidx)/2),2,7); contour( OFFS(:,:,Cidx(7)))
subplot(ceil(length(Cidx)/2),2,8); contour( OFFS(:,:,Cidx(8)))
subplot(ceil(length(Cidx)/2),2,9); contour( OFFS(:,:,Cidx(9)))
subplot(ceil(length(Cidx)/2),2,10); contour( OFFS(:,:,Cidx(10)))

%% ... SCRATCH code ....
% TODO: change it so that it takes incoming data one by one 

