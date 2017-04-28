function [ CON, COFF ] = hots(events)

% Parse the CIN input event into component parts, and truncate video 
    N = size(events,1);
    x0 = events(:,1); 
    y0 = events(:,2);
    P = events(:,3);
    ts = events(:,4); 
   
% initialize constant variables:
pixels= 128; 

    logical_pol = logical(P);
    % ON polarity events 
    x_incr = x0(logical_pol);
    y_incr = y0(logical_pol);
    ts_incr = ts(logical_pol);
    % OFF polarity events 
    x_decr = x0(~logical_pol);
    y_decr = y0(~logical_pol);
    ts_decr = ts(~logical_pol);

    res_x = pixels;
    res_y = pixels;

% ON polarity events 
    Son = nan(pixels);
    ONS=[]; 
    deltaTon = zeros(pixels);
    ts_prev=0; 
    tau = 20000 ; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% tau is a parameter!!! 
    x = x_incr+1;
    y = y_incr+1;
    % Build time surfaces for ON events 
    for i = 1:length(ts_incr)
        
        ts_current = ts_incr(i) ; 
        deltaTon = deltaTon + (ts_current - ts_prev) ; 
        ts_prev = ts_current;
        deltaTon( x(i), y(i) ) = 0;
        Son = exp( -(deltaTon)/tau )  ;  
        ONS(:,:,i) = Son ;
    end
% PROTOTYPE construction using nearest neighbor and distance metrics 
Cidx=[1:10]; % Cn=10
npixels= pixels*pixels; 
Cn = reshape( permute( ONS(:,:,Cidx(:)), [3 1 2] ), [length(Cidx), npixels]);
CON = zeros(size(Cn));
pk = 0; 
beta = -1;

% ON polarity events loop
for i=1:size(ONS,3)

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
     
   
                        % idx(1)
                        % figure
                        % subplot(1,2,1)
                        % contour(ONS(:,:,i) )
                        % title('test')
                        % subplot(1,2,2)  
                        % contour( squeeze(C(idx(1),:,:) ) ) 
                        % title( sprintf('PROTOTYPE|train %d', idx(1)) )

end

                        % figure
                        % subplot(1,2,1)
                        % contour(ONS(:,:,i) )
                        % title('test')
                        % subplot(1,2,2)  
                        % contour( squeeze(CON(idx(1),:,:) ) ) 
                        % title( sprintf('PROTOTYPE|train %d', idx(1)) )
    
% OFF polarity events 
    Soff = nan(pixels);
    OFFS=[]; 
    deltaToff = zeros(pixels);
    ts_prev=0; 
%     tau = 20000 ;
    x = x_decr+1;
    y = y_decr+1;
    
    for i = 1:length(ts_decr)
    ts_current = ts_decr(i) ; 
    deltaToff = deltaToff + (ts_current - ts_prev) ; 
    ts_prev = ts_current;
    deltaToff( x(i), y(i) ) = 0;
    Soff = exp( -(deltaToff)/tau )  ;  
    OFFS(:,:,i) = Soff ;
    end
    
% OFF polarity events loop

Cidx=[1:10]; % Cn=10
npixels = pixels*pixels; % 16384
Cn = reshape( permute( OFFS(:,:,Cidx(:)), [3 1 2] ), [length(Cidx), npixels]);
COFF = zeros(size(Cn));
pk = 0; 
beta = -1;  

for i=1:size(OFFS,3)
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

                                % figure
                                % subplot(1,2,1)
                                % contour(OFFS(:,:,i) )
                                % title('test')
                                % subplot(1,2,2)  
                                % contour( squeeze(COFF(idx(1),:,:) ) ) 
                                % title( sprintf('PROTOTYPE|train %d', idx(1)) )




% Display the surface and contounrs for ON and OFF events 
% e=2050; % event snapshot 
on = length(ts_incr);
off= length(ts_decr);

% figure
% subplot(2,2,1)
% surf(ONS(:,:,on) )
% subplot(2,2,2)  
% contour(ONS(:,:,on)) 
% title('ON events')
% subplot(2,2,3)
% surf(OFFS(:,:,off) )
% subplot(2,2,4)  
% contour(OFFS(:,:,off)) 
% title('OFF events')


%% INITIAL prototypes : 
figure; 
subplot(ceil(length(Cidx)/2),2,1); contour( ONS(:,:,Cidx(1)))
subplot(ceil(length(Cidx)/2),2,2); contour( ONS(:,:,Cidx(2)))
subplot(ceil(length(Cidx)/2),2,3); contour( ONS(:,:,Cidx(3)))
subplot(ceil(length(Cidx)/2),2,4); contour( ONS(:,:,Cidx(4)))
subplot(ceil(length(Cidx)/2),2,5); contour( ONS(:,:,Cidx(5)))
subplot(ceil(length(Cidx)/2),2,6); contour( ONS(:,:,Cidx(6)))
subplot(ceil(length(Cidx)/2),2,7); contour( ONS(:,:,Cidx(7)))
subplot(ceil(length(Cidx)/2),2,8); contour( ONS(:,:,Cidx(8)))
subplot(ceil(length(Cidx)/2),2,9); contour( ONS(:,:,Cidx(9)))
subplot(ceil(length(Cidx)/2),2,10); contour( ONS(:,:,Cidx(10)))
suptitle('ON polarity initial prototypes')
 
figure; 
subplot(ceil(length(Cidx)/2),2,1); contour( OFFS(:,:,Cidx(1)))
subplot(ceil(length(Cidx)/2),2,2); contour( OFFS(:,:,Cidx(2)))
subplot(ceil(length(Cidx)/2),2,3); contour( OFFS(:,:,Cidx(3)))
subplot(ceil(length(Cidx)/2),2,4); contour( OFFS(:,:,Cidx(4)))
subplot(ceil(length(Cidx)/2),2,5); contour( OFFS(:,:,Cidx(5)))
subplot(ceil(length(Cidx)/2),2,6); contour( OFFS(:,:,Cidx(6)))
subplot(ceil(length(Cidx)/2),2,7); contour( OFFS(:,:,Cidx(7)))
subplot(ceil(length(Cidx)/2),2,8); contour( OFFS(:,:,Cidx(8)))
subplot(ceil(length(Cidx)/2),2,9); contour( OFFS(:,:,Cidx(9)))
subplot(ceil(length(Cidx)/2),2,10); contour( OFFS(:,:,Cidx(10)))
suptitle('OFF polarity initial prototypes')

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
suptitle('ON polarity Prototypes')

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
suptitle('OFF polarity Prototypes')


end % end of function 
