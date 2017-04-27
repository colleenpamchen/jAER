function [ CON, COFF ] = hots(CIN, vstart, vend )

% Parse the CIN input event into component parts, and truncate video 
    n = size(CIN,1);
    ts = CIN(:,1); 
    x0 = CIN(:,4);
    y0 = CIN(:,5);
    xy0 = CIN(:,4:5); 
    P = CIN(:,6);
% initialize constant variables:
pixels= 128; 

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

    res_x = pixels;
    res_y = pixels;

%% ON polarity events 
    Son = nan(pixels);
    ONS=[]; 
    deltaTon = zeros(pixels);
    ts_prev=0; 
    tau = 20000 ;
    x = x_incr+1;
    y = y_incr+1;
    % Build time surfaces for ON events 
    for i = vstart:length(ts_incr)/1.25
        
        ts_current = ts_incr(i) ; 
        deltaTon = deltaTon + (ts_current - ts_prev) ; 
        ts_prev = ts_current;
        deltaTon( x(i), y(i) ) = 0;
        Son = exp( -(deltaTon)/tau )  ;  
        ONS(:,:,i) = Son ;
    end
% PROTOTYPE construction using nearest neighbor and distance metrics 
Cidx=[10:100:1049]; % N=11 
npixels= pixels*pixels; 
Cn = reshape( permute( ONS(:,:,Cidx(:)), [3 1 2] ), [length(Cidx), npixels]);
CON = zeros(size(Cn));
pk = 0; 
beta = -1;

% ON polarity events loop
for i=idxstart:size(ONS,3)

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
    
%% OFF polarity events 
    Soff = nan(pixels);
    OFFS=[]; 
    deltaToff = zeros(pixels);
    ts_prev=0; 
    tau = 20000 ;
    x = x_decr+1;
    y = y_decr+1;
    
    for i = vstart:ts_decr(end) %length(ts_decr)/1.25
    ts_current = ts_decr(i) ; 
    deltaToff = deltaToff + (ts_current - ts_prev) ; 
    ts_prev = ts_current;
    deltaToff( x(i), y(i) ) = 0;
    Soff = exp( -(deltaToff)/tau )  ;  
    OFFS(:,:,i) = Soff ;
    end

% Display the surface and contounrs for ON and OFF events 
% e=2050; % event snapshot 
on = length(ts_incr);
off= length(ts_decr);

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

% OFF polarity events loop

Cidx=[10:100:1049]; % N=11 
npixels = pixels*pixels; % 16384
Cn = reshape( permute( OFFS(:,:,Cidx(:)), [3 1 2] ), [length(Cidx), npixels]);
COFF = zeros(size(Cn));
pk = 0; 
beta = -1;  

for i=idxstart:idxend
 
   
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
     
    idx(1)
% figure
% subplot(1,2,1)
% contour(OFFS(:,:,i) )
% title('test')
% subplot(1,2,2)  
% contour( squeeze(C(idx(1),:,:) ) ) 
% title('PROTOTYPE|train') 
end

figure
subplot(1,2,1)
contour(OFFS(:,:,i) )
title('test')
subplot(1,2,2)  
contour( squeeze(COFF(idx(1),:,:) ) ) 
title( sprintf('PROTOTYPE|train %d', idx(1)) )









end
