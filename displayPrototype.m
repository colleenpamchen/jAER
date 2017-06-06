function displayPrototype(CON,COFF)

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

% View the prototype : OFF
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


end


