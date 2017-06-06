function savePrototype(CON,COFF,Cidx)

workingDir = './';
mkdir(workingDir,'current_images')

on_prototypes = permute(reshape(CON,[length(Cidx),128,128]), [2,3,1]);
off_prototypes = permute(reshape(COFF,[length(Cidx),128,128]), [2,3,1]);

for k=1:length(Cidx)
    newFig = figure(k);
    subplot(2,1,1)
    surfc(on_prototypes(:,:,k) ); colorbar
    title(sprintf('ON polarity Prototype %d',k))
    view(0,90)
    
    subplot(2,1,2)
    surfc(off_prototypes(:,:,k) ); colorbar
    title(sprintf('OFF polarity Prototype %d',k))
    view(0,90)
    
    filename = [sprintf('prototype_%d',k) '.png'];
    fullname = fullfile(workingDir,'current_images',filename);
    saveas(newFig, fullname,'png')  

end

imageNames = dir(fullfile(workingDir,'current_images','*.png'));
imageNames = {imageNames.name}';

outputVideo = VideoWriter(fullfile(workingDir,'Timesurface_prototypes.avi'));
outputVideo.FrameRate = 1;
open(outputVideo)

for ii = 1:length(imageNames)
   img = imread(fullfile(workingDir,'current_images',imageNames{ii}));
   writeVideo(outputVideo,img)
end

close(outputVideo)

shuttleAvi = VideoReader(fullfile(workingDir,'Timesurface_prototypes.avi'));


%%
% pCON=reshape(CON,[length(Cidx),128,128]); %for plotting/displaying purpose 
% 
% figure
% subplot(ceil(length(Cidx)/2),2,1)
% contour( squeeze(pCON(1,:,:) ) )  , colorbar
% title( sprintf('ON PROTOTYPE 1' ))
% subplot(ceil(length(Cidx)/2),2,2)
% contour( squeeze(pCON(2,:,:) ) ) , colorbar
% title( sprintf('ON PROTOTYPE 2' ))
% subplot(ceil(length(Cidx)/2),2,3)
% contour( squeeze(pCON(3,:,:) ) ) , colorbar
% title( sprintf('ON PROTOTYPE 3'))
% subplot(ceil(length(Cidx)/2),2,4)
% contour( squeeze(pCON(4,:,:) ) ) , colorbar
% title( sprintf('ON PROTOTYPE 4' ))
% subplot(ceil(length(Cidx)/2),2,5)
% contour( squeeze(pCON(5,:,:) ) ) , colorbar
% title( sprintf('ON PROTOTYPE 5' ))
% subplot(ceil(length(Cidx)/2),2,6)
% contour( squeeze(pCON(6,:,:) ) ) , colorbar
% title( sprintf('ON PROTOTYPE 6' ))
% subplot(ceil(length(Cidx)/2),2,7)
% contour( squeeze(pCON(7,:,:) ) ) , colorbar
% title( sprintf('ON PROTOTYPE 7' ))
% subplot(ceil(length(Cidx)/2),2,8)
% contour( squeeze(pCON(8,:,:) ) ) , colorbar
% title( sprintf('ON PROTOTYPE 8' ))
% subplot(ceil(length(Cidx)/2),2,9)
% contour( squeeze(pCON(9,:,:) ) ) , colorbar
% title( sprintf('ON PROTOTYPE 9' ))
% subplot(ceil(length(Cidx)/2),2,10)
% contour( squeeze(pCON(10,:,:) ) ) , colorbar
% title( sprintf('ON PROTOTYPE 10' ))
% suptitle('ON polarity Prototypes')
% 
% % View the prototype : OFF
% pCOFF=reshape(COFF,[length(Cidx),128,128]);
% 
% figure
% subplot(ceil(length(Cidx)/2),2,1)
% contour( squeeze(pCOFF(1,:,:) ) ) , colorbar
% title( sprintf('OFF PROTOTYPE 1' ))
% subplot(ceil(length(Cidx)/2),2,2)
% contour( squeeze(pCOFF(2,:,:) ) ) , colorbar
% title( sprintf('OFF PROTOTYPE 2' ))
% subplot(ceil(length(Cidx)/2),2,3)
% contour( squeeze(pCOFF(3,:,:) ) ) , colorbar
% title( sprintf('OFF PROTOTYPE 3'))
% subplot(ceil(length(Cidx)/2),2,4)
% contour( squeeze(pCOFF(4,:,:) ) ) , colorbar
% title( sprintf('OFF PROTOTYPE 4' ))
% subplot(ceil(length(Cidx)/2),2,5)
% contour( squeeze(pCOFF(5,:,:) ) ) , colorbar
% title( sprintf('OFF PROTOTYPE 5' ))
% subplot(ceil(length(Cidx)/2),2,6)
% contour( squeeze(pCOFF(6,:,:) ) ) , colorbar
% title( sprintf('OFF PROTOTYPE 6' ))
% subplot(ceil(length(Cidx)/2),2,7)
% contour( squeeze(pCOFF(7,:,:) ) ) , colorbar
% title( sprintf('OFF PROTOTYPE 7' ))
% subplot(ceil(length(Cidx)/2),2,8)
% contour( squeeze(pCOFF(8,:,:) ) ) , colorbar
% title( sprintf('OFF PROTOTYPE 8' ))
% subplot(ceil(length(Cidx)/2),2,9)
% contour( squeeze(pCOFF(9,:,:) ) ) , colorbar
% title( sprintf('OFF PROTOTYPE 9' ))
% subplot(ceil(length(Cidx)/2),2,10)
% contour( squeeze(pCOFF(10,:,:) ) ) , colorbar
% title( sprintf('OFF PROTOTYPE 10' ))
% suptitle('OFF polarity Prototypes')

%% isosurface ???
% V = ONS;
% [m,n,p]=size(V);
% [X,Y,Z] = meshgrid(1:n,1:m,1:p); 
% P = patch(isosurface(X,Y,Z,V));
% isonormals(X,Y,Z,V,P)
% set(P,'FaceColor','red','EdgeColor','blue');
% daspect([1,1,1])
% view(3); axis tight
% camlight
% lighting gouraud 


end


