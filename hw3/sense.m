% fully-sampled 8-channel k-space data
% measured coil sensitivity maps
load Data_Assignment3_Problem2

%% 2a
% Using the fully-sampled data and sensitivity maps, compute the coil-
% combined image. Display the magnitude of this image

[nx, ny, ncoils] = size(kspaceData);
imageData = zeros(size(kspaceData));

% convert from k-space to image space and multiply by scaling factor
for i = 1:ncoils
    imageData(:,:,i) = (ifftshift(ifft2(kspaceData(:,:,i)))).*coilmaps(:,:,i);
end

% combine info from the different coils and plot
im_fully_samp = abs(sum(imageData,3));
figure
imagesc(im_fully_samp);

%% 2b

% Retrospectively undersampled the k-space data using an acceleration
% factor of R=2 by setting every other phase encoding line equal to zero. 
% Assume that the phase encoding direction is oriented vertically. 
% Display the resulting magnitude image. You
% should see aliasing along the phase encoding direction

kspaceData_R2 = kspaceData;
for j = 2:2:ny
    kspaceData_R2(j,:,:) = 0;
end

[nx_R2, ny_R2, ncoils_R2] = size(kspaceData_R2);
imageData_R2 = zeros(size(kspaceData_R2));

% convert from k-space to image space and multiply by scaling factor
for i_R2 = 1:ncoils_R2
    imageData_R2(:,:,i_R2) = (ifftshift(ifft2(kspaceData_R2(:,:,i)))).*coilmaps(:,:,i);
end

% combine info from the different coils and plot
im_samp_R2 = abs(sum(imageData_R2,3));
figure
imagesc(im_samp_R2);

%% 2c

%SENSE R=2 Reconstruction: Using only the undersampled k-space data from Part b and the
%coil sensitivity maps, implement your own SENSE reconstruction and display the
%reconstructed magnitude image. In addition, compute the difference between the SENSE
%reconstruction and the fully-sampled image, and display the magnitude of the difference
%image.

%imageData_R2
recon_im = zeros(nx, ny);
R = 2;

% I = Cp
% calculate the locations for the values in the p matrix
for ii = 1:nx
    for jj = 1:(ny/R/2)
        p1_loc = jj; % check this, I think it is wrong
        p2_loc = jj + (ny/R); % check this, I think it is wrong

        % set up the I matrix (overlaping pixels)
        % size 8 x 1
        I = im_samp_R2(p1_loc, p2_loc, :);
        
        % set up the C matrix (coil sensitivities)
        C = [coilmaps(ii, p1_loc, :), coilmaps(ii, p2_loc, :)]; % need to check the dimensions of this
        
        % calculate the p matrix using the pseudoinverse
        p = pinv(C)*I;

        recon_im(ii,p1_loc) = p(1);
        recon_im(ii,p2_loc) = p(2);

    end
end

%%
%%%%% continue from here 

% for nx = 1:Nx
%     for ny = Ny/R/2 + 1 : Ny/R*3/2
%         % Calculate source pixel locations
%         ny1 = ny;
%         ny2 = ny + Ny/R;
%         if ny2 > Ny
%             ny2 = mod(ny2,Ny);
%         end
% 
%         % Recover pixel values via pseudoinverse
%         pixels_aliased = squeeze(imgs_R2(nx,ny,:)); % 8 x 1
%         smap_weights = [squeeze(smaps(nx,ny1,:)),...
%                         squeeze(smaps(nx,ny2,:))]; % 8 x 2;
%         pixels_unaliased = pinv(smap_weights)*pixels_aliased; % 2 x 1
%         
%         % Allocate recovered pixel values
%         img(nx,ny1) = pixels_unaliased(1);
%         img(nx,ny2) = pixels_unaliased(2);
%     end
% end
% 
% diff = img - img_fs;
% 
% figure;
% subplot(121); im(abs(img));
% subplot(122); im(abs(diff));
% sgtitle('2c. SENSE recon for R = 2. Left is reconstructed image. Right is difference image')

