%% Partial Fourier Imaging

% zero filled recon - part a

% Zero-Filled Reconstruction: Load Data_Assignment3_Problem1.mat into MATLAB. This
% contains fully-sampled k-space data from a T1-weighted brain scan. The data only has a
% single receiver coil (to keep things simple for now!). 
% Retrospectively undersample this data
% using a partial Fourier factor of 5/8. This is a bit larger than what is often used in practice –
% however, this will make it easier to see the artifacts. 
% Assume that the phase encoding
% direction is oriented vertically. Perform a zero-filled reconstruction, and display both the
% magnitude and phase of the image. Next, compute the difference between this image and
% the fully-sampled image, and display the magnitude and phase of the difference image. Note
% that you may need to adjust the windowing level to better visualize the difference image.

load Data_Assignment3_Problem1

% use the middle 125 points
[nx,ny] = size(kspaceData_SingleCoil);
num_skip = ny-(5/8)*ny;
num_skip_side = round(num_skip/2);

kspace_zero = complex(zeros(nx,ny));

kspace_zero(num_skip_side:ny-num_skip_side,:) = kspaceData_SingleCoil(num_skip_side:ny-num_skip_side,:);

im_under = ifftshift(ifft2(kspace_zero));
im_full = ifftshift(ifft2(kspaceData_SingleCoil));

figure
imagesc(abs(im_under));
title("magnitude under")

figure
imagesc(angle(im_under));
title("phase under")

diff_im = im_full - im_under;

figure
imagesc(abs(diff_im));
title("magnitude diff")

figure
imagesc(angle(diff_im));
title("phase diff")





