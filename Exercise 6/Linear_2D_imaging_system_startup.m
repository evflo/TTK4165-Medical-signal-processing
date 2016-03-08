%% Exercise 5 Startup aid
%  This script is intended used as a startup aid when solving the
%  2D linear imaging systems exercise in TTK4165 Signal Processing in 
%  Medical imaging
%  Original version: ??, 20??
%  Last revision: Sten Roar Snare, Jan, 2010

%% -----------------------Part 1-------------------------------------------

%Draw lines in paint, and save them as monochrome bitmap files.

%Loads the file "vertical_line1.bmp" into Matlab. 'line1' is a 2D matrix:
line1 = double(imread('vertical_line1.bmp')); 

%Displays the image in a new figure:
figure;
imagesc(line1); %Imagesc means that the data is scaled to use the full colormap

%Sets the current figure's colormap to grayscale, uses 256 gray levels:
colormap(gray(256));

%Sets the aspect ration:
axis image

%2D Fourier transform of the lines:
Line1FT = fft2(line1);

%Use fftshift to get the 0 frequency in the middle of the figure:
Line1FT_shifted= fftshift(abs(Line1FT));

%Normalize the data to the max value:
Line1FT_norm = Line1FT_shifted/max(max(Line1FT_shifted)); %Remember: This is a 2D matrix

%Sets dynamic range - you can try to change this to see the difference:
dynRange = 60; %Large dynamic range -> more details visible in the image

%Number of gray levels:
NGray = 255;

%Calculates the image frequency specter logarithmically
imagePowerSpecter = NGray*(1+20*log10(Line1FT_norm)/dynRange);

%Displays the image in a new figure:
figure, image(imagePowerSpecter);

%Sets the current figure's colormap to grayscale, uses NGrat gray levels:
colormap(gray(NGray));

%Displays a colorbar on the side of the image if you want that
colorbar

%% --------------------Part 2------------------------------------------------

%--------------1: Load and display the image------------------------------
dynRangeImage = 100;
dynRangeMask = 50;
dynRangeRes = 100;

%Loads the image into matlab. 
bilde = double(imread('bilde.jpg'))+1;

% Get number of rows and columns in image.
sizeImage = size(bilde);

% Set the image dimensions to be 30x30cm
propX = 0.30;
propY = 0.30;
    
% Create X- and Y-axis
Xaxis = (0:sizeImage(2)-1)/sizeImage(2)*propX*100; % Creates X-axis from 0 to 30 cm.
Yaxis = (0:sizeImage(1)-1)/sizeImage(1)*propY*100; % Creates X-axis from 0 to 30 cm.
Yaxis = fliplr(Yaxis); % flip to get zero in lower left corner.

% Find metres per pixel.
deltaX = propX/sizeImage(2);
deltaY = propY/sizeImage(1);

% Calculate spatial samplingrate to make spatial frequency axis.
spatialFsX = 1/deltaX;
spatialFsY = 1/deltaY;

% Define the frequency axis from -Fs/2 to Fs/2 - Try to understand how this
% works!!
frekAxisX = (0:sizeImage(2)-1)/sizeImage(2)*spatialFsX - spatialFsX/2;
frekAxisY = (0:sizeImage(1)-1)/sizeImage(1)*spatialFsY - spatialFsY/2;

% Subplot the image. Set axis to 'image' and 'xy' to get proper scaling in normal carthesian system.
figure, image(Xaxis,Yaxis,bilde); axis image; axis xy;
colormap(gray(NGray));

% write out axis labels.
xlabel('cm');
ylabel('cm');


%----------------2: Calculate and plot the 2D FT of the image--------------

%Here you can look at the code from part 1 and try to use some of it
%again. First you use fft2 to FT the image, then you have to shift,
%normalize and scale it. Remember to use spatial frequencies on the axes.
%(frekAxisX and frekAxisY). 


%--------------3: Create a filter------------------------------------------

%Create a low pass filter - a matrix that contains the weighting of the
%pixel values of the pixels around the centre pixel. 


%--------------4: Calculate and plot the 2D FT of the filter matrix--------

%Again - look at the code from part 1. Important: When you calculate the
%fft of the filter, write fft2(filteretDitt, sizeImage(1), sizeImage(2)) to
%zero pad the FT to the same size as the image. 


%--------------5: Filter the image by using the filter from point 3--------

%This can be done in two ways; time domain or frequency domain. In time
%domain, use the function conv2(A, B,'same'); 'same' -
%returns the central part of the convolution that is the same size as A.
%You can also use the function filter2. In frequency domain, multiply the
%FT and use the inverse transform ifft. !!!IMPORTANT: Filter the original FT
%of the image, not the scaled, shifted and normalized one!!!


%---------------6: Calculate and plot the 2D FT of the filtered image------

%Use different sizes of this matrix and comment the resulting images. 