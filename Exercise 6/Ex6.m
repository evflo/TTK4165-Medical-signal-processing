function Ex6(part1,part2)


if part1
    fileNames = {'OneLineHorizontal.bmp','OneLineHorizontal2.bmp','OneLineVertical.bmp',...
        'OneLineDiagonal.bmp','OneLineVertical2.bmp'};

    for i = 1:length(fileNames)
        Image = LoadImage(fileNames{i});
        titleName = ['Image plot of ' fileNames{i}];
        figure(i),subplot(1,2,1),imagesc(Image),title(titleName),colormap(gray);
        axis image

        imagefft = fft2(Image);

        imagefftShift = fftshift(abs(imagefft));

        imagefftShiftNorm = abs(imagefftShift)/max(max(imagefftShift));

        dynRange = 60;
        nGray = 255;

        imagePowerSpectrum = nGray*(1+20*log10(imagefftShiftNorm)/dynRange);
        titlePowerName = ['fft power spectrum of ' fileNames{i}];
        figure(i),subplot(1,2,2),image(imagePowerSpectrum),title(titlePowerName);
        colormap(gray(nGray));
        axis image
    end
end

if part2
   lenaImage = imread('bilde.jpg');
   lenaImage = double(lenaImage);
   % Get number of rows and columns in image.
    sizeImage = size(lenaImage);

    propX = 0.30;
    propY = 0.30;
    
    Xaxis = (0:sizeImage(2)-1)/sizeImage(2)*propX*100; % Creates X-axis from 0 to 30 cm.
    Yaxis = (0:sizeImage(1)-1)/sizeImage(1)*propY*100; % Creates X-axis from 0 to 30 cm.
    Yaxis = fliplr(Yaxis); % flip to get zero in lower left corner.

    deltaX = propX/sizeImage(2);
    deltaY = propY/sizeImage(1);



    figure, image(Xaxis,Yaxis,lenaImage); axis image; axis xy;
    colormap(gray(255));

    xlabel('cm');
    ylabel('cm');
   figure,subplot(1,2,1),imagesc(lenaImage),colormap(gray);
   axis image;
   lenaFFT = fft2(lenaImage); 
   lenaFFTShift = fftshift(abs(lenaFFT));
   lenaFFTShiftNorm = lenaFFTShift/max(max(lenaFFTShift));
   dynRange = 60;
   nGray = 255;
   lenaPowerSpectrum = nGray*(1+20*log10(lenaFFTShiftNorm)/dynRange);
   subplot(1,2,2),image(lenaPowerSpectrum);
   axis image;
   keyboard
end