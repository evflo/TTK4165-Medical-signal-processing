function image = LoadImage(fileName)


image = imread(fileName);
image = double(image);
image = image(:,:,1);