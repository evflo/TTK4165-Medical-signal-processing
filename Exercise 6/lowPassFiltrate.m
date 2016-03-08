function filter = lowPassFiltrate(image)
filter = zeros(size(image,1),size(image,2));
for i = 1:size(image,1)
    for j = 1:size(image,2)
        neighbouringElements = [];
        if i>1 && i< size(image,1) && j> 1 && j<size(image,1)
            neighbouringElements = image(i-1:i+1,j-1:j+1);
        elseif i == 1 && j == 1
            neighbouringElements = image(i:i+1,j:j+1);
        elseif i==1 && j == size(image,1)
            neighbouringElements = image(i:i+1,j-1:j);
        elseif i == 1
            neighbouringElements = image(i:i+1,j-1:j+1);
        elseif j == 1 && i == size(image,1)
            neighbouringElements = image(i-1:i,j:j+1);
        elseif j == 1
            neighbouringElements = image(i-1:i+1,j:j+1);
        elseif i == size(image,1) && j == size(image,2)
            neighbouringElements = image(i-1:i,j-1:j);
        end
        neighbouringElements(i,j) = 0;        
        filter(i,j) = round(1/(size(neighbouringElements,1)*size(neighbouringElements,2)-1)*...
            sum(sum(neighbouringElements)));
    end
end
