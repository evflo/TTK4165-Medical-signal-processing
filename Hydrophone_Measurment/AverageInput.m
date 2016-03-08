function AveragedInput = AverageInput(Input,Step)

AveragedInputLength = floor(length(Input)/Step);
AveragedInput = zeros(1,AveragedInputLength);
InputPosition = 1;
for i = 1:AveragedInputLength
    AveragedInput(i) = mean(Input(InputPosition:InputPosition+9));
    InputPosition = InputPosition+10;
end
