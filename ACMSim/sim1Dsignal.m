function x=sim1Dsignal(R,cols,rows);
% x=sim1Dsignal(R,rows,cols);
% R: 1D autocorr. function 
% Simulation by circular embedding 
% [1] percival2005, "Exact simulation of complex-valued Gaussian stationary processes via circulant embedding"
%
N = length(R); M = 2*N - 1;
Y = [conj(R) fliplr(R(2:end))].'; % Circulant embedding
%Y = [fliplr(conj(R)) R(2:end)].'; % This one works too, straight forward, needs ifft!
Sk = real(fft(Y)); % Should be real (and positive definite)...
Sk = repmat(Sk,1,rows);
Z = randn(M,rows)+i*randn(M,rows);
z = fft(Z.*sqrt(Sk/(2*M)));
%z = ifft(Z.*sqrt(Sk*M/2));
x = z(1:cols,1:rows);

if(0)
    figure(106);
    clf;
    plot(20*log10(abs(fftshift(fft(z(:,1))))));
    %Rx = 1/rows*x.'*conj(x);
end

