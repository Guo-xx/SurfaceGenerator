function surface_generation(rL, N, sigma, sk, ku, Bx, By)

%% Parameters

% rL     length of surface
% N      number of nodes
% sigma  standard deviation of rough surface
% sk     skewness of rough surface
% ku     kurtosis of rough surface
% Bx     autocorrelation length in x
% By     autocorrelation length in y

%% Surface Generation

% Check Skewness and Kurtosis
if (ku - sk^2 - 1 <= 0)
    error('Invalid Skewness and Kurtosis. Skewness and Kurtosis must fulfil ku - sk^2 - 1 > 0.')
end

x = linspace(-rL/2,rL/2,N);
y = linspace(-rL/2,rL/2,N);
[X,Y] = meshgrid(x,y);

eta = randn(N,N);

% Autocorrelation Function
Rzz = sigma^2*exp(-2.3*sqrt((X/Bx).^2+(Y/By).^2));

% FFT Filter
Z_G = ifft2(fft2(eta).*sqrt(abs(fft2(Rzz))));
Z_G = Z_G./std(Z_G(:));

% Non-Gaussian Transformation
result = f_johnson_M(0,sigma,sk,ku);

Z_lin_G = reshape(Z_G,N*N,1);
Z_lin_NG = f_johnson_z2y(Z_lin_G,result.coef,result.type);
Z_NG = reshape(Z_lin_NG,N,N);

%% Write to File

fileID = fopen('surface.txt','w');
formatSpec = '%15.8f';

fprintf(fileID,formatSpec,x);
fprintf(fileID,'\r\n');
fprintf(fileID,formatSpec,y);
fprintf(fileID,'\r\n');

for i=1:N
    fprintf(fileID,formatSpec,Z_NG(i,:));
    fprintf(fileID,'\r\n');
end

fclose(fileID);

end