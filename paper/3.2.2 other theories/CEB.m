function [h,d,pC] = CEB(K,sigma_y,E1,E2,nu1,nu2,eta,beta,sigma,sk,ku,m2,m4)

%   Input:

%   K        maximum contact pressure factor
%   sigma_y  yield strength
%   E1, E2   Young's modulus
%   nu1,nu2  Poisson's ratio
%   eta      equivalent areal density of asperities
%   beta     equivalent asperity radius of curvature
%   sigma    equivalent standard deviation of surface heights
%   sk       equivalent skewness of surface heights
%   ku       equivalent kurtosis of surface heights
%   m2, m4   equivalent spectral moments

%   Output:

%   h        separation based on surface heights
%   d        separation based on asperity heights
%   pC       asperity contact pressure

%% Calculation

m0     = sigma^2;
alpha  = m0 * m4 / m2^2;
yS     = 4 * sqrt(m0/pi/alpha) / sigma;
sigmaA = sqrt(sigma^2 - 3.717E-4 / (eta^2 * beta^2));

% error check
if real(sigmaA) <= 0
    error('Invalid parameters.')
end

result = f_johnson_M(0,1,sk,ku);
verteilung_surf = @(z) Distr(z,result);
verteilung_asp = @(z) sigma/sigmaA * verteilung_surf(z * sigma/sigmaA);

E = ((1 - nu1^2)/E1 + (1 - nu2^2)/E2)^(-1);

H = 2.8*sigma_y;

omega1 = (0.5*pi*K*H/E)^2 * beta / sigma;

h = 0:0.01:4;
d = h - yS;

pC = zeros(1,length(h));

for i = 1:length(h)
    
    w_e = @(z) (z-d(i)).^1.5       .* verteilung_asp(z);
    w_p = @(z) (2*(z-d(i))-omega1) .* verteilung_asp(z);
    
    pC(i) = eta * beta * sigma * E * (4/3 * sqrt(sigma/beta) * integral(w_e,d(i),d(i)+omega1) +...
        pi * K*H/E * integral(w_p,d(i)+omega1,Inf));

end

end


function y = Distr(z,result)
y = f_johnson_pdf(z',result.coef,result.type);
y(isnan(y)) = 0;
y = y';
end