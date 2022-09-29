function [h,d,pC] = KE(sigma_y,E1,E2,nu1,nu2,eta,beta,sigma,sk,ku)

%   Input:

%   sigma_y  yield strength
%   E1, E2   Young's modulus
%   nu1,nu2  Poisson's ratio
%   eta      equivalent areal density of asperities
%   beta     equivalent asperity radius of curvature
%   sigma    equivalent standard deviation of surface heights
%   sk       equivalent skewness of surface heights
%   ku       equivalent kurtosis of surface heights

%   Output:

%   h        separation based on surface heights
%   d        separation based on asperity heights
%   pC       asperity contact pressure

%% Calculation

yS     = 1 / (sqrt(48)*pi*eta*beta*sigma);
sigmaA = sqrt(sigma^2 - 3.717E-4 / (eta^2 * beta^2));

% error check
if real(sigmaA) <= 0
    error('Invalid parameters.')
end

result = f_johnson_M(0,1,sk,ku);
verteilung_surf = @(z) Distr(z,result);
verteilung_asp = @(z) sigma/sigmaA * verteilung_surf(z * sigma/sigmaA);

E = ((1 - nu1^2)/E1 + (1 - nu2^2)/E2)^(-1);

K_KE = 0.454 + 0.41 * nu1;
H = 2.8 * sigma_y;

omega1 = (0.5*pi*K_KE*H/E)^2 * beta / sigma;

h = 0:0.01:4;
d = h - yS;

pC = zeros(1,length(h));

for i = 1:length(h)
    
    w_e   = @(z)           ((z-d(i))./omega1).^1.5   .* verteilung_asp(z);
    w_ep1 = @(z) 1.03   .* ((z-d(i))./omega1).^1.425 .* verteilung_asp(z);
    w_ep2 = @(z) 1.4    .* ((z-d(i))./omega1).^1.263 .* verteilung_asp(z);
    w_p   = @(z) 3/K_KE .* ((z-d(i))./omega1)        .* verteilung_asp(z);
    
    pC(i) = 2/3 * H * pi * eta * beta * sigma * K_KE * omega1 *...
        (integral(w_e,d(i),d(i)+omega1) +...
        integral(w_ep1,d(i)+omega1,d(i)+6*omega1) +...
        integral(w_ep2,d(i)+6*omega1,d(i)+110*omega1) +...
        integral(w_p,d(i)+110*omega1,Inf));

end

end


function y = Distr(z,result)
y = f_johnson_pdf(z',result.coef,result.type);
y(isnan(y)) = 0;
y = y';
end