function [h,d,pC] = JG(sigma_y,E1,E2,nu1,nu2,eta,beta,sigma,sk,ku)

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

yS     = 0.045944 / (eta * beta * sigma);
sigmaA = sqrt(sigma^2 - 3.717E-4 / (eta^2 * beta^2));

% error check
if real(sigmaA) <= 0
    error('Invalid parameters.')
end

result = f_johnson_M(0,1,sk,ku);
verteilung_surf = @(z) Distr(z,result);
verteilung_asp = @(z) sigma/sigmaA * verteilung_surf(z * sigma/sigmaA);

E = ((1 - nu1^2)/E1 + (1 - nu2^2)/E2)^(-1);

e_y = sigma_y / E;
B = 0.14*exp(23*e_y);
C = 1.295*exp(0.736 * nu1);
Fc = 4/3 * (beta/E)^2 * (0.5*C*pi*sigma_y)^3;

omega1 = (pi*C*sigma_y/(2*E))^2 * beta / sigma;

h = 0:0.01:4;
d = h - yS;

pC = zeros(1,length(h));

for i = 1:length(h)
    
    w_e  = @(z) (z-d(i)).^1.5 .* verteilung_asp(z);
    w_ep = @(z) ((exp(-0.25 .* ((z-d(i))./omega1).^(5/12))) .* ((z-d(i))./omega1).^1.5 +...
        4/C .* 2.84 .* (1 - exp(-0.82.*(sqrt((z*sigma-d(i)*sigma)/beta) .* ((z-d(i))./(1.9*omega1)).^(B/2)).^(-0.7))) .*...
        (1 - exp(-0.04 .* ((z-d(i))./omega1).^(5/9))) .* ((z-d(i))./omega1)) .* verteilung_asp(z);
    
    pC(i) = eta * E * beta * sigma * 4/3 * sqrt(sigma/beta) * integral(w_e,d(i),d(i)+1.9*omega1) +...
        Fc * eta * integral(w_ep,d(i)+1.9*omega1,Inf);

end

end


function y = Distr(z,result)
y = f_johnson_pdf(z',result.coef,result.type);
y(isnan(y)) = 0;
y = y';
end