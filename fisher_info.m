%% Ashley Lyons, Heriot-Watt University Physics Department, created on 21/02/17, v1.0
% This code calculates the expected Fisher information of a HOM dip or any
% fucntion for which the mean number of discrete events is known as a
%  single parameter is changed. The only noise assumed is that of the event
%  counting statistics i.e. shot noise

clear
close all
clc

N = 10000;      % Number of delay positions
x = linspace(1,N,N);        % Array for delay values

%% Mean Coincidence Count Distribution - Function for which the Fisher informatin is analysed
sig = N/6;      % Gaussian Width
m_c = 100;      % Max coincidences
n_c = 0.1*m_c;  % Min counts
% y = (m_c-(m_c-n_c)*exp((-((x-(N/2)).^2))/sig^2));        % Gaussian distribution
% y = (m_c+(m_c-n_c)*exp((-((x-(N/2)).^2))/sig^2));        % Gaussian distribution - positive
y = m_c - (m_c-n_c)*triangularPulse((N/2)-3000,(N/2)+3000,x);   % Triangular dip
% y = (m_c-(m_c-2*n_c)*exp((-((x-(N/2)).^2))/sig^2)) - (m_c/20)*(exp(((-((x-(N/2)).^2)))/sig^2)).*(cos(0.000002*N*x + N/2 + 1.5*N/50));        % Wiggly Gaussian distribution

%% Build Probability Density Function assuming Poissonian distribution & Calculate Fisher information
c_lim = 2.4*m_c;    % Max counts to build pdf for
pdf = zeros(c_lim,N);
for aa = 1:N
    for kk = 1:c_lim
        lambda = y(aa);
        %pdf(kk,aa) = ((lambda^kk)*exp(-lambda))/factorial(kk);
        pdf(kk,aa) = poisspdf(kk,lambda);
    end
end

% Find the Score (log likelihood)
score = zeros(c_lim,N-1);
for ll = 1:c_lim
    score(ll,:) = diff(log(pdf(ll,:)));
end

% Calculate Fisher Information
fi = zeros(1,N-1);
for bb = 1:N-1
    fi(bb) = sum(score(:,bb).^2.*(pdf(1:end,bb)));
end


%% Plots

% Plot pdf and mean coincidences
figure(1)
imagesc((x-500)/100,[1:c_lim],pdf)
hold on
plot((x-500)/100,y)

% Plot mean coincidences and Fisher information
figure(2)
plotyy(x,y,x(1:end-1),fi)
