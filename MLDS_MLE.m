function [Estimate,ExitFlag,Likelihood]=MLDS_MLE(StimList,R,ModelDimen)
%[Estimate,ExitFlag]=MLDS_MLE(StimList, R, ModelDimen);
%
%
% This is the core optimization function.
% You need to optimize your Likelihood calculation , because otherwise you
% will end up with very small values, which are difficult to interpret...
%
% STIMLIST is the list of the stimulus indices. R is the binary response variable. MODELDIMEN can be 1, 2 or 3.
%
% ESTIMATE contains the coordinates of the perceptual scale values in MODELDIMENx8 
% matrix.
%
% Depends on MLDS_TriadLikelihood


options   = optimset('Display','off','maxfunevals',100000,'tolfun',10^-8,'tolx',10^-8,'tolcon',10^-8,'maxiter',8000);
%Function that will be optimized... It is a wrapper aroudn
%TRIADLIKELIHOOD so that the initial values in [P0 sigm] are passed
%to it easily
funny     = @(x) MLDS_TriadLikelihood(StimList, R, x);
%
X    = [];

%define the upper and lower boundaries for the optimization as well
%as the initial parameters
if ModelDimen == 1
%     p0          = linspace(0,1,8)+rand(1,8)*0.05;
%     p0          = p0./max(p0);
%     p0          = p0(:);%make a column vector
%     %p0         = linspace(0,1,8)*2*pi;
%     p0          = p0(2:end);%the first entry is inserted in the likelihood function and kept 0 position.           
    %
% %     LB          = [zeros(1,7)-0.1 0.01];%7 instead of 8 because we fix the first one to 0
% %     UB          = [ones(1,7)+0.1 2];
 p0 = randn(6, 1);
 
elseif ModelDimen == 2
    %vielleicht vernünftige startwerte mit MDS
    %function mdscale(X)
    %X = matrix, NxN, N = Anzahl stimuli
    %X(i, j) = Dissimilarity zwischen stimulus i und j
    %abs(rotation(i) - rotation(j))


    %p0          = [cos(linspace(0,2*pi-2*pi/8,8)) ; sin(linspace(0,2*pi-2*pi/8,8))];
    %p0          = p0(:,2:end-1);%the first entry is inserted in the likelihood function and kept at [1 0] position.
    %p0          = p0(:);%make it a column vector
    
    %12 random werte aus standardnormalverteilung (mean = 0, sd = 1)
    p0 = randn(12, 1); % R = randn(N,M) returns an N-by-M matrix (double) containing pseudorandom values drawn from the standard normal distribution
    
    %
%     LB          = -ones(1,7)-0.1;
%     LB          = [LB(:) ;0.01];
%     UB          = ones(1,14)+0.1;
%     UB          = [UB(:) ;10];
elseif ModelDimen == 3;
    p0          = [cos(linspace(0,2*pi-2*pi/8,8)) ; sin(linspace(0,2*pi-2*pi/8,8)) ; rand(1,8)];
    p0          = p0(:,2:end);%the first entry is inserted in the likelihood function and kept 0.    
    p0          = p0(:);%make it a column vector   
    %
% %     LB          = repmat([-ones(1,7) 0.01],3,1);
% %     UB          = repmat([ones(1,7) 2],3,1);
end
%% estimate initial sigma value
c=0;
l = nan(1,100);
x = linspace(0.01,6,100);
for s = x
    c=c+1;
    l(c) = [funny([p0 ;s])];
end
[mi i] = min(l);
%figure(100);plot(x,l)
sigma  = x(i);
%fprintf('Estimated Sigma is %.2f\n', sigma);
sigma = 1;
%ML Estimate of perceptual scale values
[Estimate,Likelihood,ExitFlag]  = fminsearch(funny , [p0(:) ;sigma], options); %funny receives the STIMLIST AND Response vector automatically.

sigma = Estimate(end);
p = Estimate(1:end-1);

%val1 = 3; 
%val2 = 6; % index to insert

%fix1 = 0;
%fix2 = 1;
    %1D
if ModelDimen == 1
    p = reshape(p, 1, 6);
    %p = [p(1:val1-1) fix1 p(val2:end)] % insert
    p = [p(1:3), 0, p(4:5), 1, p(end)];
    
    %2D
elseif ModelDimen == 2
    p = reshape(p, 2, 6);
    p = [p(:,1:3), [1;0], p(:, 4:5), [1;1], p(:, end)]; % fixate position 4 and 7 (GS2 and GS5) 2D
end;

%Add the first point that is immobile

if ModelDimen == 2
   % Estimate = [1;0;p;1;1;sigma];
    Estimate = [p(:);sigma];
elseif ModelDimen == 1
%     Estimate = [0;p;1;sigma];
    Estimate = [p(:);sigma];
end;
    
    
    
    
    
    
    
    
    
    
    