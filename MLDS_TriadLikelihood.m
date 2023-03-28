function [L,P,D]=MLDS_TriadLikelihood(stim,response,p)
%[L,P,D]=TriadLikelihood(stim,response,p)
% 
%
% Computes the likelihood L of observing the responses (binary RESPONSE vector) to stimuli 
% indexed in the STIM vector given the perceptual scale values in P. 
%
% L is the average -log(probability), and P is the probability. D represents distances between the perceptual 
% values. 
%
% Uncomment the visualization part where perceptual scale values are monitored.

plotting = 0; % if plotting = 1, dann erstellst du auch plots

sigma         = p(end); % p is a vector containing different values on which function checks for optimization
p(end)        = [];
paramperdimen = 6;%the first point is located somewhere arbitrary, doesnt needt to float, can be fixed and it is fixed.
                %we have 8 stimuli of which two are fixed
                %([1;0;p;1;1;sigma], here our start value 1|0 and end value
                %1|1 are fixed
tdimen        = length(p)./paramperdimen; %divide each value of p by 6
                % p is set to 12, because we have 6 parameters in two dimensionsional space 
p             = reshape(p,tdimen,length(p)/tdimen);

% One point never moves around, we added it here
if tdimen == 1;
    %fix the first one to 0
    p     = [0 p 1];
elseif tdimen == 2%if we have more than any one dimensions, than fix it to point [0 1]...
    p     = [[1;0] p [1;1]]; %two dimensions expressed by 1|0 / 1|1
elseif tdimen == 3
    p     = [[1;0;0] p]; %three dimensions expressed by 1|0|0
end

%whether this is quadruplet or triad, we have different distance functions
%accordingly.
if size(stim,2) == 4
    for i = 1:size(stim,1)
        D(i) = norm(p(:,stim(i,2))'-p(:,stim(i,1))') - norm(p(:,stim(i,4))'-p(:,stim(i,3))');        
    end
elseif size(stim,2) == 3
    %D = abs(diff(p([stim(:,2) stim(:,1)]),1,2)) - abs(diff(p([stim(:,3) stim(:,2)]),1,2));
    for i = 1:size(stim,1)
        D(i) = norm(p(:,stim(i,1))'-p(:,stim(i,2))') - norm(p(:,stim(i,1))'-p(:,stim(i,3))');        
    end
end
%
P = normcdf(D(:),0,sigma);
chosenP = (P).^response(:).*(1-P).^(1-response(:));
L = -log(  chosenP );
PM = mean(chosenP); % average likelihood
L = mean(L); % calculates average of neg.log.likelihood

if plotting
    %plot the stuff online
    figure(6);
    if tdimen == 2
        plot(p(1,:),p(2,:),'o-')
    elseif tdimen == 1
        plot(p,'o-');
    else tdimen == 3
        plot(p(1,:),p(2,:),'o-')
    end
    axis square
    drawnow;
end