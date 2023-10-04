close all, clear vars, clc, format long, clf
%% Data
LB = -50; % lower bound of the space domain
UB = 50;  % upper bound of the space domain   
% Parameters of the initial PDF
mu_f0 = 0;
sigma_f0 = 5;
T = 1; % final time   
NT = 500; % #of intervals in the time grid
% Rate and parameter distribution of the jumps
lambda = 5;
mu_jumps = 0; 
sigma_jumps = 3;  
% Temporal mesh
dt =  T / NT;
timeArray = dt : dt : T;

u = 0 .* timeArray;

% Number of Monte Carlo runs
n_MCruns = 4;
%% Monte Carlo runs
waitbar_MC = waitbar(0,'Monte Carlo runs');
for k = 1:n_MCruns % for each run, W and P are simulated
    Xiniz = 0;
    Xtemp = Xiniz;
    Xvalues = zeros(1,NT);
    waitbar(k / n_MCruns)
    % 1 Brownian increments
    dW = sqrt(dt) * randn(1, NT);    
    % 2 Poisson process (alternatively, one can simulate only the increments)   
    
    % if MATLAB poissrnd available (Stat. and Machine Learning Toolbox) use this
    % nJumps = poissrnd(lambda*(T)); % number of jumps in the interval [0, T]
    % Otherwise: 
    
    nJumps = funPoissonRecursive(lambda*(T)); % number of jumps in the interval [0, T]

    PoisProc = zeros (1, NT);
    if nJumps>0
        jumpAmplitudes = zeros(1,nJumps);
        Bjumps = zeros(nJumps, NT);
        jumpTimes = sort ( unifrnd0(0, T, 1, nJumps) ); % sorted array containing the jump times
        for cont_jumps = 1:nJumps   
            jumpAmplitudes(cont_jumps) = trunc_norm(mu_jumps, sigma_jumps, LB-UB, UB-LB);
        end  
        for cont_jumps=1:nJumps
            Bjumps(cont_jumps,:) =  jumpAmplitudes(cont_jumps) * stepfun( timeArray, jumpTimes(cont_jumps));
        end    
        Bjumps = cumsum(Bjumps,1);    
        PoisProc = Bjumps(nJumps,:);
    end     
    % Euler-Maruyama method
    Pinc = PoisProc(1);
    for j = 1:NT
        while j>1
            Pinc = PoisProc(j) - PoisProc(j-1) ;
            break
        end       
        Xtemp = Xtemp + b_drift(Xtemp, u(j) ) * dt + ... 
            sigma_diff(Xtemp, u(j)) * dW(j) + ... 
            cjump(Xtemp,timeArray(j)-dt ) * Pinc; 
        % Reflecting boundary condition
        while (Xtemp < LB || Xtemp > UB)  
            Xtemp = (2*UB - Xtemp) * (Xtemp > UB) + ...
                    (2*LB - Xtemp) * (Xtemp < LB);
        end
        Xvalues(j) = Xtemp;
    end        
    figure(1)
    plot([0, timeArray], [Xiniz, Xvalues],'LineWidth',1.6)
    %axis([0 T -20 20])
    hold on
end
close(waitbar_MC)

ybottom = 0;
yline(ybottom,'b','Linewidth',1);

print('-depsc2', 'JD01.eps','-b0'); 
%print('-dpdf', 'JD01.pdf','-b0');


%% function truncated normal
function [y] = trunc_norm(mu, sigma, a, b)
    pb = normcdf0( (b - mu) / sigma );
    pa = normcdf0( (a - mu) / sigma );
    diff_val = pb - pa;
    y = mu + sigma * norminv(diff_val * rand + pa);
end
%% Function b(x,t) (drift coefficient - derivative of C)
function y = b_drift(x,u)
    y = -4.*x ;
end
%% Function sigma(x,t) (diffusion coefficient)
function y =  sigma_diff(x,t)
 %   y = x + t;
        y = 2.;
end

%% Function cjump(x,t) (diffusion coefficient)
function y =  cjump(x,t)
 %   y = 2.+ sin(2*pi*t);
        y = 1.0 ;
end


