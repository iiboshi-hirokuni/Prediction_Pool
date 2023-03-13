function [ re_lambda, log_lik,mean_lambda ]=...
    ParticleFilter_par(PredDen_SW, PredDen_KK, T0_Forecast, h_Forecast, particles, rho, mu, sigma, rv) 
%     ParticleFilter    (PredDen_SW, PredDen_KK, T0_Forecast, h_Forecast, particles, rho,mu, sigma) 
    
% particles = nsim;

t       = size(PredDen_SW,1);                 % Number of period
xt      = zeros(t, particles);                  % Particles of factors (3*1) at period t
lambda  = zeros(t, particles);
et      = zeros(t, particles);                  % Shock e_t ~ N(0,1) coresponding to each particles h_t
wt      = zeros(t, particles);                  % Weight of each particles
xt_temp = zeros(1,particles);
      
% Generation of perticles of xt, et, wt at period 1 
parfor j = 1:1:particles
    %   Transition Eq.: x_t = rho* x_t-1 + sqrt(1-rho^2)*v_t,   v_t -- N(0, 1)
%     xt_1            =  (1-rho)*mu+rho*rand(1,1) + sqrt(1 - rho^2)*sigma*randn(1,1); 
     xt_1            =  (1-rho)*mu+rho*rand(1,1) + sqrt(1 - rho^2)*sigma*rv(j,1); 
    
    %   Measurement Eq.
    pd0              =   'Normal';
    lambda_1        =   cdf(pd0,xt_1,0,1);
    xt_temp(:,j)    =   xt_1;                    % xt(j) ->  particles( j )
    
    %  likelihood function   
    wt(T0_Forecast, j) = (1-lambda_1)*PredDen_SW(T0_Forecast) + lambda_1*PredDen_KK(T0_Forecast) ;  % wt(j) -> weight of particles(i); xt(j) 
end

% Calculating likelihood at period 1
prob_x  =  wt(T0_Forecast,:);
lik_1   = mean(wt(T0_Forecast,:));        
log_lik = log(lik_1);

%  Resampling xt(j) with Prob (prob_x(j)) at period 1
parfor j=1:1:particles                            
    xt(T0_Forecast,j)       = mnrnd(1, prob_x/sum(prob_x))*xt_temp';
    pd1              =   'Normal';
    lambda(T0_Forecast,j)   = cdf(pd1,xt(T0_Forecast,j),0,1);
end 

% Generation of perticles of xt, wt at period 2 through T
% t is period from 1 to T 
for i = T0_Forecast+1:1:t-h_Forecast                   
    parfor j=1:1:particles     % j is j-th of particles
    % Step 1:  generating particle x_t        
    %   Transition Eq. - x_t = rho* x_t-1 + sqrt(1-rho^2)*v_t,   v_t -- N(0, 1)
    xt_1        = xt(i-1,j);                                       % factors at period t-1         
%     xt_t        = (1-rho)*mu+rho*xt_1' + sqrt(1-rho^2)*sigma*randn(1,1);            % generation of factors at period t
    xt_t        = (1-rho)*mu+rho*xt_1' + sqrt(1-rho^2)*sigma*rv(j,i);    
    pd2              =   'Normal';
    lambda_t    = cdf(pd2,xt_t,0,1);
    xt_temp(:,j) = xt_t;

    % Step 2: Calculating Probabilities corresponding to particles 
    %   Measurement Eq.- weight of likelihood p(y_t|x_t) for each particle h_t
    wt(i, j) = (1 - lambda_t)*PredDen_SW(i) + lambda_t*PredDen_KK(i); 
    end
    prob_x =  wt(i,:);     % Set probability of each particle
    
    % Step 3: Resampling xt(j) with Prob (prob_x(j))
    parfor j=1:1:particles             
        xt(i,j)     = mnrnd(1, prob_x/sum(prob_x)) * xt_temp';
        pd3              =   'Normal';
        lambda(i,j) = cdf(pd3,xt(i,j),0,1);
    end
    
    % Step 4: Calculating likelihood
    lik_t   = mean(wt(i,:));         
    log_lik = log_lik + log(lik_t);  % likelihood of full period
    
     %==     display iteration number ===============     
%      itr = num2str(i) ;
%      if mod(i, 5)==0 
%              disp([ itr 'th  period.......']);
%      end        
     %===============================================        
end   
  
% Sampling one partilce from the distribution   
re_lambda = lambda(:,:) ;

mean_lambda = mean(lambda(:,:),2);
  
 
  
  

  
  
  
  
  