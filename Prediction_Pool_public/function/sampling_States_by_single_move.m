% -----------------------------------------------------------------------    
%   MH algorithm of Markov Switching of Lambda
% -----------------------------------------------------------------------    
    
% For common component
QPR = prob_old(1);    % Pr[St=0/St-1=0]@
PPR = prob_old(2);    % Pr[St=1/St-1=1]@

% Transition Probability Matrix    
% for one month
PR_TR = [QPR,  (1-PPR); (1-QPR),  PPR];

for i= T0_Forecast:(Tobs-h_Forecast)  % i = period    

    if (i < Tobs-h_Forecast)   
        if (St_old(i+1)==0)
            prob1_0 = QPR;              % prob(St+1=0|St=0)            
            prob1_1 = 1-PPR;            % prob(St+1=0|St=1)
        elseif (St_old(i+1)==1)
            prob1_0 = 1-QPR;            % prob(St+1=1|St=0)            
            prob1_1 = PPR;              % prob(St+1=1|St=1)
        end
        else
        prob1_0 = 1;                    % prob(St+1=0|St=0)            
        prob1_1 = 1;                    % prob(St+1=0|St=1)
    end   

    if (i > 1)
        if (St_old(i-1)==0)
            prob0_0 = QPR;              % prob(St=0|St-1=0)            
            prob0_1 = 1-QPR;            % prob(St=1|St-1=0)
        elseif (St_old(i-1)==1)
            prob0_0 = 1-PPR;            % prob(St=0|St-1=1)            
            prob0_1 = PPR;              % prob(St=1|St-1=1)
        end
        else
        prob0_0 = 1;                    % prob(St=0|St-1=0)            
        prob0_1 = 1; 
    end                        
    lik_0       =   (1-lambda0_old)*PredDen_SW(i) + lambda0_old*PredDen_KK(i) ;  % priorSt0 = 0; % mean_lik = mean([post_St0 post_St1]);
    post_St0    =   prob1_0*prob0_0*lik_0;
    lik_1       =   (1-lambda1_old)*PredDen_SW(i) + lambda1_old*PredDen_KK(i) ;  % priorSt1 = 0; % mean_lik = mean([post_St0 post_St1]);
    post_St1    =   prob1_1*prob0_1*lik_1;
    rate        =   min(1, post_St0/(post_St0+post_St1) );
    if (rand < rate)
        St_old(i) = 0;
    else
        St_old(i) = 1;              
    end                          

end 
     
save_prob(:,j) =[post_St0 post_St1 lik_0  lik_1 rate ];
     