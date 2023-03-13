function   [wt,xt_temp]= sub_parallel_2model(P, PredDen, shock,xt_Resamp,particles,t,nn)  

  wt      =  zeros(1,particles); 
  xt_temp =  zeros(nn,particles); 
  num     =  particles; 

 rho = P.rho; 
 mu  = P.mu;
 sigma = P.sigma;
 %  rv  = P.rv;
% rho_mat = [rho ];
 rho0_mat = [sqrt(1 - rho^2)];
 
 for j=1:num 
    % Step 1:  generating particle x_t        
    %   Transition Eq. - x_t = rho* x_t-1 + sqrt(1-rho^2)*v_t,   v_t -- N(0, 1)
    xt_1        = xt_Resamp(1,j);    
         
    xt_t        = (1-rho)*mu+rho*xt_1' + sqrt(1-rho^2)*sigma*shock(j);  
    xt_t        = real(xt_t);
    pd2              =   'Normal';
    lambda_t    = cdf(pd2,xt_t,0,1);
            
    xt_temp(1,j) = xt_t;               % xt(j) ->  particles( j )
    xt_temp(2,j) = lambda_t; 
 
    % Step 2: Calculating Probabilities corresponding to particles 
    %   Measurement Eq.- weight of likelihood p(y_t|x_t) for each particle h_t
     wt(j) = (1 - lambda_t)*PredDen(t,1) + lambda_t*PredDen(t,2); 
  
     if isfinite(wt(j))==0
          wt(j) = 0;
     end    
       
  end
   
 
 
     