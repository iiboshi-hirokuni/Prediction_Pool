%@=========================================================================@
%@====PROCEDURE THAT GENERATES ST==========================================@
%@=========================================================================@

function [S_T]= HamiltonFilter( PredDen_Model_A, PredDen_Model_B, lam0, lam1, ...
                                prob_old,T0_Forecast,h_Forecast,h)

% YSTAR=YTT;
   Tstar = size(PredDen_Model_A,1);     % number of sample
   FLT_PR = zeros(2,Tstar);
   
   DMNSION = 4;   %  2 state (t period) time 2 state (t-1 period)

          % For common component
                  QPR = prob_old(1);    % Pr[St=0/St-1=0]@
                  PPR = prob_old(2);    % Pr[St=1/St-1=1]@
         
         % Transition Probability Matrix            
                     PR_TR= [ QPR    (1-PPR);
                             (1-QPR)  PPR   ];
                         
                     PR_TR_1 =   PR_TR;  
                 if h >1
                     for i = 2:h    
                       PR_TR =  PR_TR_1*PR_TR;  
                     end
                 end    

   coef_SW = [1-lam0; 1-lam1; 1-lam0; 1-lam1 ];
   coef_KK = [lam0; lam1; lam0; lam1 ];

   %@<<<<<<<<<<<<<<<<<<<START FILTERING>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>@

           PR_TRF = vec(PR_TR);

           % initialization  
           A = [(eye(2) - PR_TR); ones(1,2)];
           EN = [0;0;1];           
           PROB__T = inv(A'*A)*A'*EN;    % PR[S_t=0]|PR[S_t=1],
                                        % 2x1 steady-state PROBABILITIES@
           PROB__T0 = PROB__T;
           PROB__= vec([PROB__T PROB__T]);

   %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++@
   %  Hamilton Filter     ++++++++++++++++++++++++++++++++++++++++++++++++@
   %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++@

    for J_ITER = T0_Forecast:(Tstar - h_Forecast)

          f_cast1 = PredDen_Model_A(J_ITER,1)*ones(DMNSION,1).*coef_SW ...
                  + PredDen_Model_B(J_ITER,1)*ones(DMNSION,1).*coef_KK ;

          prob_dd=PR_TRF .* PROB__;
          
          likelihood =  f_cast1;

          pr_vl=likelihood.*prob_dd;       % PR[St-m,...,St-1,St|Y_t] @    

          pr_val = sum(pr_vl,1);

          PRO_=pr_vl/pr_val;

          PROB__T = PRO_(1:DMNSION/2,1)+PRO_(DMNSION/2+1:DMNSION,1);
                                            % Pr[St-m+1,...,St/Yt]@

          PROB__ = vec( [PROB__T PROB__T] );

          FLT_PR(:,J_ITER) = PROB__T;        % Pr[St|Yt], 2x1@

%    J_ITER = J_ITER+1;
    end

  

 %=============== GENERATE S_T =====================================@
 %   Smoothing  
 %
   
  % Set State S_t  at terminal period T
         S_T = zeros(Tstar,1);
         S_T(Tstar,1)=   bingen( FLT_PR(1,Tstar), FLT_PR(2,Tstar),1);

  % Set State S_t  from period T-1 to period 1 
   
     J_ITER=Tstar-1;
     while (J_ITER >= 1);

            if S_T(J_ITER+1,1)==0
                    P0=QPR*FLT_PR(1,J_ITER);
                    P1=(1-PPR)*FLT_PR(2,J_ITER);
                    
            else S_T(J_ITER+1,1)==1;
                    P0=(1-QPR)*FLT_PR(1,J_ITER);
                    P1=PPR*FLT_PR(2,J_ITER);
            end;

            S_T(J_ITER,1) = bingen(P0,P1,1);

         J_ITER=J_ITER-1;
     end;

     
   %--------------------------------------------------------------
   
    function [b]=vec(a)
        
     [m,n]=size(a);
     b=[];
     
      for i =1:1:m
        b = [b ; a(i,:)'];
      end  
      
      %-----------------------------------------------------------
    function [s] = bingen(p0,p1,m);

      pr0 = p0/(p0+p1);       %/* prob(s=0) */
      u = rand(m,1);
      
      if u > pr0
           s = 1;
      else
           s = 0;
      end     
