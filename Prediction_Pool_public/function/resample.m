function [ re_xt,  re_prob  ] = resample(xt, prob_xt, nparticles)

 w = prob_xt/sum(prob_xt);
% w = prob_xt;

total = 0;

re_xt = xt;
re_prob = prob_xt;

for i = 1:nparticles
    
    n = floor(w(i)*nparticles);
    u = (w(i)*nparticles)-floor(w(i)*nparticles);
    if u > rand(1)
        n = n+1;
    end    
    
    if (n>0)&&(total<nparticles)
        
        if (total + n) > nparticles
            total_last = nparticles;
        else
            total_last = total + n;
        end
        
        for j = total+1:total_last
            re_xt(:,j)   = xt(:,i);
            re_prob(:,j) = prob_xt(:,j);
        end
    end
    
    total = total + n;
%     total = total_last ;
end