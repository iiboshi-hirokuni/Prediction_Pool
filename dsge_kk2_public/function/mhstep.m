
function mhstep(sigmult,hessian,nblock,nburnin_block,nsim,cc0,cc,paras,vst,runpath,...
    appendmode,lastGoodBlock,resupath)

npara = length(paras);
para = zeros(npara,1);
for i = 1:npara
    para(i,1) = paras(i);
end
    

% best to display some model identification stuff here?

sigscale = sigmult;
sigpropinv = hessian;

sigpropdim = sum(vst.pmaskinv);

[u,s,v] = svd(sigpropinv);
sigproplndet = 0.0;

% for i = 1:npara
%     if i <= sigpropdim
%         s(i,i) = 1.0/s(i,i);
%         sigproplndet = sigproplndet + log(s(i,i));
%     end
% end

while 1
     
    paranew = para + cc0*(sigscale*randn(npara,1));
    paranew = paranew .* vst.pmaskinv + vst.pfix .* vst.pmask;
    
    [postnew,likenew] = mhstep_fcn(paranew,vst);
    
%     leva_new = ZZ_levarage*obsmean;
    
    propdens = 0;
%     propdens = -0.5*sigpropdim*log(2*pi) - 0.5*sigproplndet - ...
%        0.5*sigpropdim*log(cc0*cc0) - 0.5*(paranew-para)'*sigpropinv*(paranew-para)/(cc0*cc0);
    
    %if postnew > -10000;
    if postnew > -1E+5;
        break;
    end
    
end

newnblock = nblock - lastGoodBlock;

%tic;
for n = 1:newnblock
    
    parasim     = zeros(nsim,npara);
    likesim     = zeros(nsim,1);
    postsim     = zeros(nsim,1);
    rej         = zeros(nsim,1);
    propsim     = zeros(nsim,1);
    proppostsim = zeros(nsim,1);
    
%     leva_sim = zeros(nsim,nobs);
    
    if n == 1 && appendmode ~= 1
        
        newj = 2;
        likesim(1,1) = likenew;
        postsim(1,1) = postnew;
        propsim(1,1) = propdens;
        proppostsim(1,1) = postnew;
        parasim(1,:) = paranew';
        postold = postnew;
        likeold = likenew;
        paraold = paranew;
        
%         leva_sim(1,:) = leva_new;
%         leva_old = leva_new;
        
    elseif n == 1 && appendmode == 1
        
        newj = 1;
        postold = postnew;
        likeold = likenew;
        paraold = paranew;
        
%         leva_old = leva_new;
    
        
    else
        
        newj = 1;
        
    end
        
    
    for j = newj:nsim
        
        paranew = paraold + cc*(sigscale*randn(npara,1));
        paranew = paranew .* vst.pmaskinv + vst.pfix .* vst.pmask;
     
        [postnew,likenew] = mhstep_fcn(paranew,vst);
    
        propdens = 0;
%         propdens = -0.5*sigpropdim*log(2*pi) - 0.5*sigproplndet - ...
%             0.5*sigpropdim*log(cc*cc) - 0.5*(paranew-paraold)'*sigpropinv*(paranew-paraold)/(cc*cc);
    
        propsim(j,1) = propdens;
        proppostsim(j,1) = postnew;
        
        r = min(1,exp(postnew-postold));
        
        if rand < r
            postsim(j,1) = postnew;
            likesim(j,1) = likenew;
            parasim(j,:) = paranew';
            
            paraold = paranew;
            postold = postnew;
            likeold = likenew;
            
%            leva_sim(j,:) = leva_new;
%            leva_old = leva_new;
            
        else
            likesim(j,1) = likeold;
            postsim(j,1) = postold;
            parasim(j,:) = paraold';
            rej(j,1) = 1;
            
%             leva_sim(j,:) = leva_old;
            
        end
        
        if mod(j,nsim/2) == 0;
            disp(sprintf('block %d of %d : sim % d of %d',...
                (n+lastGoodBlock),nblock,j,nsim));
        end
        
    end
    
    disp(sprintf(' ' ));
    disp(sprintf('==================================================='));
    disp(sprintf('Rejection rate:    %f\n',mean(rej)));
    disp(sprintf('Likelihood:        %f\n',mean(likesim)));
    disp(sprintf('Posterior:         %f\n',mean(postsim)));
     % サンプリング前の調整係数の調整
    if n <= nburnin_block
       if  mean(rej) > 0.90
            cc = 0.1*cc;
       elseif  mean(rej) > 0.8
            cc = 0.3*cc;
       elseif  mean(rej) < 0.1
            cc = 2*cc;     
       elseif  mean(rej) < 0.5
            cc = 1.2*cc;
       end     
    end
    
    disp(sprintf('cc (adjustment coefficient): %f\n',cc ));
    disp(sprintf('==================================================='));
    
    fname = strcat(resupath, runpath,'/-block-',num2str(n+lastGoodBlock));
    
    save(fname, 'parasim', 'postsim', 'likesim', 'proppostsim', 'rej');
    
    
   
end

%disp(sprintf('elapsed time is %f minutes',(toc/60)));

