%/* filename:    nkbcpri1.g
%** description: The program generates draws from the prior of the 
%**              DSGE model
%*/

block_n    = 10;
block_size = 100;

parasim_out = [];

for block_it = 1:block_n
   
   % Generate THETA draws
   parasim_t = zeros(block_size,npara);
   
   for i = 1:npara
   
      if pshape(i) == 1;    % BETA Prior
         for j = 1:block_size
             a = (1-pmean(i))*pmean(i)^2/pstdd(i)^2 - pmean(i);
             b = a*(1/pmean(i) - 1);
                  parasim_t(j,i) = rndb(a,b); %betarnd(a,b);  
         end
      
      elseif pshape(i) == 2; % GAMMA PRIOR
         b = pstdd(i)^2/pmean(i);
         a = pmean(i)/b;
         %parasim_t(:,i) = b*gamrnd(a,1,block_size,1);
         for j = 1:block_size
             parasim_t(j,i) = b*rg1(a);
         end
         
      elseif pshape(i) == 3; % GAUSSIAN PRIOR
         a = pmean(i);
         b = pstdd(i);
         parasim_t(:,i) = a + b*randn(block_size,1);
      
      elseif pshape(i) == 4; % INVGAMMA PRIOR
        
         a = pmean(i);
         b = pstdd(i);
         
         for j = 1:block_size
              parasim_t(j,i) = sqrt( b*a^2/sum(randn(b,1).^2)  );
         end  
          
%          for j = 1:block_size
%             nn = 2;
%             d = pstdd(i)^2/ pi ;
%             c = rndc(nn);
%             t2 = c/d;
%             parasim_t(j,i)= sqrt( 1/t2 );
%          end
      
      elseif pshape(i) == 5; % uniform prior
         a = pmean(i);
         b = pstdd(i);
         parasim_t(:,i) = a + (b-a)*rand(block_size);
      
      elseif pshape(i) == 0; % No Prior - Fixed
         a = pmean(i);
         parasim_t(:,i) = a*ones(block_size,1);
      end
   end
   parasim_out = [parasim_out; parasim_t];
end

oparapri = strcat(priopath,runname,'pa0.csv');
util_csvwrite(oparapri, parasim_out, para_names);

%---------------------------
% filename:    lpprior2.g
% description: The program deletes the draws that lie 
%              outside of the parameter space
%---------------------------

iparapri = strcat(priopath,runname,'pa0.csv');
thetasim = csvread(iparapri, 1, 0);

nread   = size(thetasim,1);

invaliddraw = zeros(nread,1);
detsim      = zeros(nread,2);

for j = 1:nread
    % For each element of the block, compute the desired statistic
    [T1,TC,T0,RC] = dsgesolv(thetasim(j,:)');
    if (RC(1) == 1) && (RC(2)==0)
        % indeterminacy */
        detsim(j,2) = 1;
    
    elseif (RC(1) == 1) && (RC(2)==1)
        % determinacy */
        detsim(j,1) = 1;   
    else
        invaliddraw(j,1) = 1;  
    end
end

thetasim_0 = thetasim;
thetasim = delif(thetasim, invaliddraw);
detsim   = delif(detsim, invaliddraw);  
thetasim = delif(thetasim, detsim(:,2));
nvalid  = size(thetasim, 1);

oprisim = strcat(priopath,runname,'pa.csv');
util_csvwrite(oprisim, thetasim, para_names);

disp('Fraction of valid draws')
nvalid/nread




