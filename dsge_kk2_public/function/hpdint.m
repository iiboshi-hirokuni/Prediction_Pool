function hpdband = hpdint(draws,percent)

%drawdim   = cols(draws);
%ndraws    = rows(draws);
[ndraws, drawdim] = size(draws);
hpdband   = zeros(2,drawdim);
nwidth    = floor(percent*ndraws);

for i = 1:drawdim

   drawcoli = draws(:,i);
   % sort response for period i, element 1 is max
   
   drawcoli = flipud(sort(drawcoli));
   bup   = 1;
   minwidth  = drawcoli(1) - drawcoli(nwidth);
   %j = 2;
   %do until j > (ndraws-nwidth+1);
   for j = 2:(ndraws-nwidth+1)
      newwidth = drawcoli(j) - drawcoli(j+nwidth-1);
      if newwidth < minwidth
         bup = j;
         minwidth = newwidth;
      end;
   %   j = j+1;
   %endo;
   end
   hpdband(2,i) = drawcoli(bup);
   hpdband(1,i) = drawcoli(bup+nwidth-1);

end
