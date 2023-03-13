%=======================================================
%
%     Forecasting by a Bayesain VAR model 
%
%=======================================================

clear;

% case_data_rich = 4;
% 
% if case_data_rich == 1
%    YY_f =  xlsread( 'm1021FC_SW.xls');
% elseif case_data_rich == 2
%    YY_f =  xlsread( 'm1021FC_A.xls');   
% elseif case_data_rich == 3
%    YY_f =  xlsread( 'm1021FC_B.xls');
% elseif case_data_rich == 4
%    YY_f =  xlsread( 'm1021FC_C.xls');    
% end   

% var_fore = 2;   % output=1, inflation=2 %
%    
% data = csvread('data_rich_80-98-1.csv',1,1);

YYact = YY(:,1);

YYact_f = forecast_m(:,1);

[nsim,h] = size(YY_f);
nburn = 0;
n = 1;
p=0;

% [T,k]=size(XXact);  % number of sample 
% 
%  h = 12;           % horizon of forecasting
 
% YY_f = zeros(nsim-nburn,h,n);
% XX_f = zeros(nsim-nburn,h,n*p+1);


%*********************************************
%
%  �t�@���`���[�g�̍쐬
%
%*********************************************

a  =0.90;   % (1-2*rate) percentage of Credibile interval of Parameters  
  rate = (1-a)/2;
    
m = 10;   %  �t�@���`���[�g�̊K���̐�
b = a/2/m; % �t�@���`���[�g�̊e�K���̕� 

     % �t�@���`���[�g�̐F�̐ݒ�      
      col1 = 0.85; %    ���̒����@ 
      col2 = 0.75; %    ���̒����@
      col3 = 0.3; %  col3< col1  ���E���F�̒����@�傫���ƍ����Ȃ�

  
sort_forecast = zeros( nsim-nburn, h, n );


for k = 1:1:n 
  for i=1:1:h
     sort_forecast(:,i,k) = sort(YY_f(:,i,k),1);
  end
end  
  
for i =1:n             %  i = �n�񐔂̐�
  %  Graph of Forecast 
  YY_mean = mean(YY_f(:,:,i),1)';
  
  
  YY_band = zeros(h,2*m);
  sort_forecast_1 = squeeze(sort_forecast(:,:,i));
  YY_median = sort_forecast_1(round((nsim-nburn)*0.5),:)'; 
  
  for j = 1:m 
      
      YY_band(:,j) = sort_forecast_1(round((nsim-nburn)*(rate+b*(j-1))),:)'; 
      YY_band(:,2*m-j+1) = sort_forecast_1(round((nsim-nburn)*(1-(rate+b*(j-1)))),:)';
      
%       rate+b*(j-1)
%       1-(rate+b*(j-1))
  end
  
  YY_f_save = YY_band(:,1);
  for j = 1:m-1
    YY_f_save = [ YY_f_save YY_band(:,j+1)-YY_band(:,j) ];
  end  
    
    YY_f_save = [ YY_f_save YY_median-YY_band(:,m) YY_band(:,m+1)-YY_median ];
    
   for j = m:-1:2
     YY_f_save = [ YY_f_save YY_band(:,2*m+2-j)-YY_band(:,2*m+1-j) ];
   end
 
  if case_data_rich == 1  
         figure('Name','Forecast of GDP (Case SW)','Color','w');
  elseif case_data_rich ==2
         figure('Name','Forecast of GDP (Case A)','Color','w');
  elseif case_data_rich==3
        figure('Name','Forecast of GDP (Case B)','Color','w');             
  elseif case_data_rich==4
        figure('Name','Forecast of GDP (Case C)','Color','w');
  end
          
      tt = size([YYact(:,1);YY_mean],1);  
%        tt = size([ YY_mean],1); 
     ti=linspace(1980-p,1980-p+tt*0.25,tt)';
  
 
 hold on 
  
  hh = area(ti, [YYact(:,i) zeros(size(YYact(:,1),1),2*m);...
                YY_f_save ] ) ;
           
%  hh = area(ti, [    YY_f_save ] ) ;
      
 
           
      for j = 1:m    
         set(hh(j),'FaceColor',[col1*(1-(j-1)*1/m)+(1-col1) col2*(1-(j-1)*1/m)+(1-col2) col3*(1-(j-1)*1/m)+(1-col3) ])
    
         set(hh(2*(m+1)-j),'FaceColor',[col1*(1-(j)*1/m)+(1-col1) col2*(1-(j)*1/m)+(1-col2) col3*(1-(j)*1/m)+(1-col3) ])

      end
         set(hh(m+1),'FaceColor',[1-col1 1-col2 1-col3 ])
         set(hh(1),'FaceColor',[1 1 1 ])
      
      %  �r���̐ݒ� 
      set(hh,'LineStyle','none') % Set all to same value
%        set(hh,'LineStyle','none','LineWidth',0.1) % Set all to same value

       ylim([-7 7]);

      %  �\���l
      plot(ti, [YYact(:,i); YY_median ],'LineStyle','-','Color','r', 'LineWidth',2.5);    

      
      % ���ђl
       plot(ti, [YYact(:,i); YYact_f ],'LineStyle','-','Color','black', 'LineWidth',2.5);  
  
         
      hold off
end
  
 