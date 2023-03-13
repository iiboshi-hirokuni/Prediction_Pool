nobs = size(series_YT,1);
ti = 1981:0.25:1981+(nobs-1)/4; 

figure(100)
subplot(4,3,1)
plot(ti, series_YT(:,1));title('Output Gap');
subplot(4,3,2)
plot(ti, series_YT(:,2)); title('Consumption');
subplot(4,3,3)
plot(ti, series_YT(:,4)); title('Investment');
subplot(4,3,4)
plot(ti, series_YT(:,4)); title('Wage');
subplot(4,3,5)
plot(ti, series_YT(:,5)); title('Labor');
subplot(4,3,6)
plot(ti, series_YT(:,6)); title('Inflation');
subplot(4,3,7)
plot(ti, series_YT(:,7)); title('Investment Price');
subplot(4,3,8)
plot(ti, series_YT(:,8)); title('Nominal Rate');
subplot(4,3,9)
plot(ti, series_YT(:,9)); title('Potencial GDP');
subplot(4,3,10)
plot(ti, series_YT(:,10)); title('Real Borrowing');


