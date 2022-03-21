p11 = predint(fff,xaxis,0.95,'observation','off');
p12 = predint(fff,xaxis,0.95,'observation','on');
p21 = predint(fff,xaxis,0.95,'functional','off');
p22 = predint(fff,xaxis,0.95,'functional','on');

subplot(2,2,1)
plot(fff,xaxis,yaxis), hold on, plot(xaxis,p11,'m--')
title('Nonsimultaneous Observation Bounds','FontSize',9)
legend off
   
subplot(2,2,2)
plot(fff,xaxis,yaxis), hold on, plot(xaxis,p12,'m--')
title('Simultaneous Observation Bounds','FontSize',9)
legend off

subplot(2,2,3)
plot(fff,xaxis,yaxis), hold on, plot(xaxis,p21,'m--')
title('Nonsimultaneous Functional Bounds','FontSize',9)
legend off

subplot(2,2,4)
plot(fff,xaxis,yaxis), hold on, plot(xaxis,p22,'m--')

legend({'Data','Fitted curve', 'Prediction intervals'},...
       'FontSize',8,'Location','northeast')
clear title