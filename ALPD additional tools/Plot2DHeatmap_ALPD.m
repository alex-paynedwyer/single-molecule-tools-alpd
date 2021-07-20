%Plot 2D heatmap for Stoich v Diff (one colour, linked or unlinked) or Stoich v Stoich. (two colours, paired).

%ticks interval in number of molecules
interval = 10;

%Set kernel density width
%params.stoichKDW = 0.7;

%load in output first
%pairedstoichs = output.PairedStoichsList;

%'pairedstoichs' needs to be a list of stoichiometries linked to each other, like the following
pairedstoichs = [];
pairedstoichs(:,1) = trackArrayCh1(trackArrayCh1(:,8)>0,2)/params.IsingleCh1;
pairedstoichs(:,2) = trackArrayCh2(trackArrayCh1(trackArrayCh1(:,8)>0,10),2)/params.IsingleCh2;

%Paired scatterplot

scatter(pairedstoichs(:,1),pairedstoichs(:,2),'filled')
xlabel('DnaQ linked stoichiometry')
ylabel('PriC linked stoichiometry')
gradient=pairedstoichs(:,1)\pairedstoichs(:,2);
%hold on

%SHOW LINE OF BEST FIT
%plot([0;pairedstoichs(:,1);max(pairedstoichs(:,1))*1.1],[0;pairedstoichs(:,1)*gradient;max(pairedstoichs(:,1))*1.1*gradient],'--','Color','black')
%disp(strcat('Gradient = ',num2str(gradient)))

axis tight
box on

xsmax = interval*ceil(max(pairedstoichs(:,1))/interval);
ysmax = interval*ceil(max(pairedstoichs(:,2))/interval);
xlim([0,xsmax])
ylim([0,ysmax])
xticks(0:interval:xsmax)
yticks(0:interval:ysmax)
pbaspect([xsmax/ysmax 1 1]);
%pbaspect([1 1 1]);

%LABEL GRADIENT ON SCATTERPLOT
%text(max(pairedstoichs(:,1))*1,max(pairedstoichs(:,2))*1,strcat('Gradient: ',num2str(round(gradient,2))),'HorizontalAlignment','right');

% 2D Kernel density plot
figure
xhmax = 300;
yhmax = 4000;
%xhmax = interval*ceil(max(pairedstoichs(:,1))/interval);
%yhmax = interval*ceil(max(pairedstoichs(:,2))/interval);
gridx = 0:params.stoichKDW/10:xhmax+params.stoichKDW/10;
gridy = 0:params.stoichKDW/10:yhmax+params.stoichKDW/10;
[xsupp,ysupp] = meshgrid(gridx, gridy); xsupp = xsupp(:); ysupp = ysupp(:);
support = [xsupp ysupp];
reflect = 0;
%[plot2D,KDFstoichs2D,KDFx2D] = KDFplotH_2D(pairedstoichs,support,params.stoichKDW,reflect);
KDFplotH_2D(pairedstoichs,support,params.stoichKDW,reflect);
set(gca,'View',[0,90])
surf1=get(gca,'Children');
set(surf1,'EdgeColor','None')
axis tight
box on

xlim([0,xhmax])
ylim([0,yhmax])
xticks(0:interval:xhmax)
yticks(0:interval:yhmax)
pbaspect([xhmax/yhmax 1 1]);
%pbaspect([9/30 1 1]);
xlabel('DnaQ linked stoichiometry')
ylabel('PriC linked stoichiometry')
