%% Radial intensity distributions for Gaussian excitation of simple fluorophores

% R: Ratio of radius to 1/e^2 beamwaist (support)
NumR = 10000; %linearly spaced values of R between
Rmin = 0.01; % and
Rmax = 2;

%Iprime: Ratio of excitation intensity to saturation intensity
NumI = 10; %linearly spaced values of Iprime between
Imin = 0.5; %and
Imax = 5;

%Max correction factor
Cmax = 3;

R0 = linspace(Rmin,Rmax,NumR);  %support for plotting radial functions;
R = repmat(R0,NumI,1); %supports to calculate multiple instances of I'
Iprime0 = linspace(Imin,Imax,NumI);
Iprime = repmat(Iprime0',1,NumR);

fex = exp(-2.*R0.^2);  %excitation distribution is Gaussian
fem = (1./Iprime+1)./(1./Iprime./exp(-2.*R.^2)+1);  %saturated emission distribution

fexmean = (1-fex)./R0.^2/2;   %means of those distributions inside a circle
femmean = (1./Iprime+1).*(1+log(fem)./R.^2./2);

fexcov = sqrt((1-exp(-4.*R0.^2))./(4*R0.^2)./fexmean.^2-1); %coefficient of variance of the distributions inside a circle
beta = (1./Iprime+1)./femmean;
femcov = sqrt((beta-1).^2+beta.*(beta-2).*log(fem)./(2.*R.^2)-beta.*(1-fem)./(2.*femmean.*R.^2));
fexcov(fexcov.^2<0) = 0;  %protect against rounding errors
femcov(femcov.^2<0) = 0;

corrfex = 1./fex;   %correction factors for flattening distributions
corrfem = 1./fem;
corrfex(corrfex>Cmax) = NaN;
corrfem(corrfem>Cmax) = NaN;

%Plots
%colors = jet(NumI);
%colors = winter(NumI);
colors = hot(2*NumI);

figure1 = figure;
p1x = plot(R0',fex','Color','black'); hold on
p1m = plot(R0',fem'); hold on
p2x = plot(R0',fexmean','Color','black'); hold on
p2m = plot(R0',femmean'); hold on
p3x = plot(R0',fexcov','Color','black'); hold on
p3m = plot(R0',femcov'); hold on
p4x = plot(R0',corrfex','Color','black'); hold on
p4m = plot(R0',corrfem'); hold on

 for i =[p1m,p2m,p3m,p4m]
     for j = 1:NumI
         set(i(j),'Color',colors(j,:))
%         set(p1m(j),'DisplayName',num2str(Iprime0(j)))
     end
 end

%create labels
xlabel('Radius of enclosing circle (units of 1/e^{2} beam waist)')
ylabel('Emission intensity distribution')
axis tight
pbaspect([1,1,1])

% Create legend
labels = {'0'};
for x = 1:NumI
    labels(x+1)={num2str(Iprime0(x))};
end
hold off
legend('on')
legend1 = legend(gca,labels);    
set(legend1,'Position',[0.646488294638319,0.586377339597708,0.107023410723361,0.274864368423631],...
    'FontSize',7);

% Create textboxes with annotation
annotation(figure1,'textbox',...
    [0.75 0.81 0.1 0.1],...
    'String',"\it{I_P/I_S}",...
    'FitBoxToText','off',...
    'EdgeColor','none','FontName','Cambria Math');

annotation(figure1,'textbox',...
    [0.37 0.68 0.06 0.23],...
    'String','\it{C_{em}}',...
    'FitBoxToText','off',...
    'EdgeColor','none','FontName','Cambria Math');

annotation(figure1,'textbox',...
    [0.74 0.20 0.08 0.07],...
    'String',{'<\it{f_{em}}>'},...
    'FitBoxToText','off',...
    'EdgeColor','none','FontName','Cambria Math');

annotation(figure1,'textbox',...
    [0.66 0.51 0.06 0.03],...
    'String','\alpha_{em}',...
    'FitBoxToText','off',...
    'EdgeColor','none','FontName','Cambria Math');

annotation(figure1,'textbox',...
    [0.28 0.29 0.06 0.03],...
    'String','\it{f_{em}}',...
    'FitBoxToText','off',...
    'EdgeColor','none','FontName','Cambria Math');
