% Erick Jorge Canales-Rodriguez, 2024

clear all; clc;

ROIs = {'ROI2', 'ROI5', 'ROI8', 'ROI10'};

figure('units', 'normalized', 'outerposition', [0.1 0.1 0.55 0.8]);
t = tiledlayout(2, 2, "TileSpacing", "tight");
t.TileSpacing = 'compact';
t.Padding = 'compact';

%D = 0.3; % um^2/ms
D = 0.5;
%D = 0.8;

for i=1:4
nexttile
ROI = ROIs{i};

if strcmp (ROI,'ROI2')
    rad   = load('Data/Hist/rad_ROI2_Giorgio.mat');
    if D == 0.3
        aeff = 0.6008016;
    elseif D == 0.5
        aeff = 0.52224449;
    elseif D == 0.8
        aeff = 0.39458918;
    end
    legend_ROI = "CC: Prefontal cortex";
end
if strcmp (ROI,'ROI5')
    rad   = load('Data/Hist/rad_ROI5_Giorgio.mat');
    if D == 0.3
        aeff = 0.9739479;
    elseif D == 0.5
        aeff = 1.00340681;
    elseif D == 0.8
        aeff = 0.98376754;
    end
    legend_ROI = "CC: Motor cortex";
end
if strcmp (ROI,'ROI8')
    rad   = load('Data/Hist/rad_ROI8_Giorgio.mat');
    if D == 0.3
        aeff = 0.64008016;
    elseif D == 0.5
        aeff = 0.57134269;
    elseif D == 0.8
        aeff = 0.45350701;
    end
    legend_ROI = "CC: Parietal cortex";
end
if strcmp (ROI,'ROI10')
    rad   = load('Data/Hist/rad_ROI10_Giorgio.mat');
    if D == 0.3
        aeff = 0.8757515;
    elseif D == 0.5
        aeff = 0.84629259;
    elseif D == 0.8
        aeff = 0.76773547;
    end
    legend_ROI = "CC: Visual cortex";
end
rad   = rad.rad;

%nbins = 33;
%X     = linspace(min(rad)*0.8, max(rad)*1.2, nbins);
X     = 0.05:0.1:5.05; % The center of these bins are located at 0.1:0.1:5.0.

H     = histogram(rad, X, 'Normalization','probability');
Values = H.Values;
save (['Data' filesep 'Hist' filesep 'Histvalues_' ROI], 'Values');  

[muhat,muci] = gamfit(rad);
a = muhat(1);
b = muhat(2);
P_inner_axon = gampdf(X,a,b);
hold on
plot(X, P_inner_axon/sum(P_inner_axon), 'blue', 'LineWidth', 1.5)
%plot(X, P_inner_axon, 'blue', 'LineWidth', 1.5)


k     = 1/b;
alpha = a;

g     = 0.7;
% bad old result (unnormalized) Pm    = ( k/gamma(alpha) ) * ( g/(1-g) ) * ( igamma(alpha-1, X*g*k) - igamma(alpha-1, X*k/g)  );
Pm    = ( k/gamma(alpha) ) * ( g/(1-g^2) ) * ( igamma(alpha-1, X*g*k) - igamma(alpha-1, X*k/g)  );

hold on
plot(X, Pm/sum(Pm), 'color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5)
%plot(X, Pm, 'color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5)


% mean
Pm     = Pm/sum(Pm);
mean   = sum(X.*Pm);
mean21 = sum((X.^2).*Pm)/mean;
mean31 = sqrt(sum((X.^3).*Pm)/mean);

x1 = xline(mean, ':k');

x1.LineWidth = 1.6;
x1.LabelHorizontalAlignment='center';
x1.LabelVerticalAlignment = 'top';
x1.LabelOrientation = 'horizontal';
x1.FontSize                 = 12;
x1.FontAngle                = 'italic';

x2 = xline(mean21, '--k');

x2.LineWidth = 1.6;
x2.LabelHorizontalAlignment = 'center';
x2.LabelOrientation         = 'horizontal';
x2.FontSize                 = 12;
x2.FontAngle                = 'italic';

x3 = xline(mean31, 'k');
x3.LineWidth = 1.6;
x3.LabelHorizontalAlignment = 'center';
x3.LabelVerticalAlignment   = 'top';
x3.LabelOrientation         = 'horizontal';
x3.FontSize                 = 12;
x3.FontAngle                = 'italic';

x4 = xline(aeff, '--r');
x4.LineWidth = 1.6;
x4.LabelHorizontalAlignment = 'center';
x4.LabelVerticalAlignment   = 'middle';
x4.LabelOrientation         = 'horizontal';
x4.FontSize                 = 12;
x4.FontAngle                = 'italic';

xlim([0.1, 3.1])

%lgd = legend('Histogram of axon radius','Gamma distribution of axon radius', 'Distribution of myelin sheath radius', '$<r>$', '$\frac{<r^2>}{<r>}$', '$(\frac{<r^3>}{<r>})^{\frac{1}{2}}$', '$a_{eff}$', 'Interpreter','latex', 'FontSize',11);
%title(lgd,'Legend');
if i==1
    legend('Histogram of axon radius','Gamma distribution of axon radius', 'Distribution of myelin sheath radius', '$<a>$', '$<a^2>/<a>$', '$(<a^3>/<a>)^{1/2}$', '$a_{eff}$', 'Interpreter','latex', 'FontSize',12);
    legend('show');
end

xlabel('$radius \ ( \mu m)$', 'Interpreter','latex', 'FontSize',13);
ylabel('Relative Frequency', 'Interpreter','latex', 'FontSize',13)
title(legend_ROI)

xticks([0.1, 0.5:0.5:3.0])

ax = gca;
ax.XTickMode = 'manual';
ax.YTickMode = 'manual';
ax.ZTickMode = 'manual';
ax.XLimMode  = 'manual';
ax.YLimMode  = 'manual';
ax.ZLimMode  = 'manual';
set(gca, 'DefaultAxesLineWidth', 1.5)

end

exportgraphics(t, ['Fig6_Effective_radius_ROIs_D' num2str(D) '.png'], 'Resolution', 600, 'BackgroundColor','white')
exportgraphics(t, ['Fig6_Effective_radius_ROIs_D' num2str(D) '.pdf'], 'Resolution', 600, 'BackgroundColor','white', 'ContentType','vector')
%exportgraphics(t, ['Fig6_Effective_radius_ROIs_D' num2str(D) '.tiff'],'Resolution', 600, 'BackgroundColor','white', 'ContentType','vector')
