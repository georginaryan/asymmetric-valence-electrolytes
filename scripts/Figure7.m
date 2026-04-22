close all
clear
% Define parameters
V  = 1;    % adjust as needed
r = 2;    % set zp
ymax=1.5;

% Define x and y ranges
x = linspace(-(1), (1), 800);  % stay within domain of log
y = linspace(-0.1, ymax, 800);
[J, II] = meshgrid(x, y);

% Compute the upper bound: y < xV / log(...)
denominator = log((1 + J) ./ (1 - J));
upperBound = 2*J .* V./ denominator;


% Inequality conditions

% Inequality conditions
cond1 = II > 0;
cond2 = II < upperBound; % region less than transition current (Eqn 79)
cond3 = II > -2*J/r; %Eqn 77 (physical constraint)
cond4 = II >  2*J; %Eqn 77 (physical constraint)
cond5 = J <= 1; % Eqn 77 (limiting current)
cond6 = J >= -1; %Eqn 77 (limiting current)

% near-Gouy-Chapman Region
regionGC = cond1 & cond2 &cond3 &cond4;

% near-Limiting Current Region
regionLC = cond1 & cond3 & cond4 & cond5 & cond6; 

% Plot the shaded region
figure();
hold on

colors = [
    44/255, 123/255, 182/255;  % "#2C7BB6" - dark blue
    116/255, 173/255, 209/255; % "#74ADD1" - medium blue
    171/255, 217/255, 233/255; % "#ABD9E9" - light blue
    189/255, 189/255, 189/255; % "#BDBDBD" - light grey centre
    215/255, 181/255, 216/255; % "#D7B5D8" - light purple
    175/255, 141/255, 195/255; % "#AF8DC3" - medium purple
    118/255, 42/255, 131/255   % "#762A83" - dark purple
];
h2 = imagesc(x, y, regionLC);
set(h2, 'AlphaData', regionLC * 0.2);  % light blue for LC region
h1 = imagesc(x, y, regionGC);
set(h1, 'AlphaData', regionGC * 0.8);
colormap(gca, [1 1 1; colors(2,:)]); % medium blue for GC region

set(gca, 'YDir', 'normal');

% ==== boundaries =====
x_line = linspace(-(1+1), (1+1), 1000);
% II = -2J/r (Eqn 77)
plot(x_line, -2/r*x_line, 'LineWidth', 3, 'Color', colors(6,:)); 
% I=2J (Eqn 77)
plot(x_line, 2*x_line, 'LineWidth', 3, 'Color', colors(6,:));
% II = 0
plot(x_line, zeros(size(x_line)), '--', 'LineWidth', 3, 'Color', colors(4,:)); 
% II = V
plot(x_line, V*ones(size(x_line)), '--', 'LineWidth', 3, 'Color', colors(4,:)); 
% II = JV / log((1+J)/(1-J)) (Eqn 79)
safe_x = x_line(abs(x_line) < 1-1e-4); % stay within log domain
y_curve = (2*safe_x .* V) ./ log((1 + safe_x) ./ (1 - safe_x));
plot(safe_x, y_curve, 'LineWidth', 3, 'Color', colors(1,:));
% Intersections
xline(tanh(V/2), '--','Color', colors(4,:), 'LineWidth', 3); 
xline(-tanh(r*V/2), '--','Color', colors(4,:), 'LineWidth', 3); 
xline(1, 'Color', colors(7,:), 'LineWidth', 3);
xline(-1, 'Color', colors(7,:), 'LineWidth', 3); 

%% == Plotting Settings ==
axis tight;
grid off;
ax = gca;
% Margins
margin = 0.03;
ax.Position = [margin, margin, 1-2*margin, 1-2*margin];
% Domain/Range
ylim([-0.1 ymax]);
xlim([-(1.2) (1.2)]);
x_limits = xlim;
y_limits = ylim;
% Ticks
ax.XTickLabel = [];
ax.YTickLabel = [];
ax.XTick = [];
ax.YTick = [];
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
%Draw a rectangle for the frame to match Mathematica plots
rectangle('Position', [x_limits(1), y_limits(1), ...
                      x_limits(2)-x_limits(1), y_limits(2)-y_limits(1)], ...
          'LineWidth', 2, ...
          'EdgeColor', 'black', ...
          'LineStyle', '-', ...
          'FaceColor', 'none');
uistack(gca, 'top');

% Image export
exportgraphics(gcf, ...
    fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figures', 'Figure7.pdf'), ...
    'ContentType', 'vector', 'Resolution', 600);
%Figure annotations made in Canva software