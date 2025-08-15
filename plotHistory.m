function plotHistory(compliance,volumefrac)
iter = numel(compliance);
color_left = 'b';
color_right = 'r';
figure;
yyaxis left;
plot(1:iter, compliance, color_left, 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Compliance', 'FontSize', 12,'Color', color_left); 
ax = gca;
ax.YAxis(1).Color = color_left;  
ax.YAxis(1).FontSize = 10;
grid on;
yyaxis right;
plot(1:iter, volumefrac, color_right, 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Volume Fraction', 'FontSize', 12,'Color', color_right);
ax = gca;
ax.YAxis(2).Color = color_right;
ax.YAxis(2).FontSize = 10;
ylim([0 1]);
yticks(0:0.2:1);
xlabel('Optimization Step', 'FontSize', 12);
legend('Compliance', 'Volume Fraction', 'Location', 'northeast');
set(gca, 'FontSize', 10, 'Box', 'on');
exportgraphics(gcf, 'MBB_history.png', 'Resolution', 1000);
end