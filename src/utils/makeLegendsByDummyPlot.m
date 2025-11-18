function makeLegendsByDummyPlot(labelList,colorList,options)
arguments
    labelList (:,1) string
    colorList (:,3) {mustBeNumeric}
    options.Marker = 'o';
    options.MarkerEdgeColor = 'none';
    options.MarkerSize = 12;
    options.LineStyle = 'none';
    options.FontSize = 12;
end
NList = length(labelList);
hold on
for n = 1:NList
    h(n) = plot(NaN,NaN, ...
        'Marker',options.Marker, ...
        'MarkerFaceColor',colorList(n,:), ...
        'MarkerEdgeColor',options.MarkerEdgeColor, ...
        'MarkerSize',options.MarkerSize, ...
        'LineStyle',options.LineStyle, ...
        'DisplayName',labelList(n));
end
lgd = legend(h,'Location','northwest');
lgd.FontSize = options.FontSize;

axis off
end

