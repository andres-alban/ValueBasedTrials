function [ rval ] = UtilStdizeFigureAbs( fignum )
% UtilStdizeFigure standardizes the size of a given figure in order to help
% make printing from figures more standardized from one figure to the next.
% 
% Revision: AA 6-11-2017

fracheight = 0.95;

scnsize = [1 1 1280 720];
h = figure(fignum);
%axhand = axes();

position = get(h,'Position');
outerpos = get(h,'OuterPosition');
borders = outerpos - position;

axis('square');
tmp = min(scnsize(3),scnsize(4));
tmp2 = fracheight*tmp;
pos1 = [scnsize(3) * (1 - fracheight)/2,...
    scnsize(4) * (1 - fracheight)/2,...
    scnsize(3) * (1 - fracheight)/2 + tmp2,...
    scnsize(4) * (1 - fracheight)/2 + tmp2] ;


set(h,'OuterPosition',pos1)
ax = get(h,'CurrentAxes');
set(ax,'fontsize',14)
lgd = findall(h,'type','legend');
if ~isempty(lgd)
    lgd.FontSize = 14;
end

rval = true;

end

