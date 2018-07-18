path = '/Users/donguk.kim/projects/dukim_work/CHiMaD3_dendritic_growth/contour_line_sort/';
fname= 'phase_field_1500.csv';
%fname= 'CHiMaD3_DK_contour_test.csv';
numrow = 1305;
numcol = 2;
%numcol = 6;
colIDforx = 1;
colIDfory = 2;

path_fname = [path fname];

dataval = dlmread(path_fname, ',' ,[1 0 numrow-1 numcol-1]);

xdat = dataval(1:numrow-1,colIDforx);
ydat = dataval(1:numrow-1,colIDfory);

plot(xdat,ydat, '-k', 'linewidth',1.5);
xlabel('x', 'fontsize', 22);
ylabel('y', 'fontsize', 22);
axis equal;
set(gca,'fontsize', 22);
xlim([0 960]);
ylim([0 960]);