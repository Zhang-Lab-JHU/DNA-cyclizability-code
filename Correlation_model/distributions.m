clear all
clc

lw1=1.5;
lw2=1;
lw3=0.75;
fs=15;
ms=4;
Position1=[100 100 600 450];
Position2=[100 100 600 420];
Position3=[100 100 600 500];
Position4=[100 100 600 667];
Position5=[100 100 680 680];

Color1=[47,156,91]/255;
Color2=[200,100,75]/255;


Rdata = readtable('/Users/rita/Downloads/Random-c0free-readCounts (1).txt', 'Delimiter', '\t');

RandomC0 = table2array(Rdata(:, "C0free"));
[RPC0,RC0axis] = hist(RandomC0,50);
RdC0 = mean(diff(RC0axis));


NN = height(Rdata);

RPC0 = RPC0/(NN*RdC0);

for kk=1:100
    idx = randperm(NN);
    idx = idx(1:round(NN/2));
    RC0_half(kk,:) = hist(RandomC0(idx),RC0axis)/(length(idx)*RdC0);
end



Tdata = readtable('/Users/rita/Downloads/Tiling-c0free-readCounts (1).txt', 'Delimiter', '\t');

TilingC0 = table2array(Tdata(:, "C0free"));
[TPC0,TC0axis] = hist(TilingC0,50);
TdC0 = mean(diff(TC0axis));
NN = height(Tdata);

TPC0 = TPC0/(NN*TdC0);

for kk=1:100
    idx = randperm(NN);
    idx = idx(1:round(NN/2));
    TPC0_half(kk,:) = hist(TilingC0(idx),TC0axis)/(length(idx)*TdC0);
end

fig1 = figure(1);
set(fig1, 'Position', Position1)

% Plot the random library data
errorbar(RC0axis, mean(RC0_half), std(RC0_half), 'o-', 'LineWidth', lw1, 'Color', Color1, 'MarkerFaceColor', Color1, 'MarkerSize', 5);
hold on

% Plot the tiling library data
errorbar(TC0axis, mean(TPC0_half), std(TPC0_half), 'o-', 'LineWidth', lw2, 'Color', Color2, 'MarkerFaceColor', Color2, 'MarkerSize', 5);

xlabel('C_0')
ylabel('probability density')
axis([-1.5 2.25 0 2])
xticks(-1:2)
yticks(0:0.5:2)
set(gca, 'FontSize', fs, 'Box', 'Off', 'TickDir', 'Out', 'LineWidth', lw3)

legend('random library', 'tiling library', 'Location', 'best', 'Box', 'off')
hold off


