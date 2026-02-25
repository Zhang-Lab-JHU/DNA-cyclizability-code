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


%Least Squares for the Random Library
Tdata = readtable('/Users/rita/Downloads/Random-c0free-readCounts.txt', 'Delimiter', '\t');

Seq = table2array(Tdata(:, "Sequence"));

NN = height(Tdata);
SS = zeros(NN, 200); % Pre-allocate SS matrix

TilingC0 = table2array(Tdata(:, "C0free"));

[TPC0,TC0axis] = hist(TilingC0,50);
TdC0 = mean(diff(TC0axis));
TPC0 = TPC0/(NN*TdC0);

for kk=1:100
    idx = randperm(NN);
    idx = idx(1:round(NN/10));
    TPC0_half(kk,:) = hist(TilingC0(idx),TC0axis)/(length(idx)*TdC0);
end


for n = 1:NN
    seq = char(Seq(n));    
    [~, x] = ismember(seq(1:50), 'ATCG');

    for k = 1:50
        SS(n, (k-1)*4+1:k*4) = x(k) == [1 2 3 4];
    end

end


XX = TilingC0-mean(TilingC0);
ss = SS - ones(NN,1)*mean(SS);


for kk=1:100
    idx = randperm(NN);
    idx = idx(1:round(NN/10));
    triggered_mean(kk,:) = pinv(ss(idx, :))*XX(idx);
    %triggered_mean(kk,:) = XX(idx)'*ss(idx,:)/length(idx); %only RL
end
for kk=1:1000
    xx = XX(randperm(NN));
    idx = randperm(NN);
    idx = idx(1:round(NN/10));
    mm(kk,:) = pinv(ss(idx, :))*xx(idx);
    %mm(kk,:) = xx(idx)'*ss(idx,:)/length(idx); %only RL
end

fig1=figure(1);
set(fig1,'Position',Position1)
scatter([1:200],mean(triggered_mean),15,'filled'); hold on
errorbar(mean(triggered_mean),std(triggered_mean),'.','LineWidth',lw3,'Color',Color1)
plot([1:200],std(mm),'-','LineWidth',lw1,'Color',Color2)
plot([1:200],-std(mm),'-','LineWidth',lw1,'Color',Color2)
hold off
xlabel('i\alpha')
ylabel('mean(C_0{\times}S_i^\alpha)_c')
set(gca,'FontSize',fs,'Box','Off','TickDir','Out','LineWidth',lw3)



%Least Squares for the Tiling Library
Tdata = readtable('/Users/rita/Downloads/Tiling-c0free-readCounts (1).txt', 'Delimiter', '\t');

Seq = table2array(Tdata(:, "Sequence"));

NN = height(Tdata);
SS = zeros(NN, 200); % Pre-allocate SS matrix

TilingC0 = table2array(Tdata(:, "C0free"));

[TPC0,TC0axis] = hist(TilingC0,50);
TdC0 = mean(diff(TC0axis));
TPC0 = TPC0/(NN*TdC0);

for kk=1:100
    idx = randperm(NN);
    idx = idx(1:round(NN/10));
    TPC0_half(kk,:) = hist(TilingC0(idx),TC0axis)/(length(idx)*TdC0);
end


for n = 1:NN
    seq = char(Seq(n));    
    [~, x] = ismember(seq(1:50), 'ATCG');

    for k = 1:50
        SS(n, (k-1)*4+1:k*4) = x(k) == [1 2 3 4];
    end

end


XX = TilingC0-mean(TilingC0);
ss = SS - ones(NN,1)*mean(SS);


for kk=1:100
    idx = randperm(NN);
    idx = idx(1:round(NN/10));
    triggered_mean(kk,:) = pinv(ss(idx, :))*XX(idx);
    %triggered_mean(kk,:) = XX(idx)'*ss(idx,:)/length(idx); %only RL
end
for kk=1:1000
    xx = XX(randperm(NN));
    idx = randperm(NN);
    idx = idx(1:round(NN/10));
    mm(kk,:) = pinv(ss(idx, :))*xx(idx);
    %mm(kk,:) = xx(idx)'*ss(idx,:)/length(idx); %only RL
end

fig2=figure(2);
set(fig2,'Position',Position1)
scatter([1:200],mean(triggered_mean),15,'filled'); hold on
errorbar(mean(triggered_mean),std(triggered_mean),'.','LineWidth',lw3,'Color',Color1)
plot([1:200],std(mm),'-','LineWidth',lw1,'Color',Color2)
plot([1:200],-std(mm),'-','LineWidth',lw1,'Color',Color2)
hold off
xlabel('i\alpha')
ylabel('mean(C_0{\times}S_i^\alpha)_c')
set(gca,'FontSize',fs,'Box','Off','TickDir','Out','LineWidth',lw3)



%correlation function approach to calculating W for the Random Library
Tdata = readtable('/Users/rita/Downloads/Random-c0free-readCounts.txt', 'Delimiter', '\t');

Seq = table2array(Tdata(:, "Sequence"));

NN = height(Tdata);
SS = zeros(NN, 200); % Pre-allocate SS matrix

TilingC0 = table2array(Tdata(:, "C0free"));

[TPC0,TC0axis] = hist(TilingC0,50);
TdC0 = mean(diff(TC0axis));
TPC0 = TPC0/(NN*TdC0);

for kk=1:100
    idx = randperm(NN);
    idx = idx(1:round(NN/10));
    TPC0_half(kk,:) = hist(TilingC0(idx),TC0axis)/(length(idx)*TdC0);
end


for n = 1:NN
    seq = char(Seq(n));    
    [~, x] = ismember(seq(1:50), 'ATCG');

    for k = 1:50
        SS(n, (k-1)*4+1:k*4) = x(k) == [1 2 3 4];
    end

end


XX = TilingC0-mean(TilingC0);
ss = SS - ones(NN,1)*mean(SS);


for kk=1:100
    idx = randperm(NN);
    idx = idx(1:round(NN/10));
    %triggered_mean(kk,:) = pinv(ss(idx, :))*XX(idx);
    triggered_mean(kk,:) = XX(idx)'*ss(idx,:)/length(idx); %only RL
end
for kk=1:1000
    xx = XX(randperm(NN));
    idx = randperm(NN);
    idx = idx(1:round(NN/10));
    %mm(kk,:) = pinv(ss(idx, :))*xx(idx);
    mm(kk,:) = xx(idx)'*ss(idx,:)/length(idx); %only RL
end

fig3=figure(3);
set(fig3,'Position',Position1)
scatter([1:200],mean(triggered_mean),15,'filled'); hold on
errorbar(mean(triggered_mean),std(triggered_mean),'.','LineWidth',lw3,'Color',Color1)
plot([1:200],std(mm),'-','LineWidth',lw1,'Color',Color2)
plot([1:200],-std(mm),'-','LineWidth',lw1,'Color',Color2)
hold off
xlabel('i\alpha')
ylabel('mean(C_0{\times}S_i^\alpha)_c')
set(gca,'FontSize',fs,'Box','Off','TickDir','Out','LineWidth',lw3)


