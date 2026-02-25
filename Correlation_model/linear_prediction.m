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

Tdata = readtable('/Users/rita/Downloads/Tiling-c0free-readCounts (1).txt', 'Delimiter', '\t');
%Tdata = readtable('/Users/rita/Downloads/Random-c0free-readCounts.txt', 'Delimiter', '\t');

Seq = table2array(Tdata(:, "Sequence"));

NN = height(Tdata);
SS = zeros(NN, 200); 

TilingC0 = table2array(Tdata(:, "C0free"));

N1 = round(NN / 10);  
N0 = NN - N1; 



for n = 1:NN
    seq = char(Seq(n));    
    [~, x] = ismember(seq(1:50), 'ATCG');

    for k = 1:50
        SS(n, (k-1)*4+1:k*4) = x(k) == [1 2 3 4];
    end

end

for kk=1:8
    kk
    ss = SS - ones(NN,1)*mean(SS);
    
    idx = randperm(NN);
    idx0 = idx(1:N0); %data
    idx1 = idx(N0+(1:N1)); %test
    
    
    C0data = TilingC0(idx0);
    C0actual = TilingC0(idx1);
    
    XXdata = C0data-mean(C0data);
    ssdata = ss(idx0, :);
    

    %%%% correlation function model
    %cc = zeros(1);
    %for n=1:length(idx0)
    %        cc = cc + XXdata'*ssdata/length(idx0);
    %end

    %C0predicted0 = zeros(N1,1); 
    %for n=1:length(idx1)
    %    C0predicted0(n,1) = cc*ss(idx1(n),:)'+mean(C0data);
    %end
    %%%%


    %%%% least-squares 
    cc = pinv(ssdata)*XXdata;
    C0predicted0 = ss(idx1,:)*cc+mean(C0data);
    %%%%


    Corr0=corrcoef(C0actual,C0predicted0);
    corr_array(kk) = Corr0(1,2);
end


fig2=figure(1);
set(fig2,'Position',Position1)
plot(C0actual,C0predicted0,'.');

hold off
xlabel('i\alpha')
ylabel('mean(C_0{\times}S_i^\alpha)_c')
set(gca,'FontSize',fs,'Box','Off','TickDir','Out','LineWidth',lw3)

mean(corr_array)
std(corr_array)
