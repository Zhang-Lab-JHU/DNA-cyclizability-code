clear all
clc

lw1=1.5;
lw2=1;
lw3=0.75;
fs=15;
ms=4;
Position1=[100 100 600 450];
Color1=[47,156,91]/255;

Tdata = readtable('/Users/rita/Downloads/Random-c0free-readCounts (1).txt', 'Delimiter', '\t');

Seq = table2array(Tdata(:, "Sequence"));
TilingC0 = table2array(Tdata(:, "C0free"));

NN = height(Tdata);
N1 = round(NN / 10);  
N0 = NN - N1;  

idx = randperm(NN);
idx0 = idx(1:N0); %data
idx1 = idx(N0+(1:N1)); %test

cc = make_cc(Seq, TilingC0, idx0, idx1, NN);



[MM1,EE1] = eig(cc);
[ee1,ii1] = sort(real(diag(EE1)),'descend');


for n = 1:NN
    seq = char(Seq(n));    
    [~, x] = ismember(seq(1:50), 'ATCG');

    for k = 1:50
        SS(n, (k-1)*4+1:k*4) = x(k) == [1 2 3 4];
    end

end

half1 = idx0;
half2 = idx1;

S1_half2 = SS(half2,:)*MM1(:,ii1(1));
S2_half2 = SS(half2,:)*MM1(:,ii1(2));
S199_half2 = SS(half2,:)*MM1(:,ii1(199));
S200_half2 = SS(half2,:)*MM1(:,ii1(200));
S1_half2 = S1_half2/std(S1_half2);
S2_half2 = S2_half2/std(S2_half2);
S199_half2 = S199_half2/std(S199_half2);
S200_half2 = S200_half2/std(S200_half2);

C0 = TilingC0;
AA1_half = polyfit(S1_half2,C0(half2),2);
AA2_half = polyfit(S2_half2,C0(half2),2);
AA199_half = polyfit(S199_half2,C0(half2),2);
AA200_half = polyfit(S200_half2,C0(half2),2);

QS1 = round(S1_half2*5);
QS2 = round(S2_half2*5);
QS199 = round(S199_half2*5);
QS200 = round(S200_half2*5);
q1range = [min(QS1):max(QS1)];
q2range = [min(QS2):max(QS2)];
q199range = [min(QS199):max(QS199)];
q200range = [min(QS200):max(QS200)];
for n=1:length(q1range)
    ii = find(QS1==q1range(n));
    s1axis(n) = mean(S1_half2(ii));
    meanC0_1_half2(n) = mean(C0(half2(ii)));
    stdC0_1_half2(n) = std(C0(half2(ii)));
end
for n=1:length(q2range)
    ii = find(QS2==q2range(n));
    s2axis(n) = mean(S2_half2(ii));
    meanC0_2_half2(n) = mean(C0(half2(ii)));
    stdC0_2_half2(n) = std(C0(half2(ii)));
end
for n=1:length(q199range)
    ii = find(QS199==q199range(n));
    s199axis(n) = mean(S199_half2(ii));
    meanC0_199_half2(n) = mean(C0(half2(ii)));
    stdC0_199_half2(n) = std(C0(half2(ii)));
end
for n=1:length(q200range)
    ii = find(QS200==q200range(n));
    s200axis(n) = mean(S200_half2(ii));
    meanC0_200_half2(n) = mean(C0(half2(ii)));
    stdC0_200_half2(n) = std(C0(half2(ii)));
end

fig1=figure(2);
set(fig1,'Position',Position1)
subplot(2,2,1)
errorbar(s1axis,meanC0_1_half2,stdC0_1_half2,'o','Color',Color1,'LineWidth',lw3,'MarkerSize',ms,'MarkerFaceColor',Color1,'MarkerEdgeColor',Color1)
hold on
plot([-5:0.01:5],AA1_half(1)*([-5:0.01:5].^2)+AA1_half(2)*[-5:0.01:5]+AA1_half(3),'k-','LineWidth',lw1)
hold off
xlabel('sequence feature 1')
ylabel('C_0')
axis([-4 4 -1 2.2])
legend('data','fit')
legend boxoff
set(gca,'FontSize',fs,'Box','Off','TickDir','Out','LineWidth',lw3)
subplot(2,2,2)
errorbar(s2axis,meanC0_2_half2,stdC0_2_half2,'o','Color',Color1,'LineWidth',lw3,'MarkerSize',ms,'MarkerFaceColor',Color1,'MarkerEdgeColor',Color1)
hold on
plot([-5:0.01:5],AA2_half(1)*([-5:0.01:5].^2)+AA2_half(2)*[-5:0.01:5]+AA2_half(3),'k-','LineWidth',lw1)
hold off
xlabel('sequence feature 2')
ylabel('C_0')
axis([-4 4 -1 2.2])
set(gca,'FontSize',fs,'Box','Off','TickDir','Out','LineWidth',lw3)
subplot(2,2,3)
errorbar(s199axis,meanC0_199_half2,stdC0_199_half2,'o','Color',Color1,'LineWidth',lw3,'MarkerSize',ms,'MarkerFaceColor',Color1,'MarkerEdgeColor',Color1)
hold on
plot([-5:0.01:5],AA199_half(1)*([-5:0.01:5].^2)+AA199_half(2)*[-5:0.01:5]+AA199_half(3),'k-','LineWidth',lw1)
hold off
xlabel('sequence feature 199')
ylabel('C_0')
axis([-4 4 -1 0.7])
set(gca,'FontSize',fs,'Box','Off','TickDir','Out','LineWidth',lw3)
subplot(2,2,4)
errorbar(s200axis,meanC0_200_half2,stdC0_200_half2,'o','Color',Color1,'LineWidth',lw3,'MarkerSize',ms,'MarkerFaceColor',Color1,'MarkerEdgeColor',Color1)
hold on
plot([-5:0.01:5],AA200_half(1)*([-5:0.01:5].^2)+AA200_half(2)*[-5:0.01:5]+AA200_half(3),'k-','LineWidth',lw1)
hold off
xlabel('sequence feature 200')
ylabel('C_0')
axis([-4 4 -1 0.7])
set(gca,'FontSize',fs,'Box','Off','TickDir','Out','LineWidth',lw3)





function cc = make_cc(Seq, TilingC0, idx0, idx1, NN)

    for n = 1:NN
        seq = char(Seq(n));    
        [~, x] = ismember(seq(1:50), 'ATCG');
    
        for k = 1:50
            SS(n, (k-1)*4+1:k*4) = x(k) == [1 2 3 4];
        end
    
    end

    XX = TilingC0-mean(TilingC0);
    ss = SS - ones(NN,1)*mean(SS);
    
    C0data = TilingC0(idx0);
    C0actual = TilingC0(idx1);
    
    XXdata = C0data-mean(C0data);
    
    cc = zeros(200);
    for n=1:length(idx0)
            cc = cc + XXdata(n)*ss(idx0(n),:)'*ss(idx0(n),:)/length(idx0);
    end
    for d=0:49
        CC = zeros(4);
        for i=1:(50-d)
            CC = CC+cc(4*(i-1)+(1:4),4*(i-1+d)+(1:4))/(50-d);
        end
        CC = RCInvariance(CC);
        for i=1:(50-d)
            cc(4*(i-1)+(1:4),4*(i-1+d)+(1:4)) = CC;
            cc(4*(i-1+d)+(1:4),4*(i-1)+(1:4)) = CC';
        end
    end

end

function CC=RCInvariance(CC)
CC(1,3)=(CC(1,3)+CC(4,2))/2;
CC(4,2)=CC(1,3);
CC(1,4)=(CC(1,4)+CC(3,2))/2;
CC(3,2)=CC(1,4);    
CC(2,3)=(CC(2,3)+CC(4,1))/2;
CC(4,1)=CC(2,3); 
CC(2,4)=(CC(2,4)+CC(3,1))/2;
CC(3,1)=CC(2,4);
end
