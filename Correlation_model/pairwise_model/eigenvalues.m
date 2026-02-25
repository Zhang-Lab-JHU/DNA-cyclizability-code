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

Tdata = readtable('/Users/rita/Downloads/Random-c0free-readCounts (1).txt', 'Delimiter', '\t');

NN = height(Tdata);
N1 = round(NN / 10);  
N0 = NN - N1;  

Seq1 = table2array(Tdata(:, "Sequence"));
TilingC0 = table2array(Tdata(:, "C0free"));

for n = 1:NN
    seq = char(Seq1(n));    
    [~, x] = ismember(seq(1:50), 'ATCG');

    for k = 1:50
        SS(n, (k-1)*4+1:k*4) = x(k) == [1 2 3 4];
    end

end



for kk=1:10    
    idx = randperm(NN);
    idx0 = idx(1:N0); %data
    idx1 = idx(N0+(1:N1)); %test
    
    
    XX = TilingC0-mean(TilingC0);
    XXr = XX(randperm(NN));
    
    ss = SS - ones(NN,1)*mean(SS);
    
    
    C0data = TilingC0(idx0);
    C0actual = TilingC0(idx1);
    
    XXdata = C0data-mean(C0data);
    XXdatar = XXr(idx0);

    cc = make_cc(Seq1, TilingC0, idx0, idx1, NN, ss, XXdata);
    eigenvectors(kk,:) = sort(real(eig(cc)),'descend')/std(XXdata);

    ccr = make_cc(Seq1, TilingC0, idx0, idx1, NN, ss, XXdatar);
    shuffled_eigenvectors(kk,:) = sort(real(eig(ccr)),'descend')/std(XXdatar);


end


%eigenvalues
fig5=figure(3);
set(fig5,'Position',Position2)
errorbar([1:200],mean(shuffled_eigenvectors),std(shuffled_eigenvectors),'o','Color',Color2,'LineWidth',lw3,'MarkerSize',ms,'MarkerFaceColor',Color2,'MarkerEdgeColor',Color2); hold on
errorbar([1:200],mean(eigenvectors),std(eigenvectors),'o','Color',Color1,'LineWidth',lw3,'MarkerSize',ms,'MarkerFaceColor',Color1,'MarkerEdgeColor',Color1)
hold off
axis([-4 205 -0.05 0.15])
yticks(-0.05:0.05:0.15)
set(gca,'FontSize',fs,'Box','Off','TickDir','Out','LineWidth',lw3)
legend('shuffled','data')
legend boxoff
xlabel('index i\alpha')
ylabel('eigenvalues of the interaction matrix')

%C0predicted0 = zeros(N1,1); %after imposing translation invariance
%for n=1:length(idx1)
%    C0predicted0(n,1) = 8*ss(idx1(n),:)*cc*ss(idx1(n),:)'+mean(C0data);
%end
%Corr0=corrcoef(C0actual,C0predicted0);
%Corr0(1,2)

%fig1 = figure(1);
%set(fig1, 'Position', Position1)
%fontsize(fig1, 14, 'points');
%plot(C0actual,C0predicted0,'.');
%xlabel('measured C_0')
%ylabel('predicted C_0')

function cc = make_cc(Seq, TilingC0, idx0, idx1, NN, ss, XXdata)
  
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

