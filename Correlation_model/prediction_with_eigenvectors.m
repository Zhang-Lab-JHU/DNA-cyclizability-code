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
    
for kk=1:100
    kk    

    idx = randperm(NN);
    idx0 = idx(1:N0); %data
    idx1 = idx(N0+(1:N1)); %test
    
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
    
    C0predicted0 = zeros(N1,1);
    for n=1:length(idx1)
        C0predicted0(n,1) = 8*ss(idx1(n),:)*cc*ss(idx1(n),:)'+mean(C0data);
    end
    Corr0=corrcoef(C0actual,C0predicted0);
    corr_coef(kk) = Corr0(1,2);
    
    
    % first two eigenvectors
    
    [MM_ti1,EE_ti1] = eig(16*cc);
    [ee_ti1,ii_ti1] = sort(real(diag(EE_ti1)),'descend');
    
    MM_ti1_sorted = MM_ti1(:, ii_ti1);    
    V1 = MM_ti1_sorted(:, 1);
    V2 = MM_ti1_sorted(:, 2);
    E1 = ee_ti1(1);
    E2 = ee_ti1(2);

    C0predicted = zeros(length(idx1), 1);
    for n=1:length(idx1)
        s = ss(idx1(n),:);
        C0predicted(n,1) = 8*s*(V1*E1*V1' + V2*E2*V2')*s'+mean(C0data);
    end

    Corr0=corrcoef(C0actual,C0predicted);
    corr_coef_eig(kk) = Corr0(1,2);

end

mean(corr_coef)
std(corr_coef)
mean(corr_coef_eig)
std(corr_coef_eig)



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
