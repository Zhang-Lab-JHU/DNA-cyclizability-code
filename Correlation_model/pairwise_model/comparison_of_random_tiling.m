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

Tdata = readtable('/Users/rita/Downloads/Tiling-c0free-readCounts (1).txt', 'Delimiter', '\t');
Seq1 = table2array(Tdata(:, "Sequence"));
TilingC01 = table2array(Tdata(:, "C0free"));

NN1 = height(Tdata);
N11 = round(NN1 / 10);  
N01 = NN1 - N11;  


idx = randperm(NN);
idx0 = idx(1:N0); %data
idx1 = idx(N0+(1:N1)); %test

cc = make_cc(Seq, TilingC0, idx0, idx1, NN);


idx1 = randperm(NN1);
idx01 = idx1(1:N01); %data
idx11 = idx1(N01+(1:N11)); %test


J = make_J_flat(Seq1, TilingC01, idx01, idx11);
J_square = square_J(J);


Corr0=corrcoef(cc, J_square);
Corr0(1,2)


fig1 = figure(7);
set(fig1, 'Position', Position1)
fontsize(fig1, 14, 'points');

Color1 = [47, 156, 91] / 255;

scatter(16*cc(:), J_square(:), 10, 'filled', 'MarkerFaceColor', Color1, 'MarkerEdgeColor', Color1);
xlabel('Tiling library J_{ij}', 'Interpreter', 'tex');
ylabel('Random library J_{ij}', 'Interpreter', 'tex');
xticks([-0.02:0.01:0.02])
yticks([-0.02:0.01:0.02])
axis square;
%axis([-0.025 0.025 -0.025 0.025]);

[MM_ti1,EE_ti1] = eig(J_square);
[ee_ti1,ii_ti1] = sort(real(diag(EE_ti1)),'descend');


%eigenvectors
fig6=figure(8);
set(fig6,'Position',Position1)
subplot(4,1,1) 
imagesc([26:75],[1:4],reshape(MM_ti1(:,ii_ti1(1)),4,50))
caxis([-0.15 0.15])
xticks(25:10:75)
ylabel('base')
title('mode 1')
yticks(1:4)
yticklabels({'A','T','C','G'})
colorbar
set(gca,'FontSize',fs,'Box','Off','TickDir','Out','LineWidth',lw3)

subplot(4,1,2)
imagesc([26:75],[1:4],reshape(MM_ti1(:,ii_ti1(2)),4,50))
colorbar
caxis([-0.15 0.15])
xticks(25:10:75)
ylabel('base')
title('mode 2')
yticks(1:4)
yticklabels({'A','T','C','G'})
set(gca,'FontSize',fs,'Box','Off','TickDir','Out','LineWidth',lw3)

subplot(4,1,3)
imagesc([26:75],[1:4],reshape(MM_ti1(:,ii_ti1(199)),4,50))
colorbar
caxis([-0.15 0.15])
xticks(25:10:75)
ylabel('base')
title('mode 199')
yticks(1:4)
yticklabels({'A','T','C','G'})
set(gca,'FontSize',fs,'Box','Off','TickDir','Out','LineWidth',lw3)

subplot(4,1,4)
imagesc([26:75],[1:4],reshape(MM_ti1(:,ii_ti1(200)),4,50))
caxis([-0.15 0.15])
xticks(25:10:75)
xticklabels({'25','35','45','55','65','75'})
yticks(1:4)
yticklabels({'A','T','C','G'})
colorbar
xlabel('sequence position')
ylabel('base')
title('mode 200')
set(gca,'FontSize',fs,'Box','Off','TickDir','Out','LineWidth',lw3)

tl_least_mode = reshape(MM_ti1(:,ii_ti1(200)),4,50);
[M,I_max] = max(tl_least_mode);
I_max



[MM_ti1,EE_ti1] = eig(16*cc);
[ee_ti1,ii_ti1] = sort(real(diag(EE_ti1)),'descend');


%eigenvectors
fig6=figure(9);
set(fig6,'Position',Position1)
subplot(4,1,1) 
imagesc([26:75],[1:4],reshape(MM_ti1(:,ii_ti1(1)),4,50))
caxis([-0.15 0.15])
xticks(25:10:75)
ylabel('base')
title('mode 1')
yticks(1:4)
yticklabels({'A','T','C','G'})
colorbar
set(gca,'FontSize',fs,'Box','Off','TickDir','Out','LineWidth',lw3)

subplot(4,1,2)
imagesc([26:75],[1:4],reshape(MM_ti1(:,ii_ti1(2)),4,50))
colorbar
caxis([-0.15 0.15])
xticks(25:10:75)
ylabel('base')
title('mode 2')
yticks(1:4)
yticklabels({'A','T','C','G'})
set(gca,'FontSize',fs,'Box','Off','TickDir','Out','LineWidth',lw3)

subplot(4,1,3)
imagesc([26:75],[1:4],reshape(MM_ti1(:,ii_ti1(199)),4,50))
colorbar
caxis([-0.15 0.15])
xticks(25:10:75)
ylabel('base')
title('mode 199')
yticks(1:4)
yticklabels({'A','T','C','G'})
set(gca,'FontSize',fs,'Box','Off','TickDir','Out','LineWidth',lw3)

subplot(4,1,4)
imagesc([26:75],[1:4],reshape(MM_ti1(:,ii_ti1(200)),4,50))
caxis([-0.15 0.15])
xticks(25:10:75)
xticklabels({'25','35','45','55','65','75'})
yticks(1:4)
yticklabels({'A','T','C','G'})
colorbar
xlabel('sequence position')
ylabel('base')
title('mode 200')
set(gca,'FontSize',fs,'Box','Off','TickDir','Out','LineWidth',lw3)

%%%%% PREDICTION
%x_test = make_x_matrix(Seq(idx1));
%C0predicted0 = x_test*J + mean(TilingC0(idx0));
%Corr0=corrcoef(C0actual,C0predicted0);
%Corr0(1,2)

%fig2=figure(1);
%set(fig2,'Position',Position1)
%plot(C0actual,C0predicted0,'.');
%xlabel('measured C_0')
%ylabel('predicted C_0')
%set(gca,'FontSize',fs,'Box','Off','TickDir','Out','LineWidth',lw3)

%%%%% EIGENVECTORS

%[MM_ti1,EE_ti1] = eig(J_square);
%[ee_ti1,ii_ti1] = sort(real(diag(EE_ti1)),'descend');


function x_vector = x_completion(sample, x_vector, i, j)
    % 'A': [1, 0, 0, 0],
    % 'T': [0, 1, 0, 0],
    % 'C': [0, 0, 1, 0],
    % 'G': [0, 0, 0, 1]
    pos_diff = j - i;
    x = x_vector(6*pos_diff+1 : 6*(pos_diff+1));

    a = sample(i);
    b = sample(j);
    
    if (a == 'A' && b == 'T')
        x(1) = x(1) - 1;
        x(2) = x(2) - 1;
        x(3) = x(3) - 1;
    elseif (a == 'T' && b == 'A')
        x(1) = x(1) - 1;
        x(4) = x(4) - 1;
        x(5) = x(5) - 1;
    elseif (a == 'C' && b == 'G')
        x(3) = x(3) - 1;
        x(5) = x(5) - 1;
        x(6) = x(6) - 1;
    elseif (a == 'G' && b == 'C')
        x(2) = x(2) - 1;
        x(4) = x(4) - 1;
        x(6) = x(6) - 1;
    elseif (a == 'A' && b == 'A') || (a == 'T' && b == 'T')
        x(1) = x(1) + 1;
    elseif (a == 'A' && b == 'C') || (a == 'G' && b == 'T')
        x(2) = x(2) + 1;
    elseif (a == 'A' && b == 'G') || (a == 'C' && b == 'T')
        x(3) = x(3) + 1;
    elseif (a == 'T' && b == 'C') || (a == 'G' && b == 'A')
        x(4) = x(4) + 1;
    elseif (a == 'T' && b == 'G') || (a == 'C' && b == 'A')
        x(5) = x(5) + 1;
    elseif (a == 'C' && b == 'C') || (a == 'G' && b == 'G')
        x(6) = x(6) + 1;
    end
    
    x_vector(6*pos_diff+1 : 6*(pos_diff+1)) = x;
end


function x_matrix = make_x_matrix(ssdata)
    x_matrix = zeros(size(ssdata, 1), 50*6);
    for sample_ind = 1:size(ssdata, 1)
        sample = ssdata(sample_ind, :);
        sample = sample{1};
        x_vector = x_matrix(sample_ind, :);
        for i = 1:49
            for j = i:50
                x_vector = x_completion(sample, x_vector, i, j);
            end
        end
        x_matrix(sample_ind, :) = x_vector;
    end
end



function J = J_completion(J, J_flat, i, j)

    pos_diff = floor((j-1)/4) - floor((i-1)/4);
    
    a = mod(i-1, 4);
    b = mod(j-1, 4);
    
    if (a == 0 && b == 1)
        J(i, j) = -(J_flat(6*pos_diff + 1) + J_flat(6*pos_diff + 2) + J_flat(6*pos_diff + 3));
        J(j, i) = J(i, j);
    elseif (a == 1 && b == 0)
        J(i, j) = -(J_flat(6*pos_diff + 1) + J_flat(6*pos_diff + 4) + J_flat(6*pos_diff + 5));
        J(j, i) = J(i, j);
    elseif (a == 2 && b == 3)
        J(i, j) = -(J_flat(6*pos_diff + 3) + J_flat(6*pos_diff + 5) + J_flat(6*pos_diff + 6));
        J(j, i) = J(i, j);
    elseif (a == 3 && b == 2)
        J(i, j) = -(J_flat(6*pos_diff + 2) + J_flat(6*pos_diff + 4) + J_flat(6*pos_diff + 6));
        J(j, i) = J(i, j);
    elseif (a == 0 && b == 0) || (a == 1 && b == 1)
        J(i, j) = J_flat(6*pos_diff + 1);
        J(j, i) = J(i, j);
    elseif (a == 0 && b == 2) || (a == 3 && b == 1)
        J(i, j) = J_flat(6*pos_diff + 2);
        J(j, i) = J(i, j);
    elseif (a == 0 && b == 3) || (a == 2 && b == 1)
        J(i, j) = J_flat(6*pos_diff + 3);
        J(j, i) = J(i, j);
    elseif (a == 1 && b == 2) || (a == 3 && b == 0)
        J(i, j) = J_flat(6*pos_diff + 4);
        J(j, i) = J(i, j);
    elseif (a == 1 && b == 3) || (a == 2 && b == 0)
        J(i, j) = J_flat(6*pos_diff + 5);
        J(j, i) = J(i, j);
    elseif (a == 2 && b == 2) || (a == 3 && b == 3)
        J(i, j) = J_flat(6*pos_diff + 6);
        J(j, i) = J(i, j);
    end
end

function J = make_J_flat(Seq, TilingC0, idx0, idx1)

    C0data = TilingC0(idx0);
    C0actual = TilingC0(idx1);
    
    XXdata = TilingC0(idx0) - mean(TilingC0(idx0));
    x_matrix = make_x_matrix(Seq(idx0));
    
    J = pinv(x_matrix)*XXdata;
end

function J_square = square_J(J)
    J_square = zeros(200, 200);
    for i = 1:200 
        for j = i:200  
            J_square = J_completion(J_square, J, i, j);
        end
    end
end



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
