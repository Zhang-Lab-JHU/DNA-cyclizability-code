%parameter fitting on all the sequences

clear all
clc

lw1=1.5;
lw2=1;
lw3=0.75;
fs=15;
ms=4;
Position1=[100 100 600 450];
Color1=[47,156,91]/255;

%Random Library

Tdata = readtable('/Users/rita/Downloads/Random-c0free-readCounts (1).txt', 'Delimiter', '\t');
Seq1 = table2array(Tdata(:, "Sequence"));
TilingC01 = table2array(Tdata(:, "C0free"));

NN1 = height(Tdata);

J_RL = make_J_flat(Seq1, TilingC01);
J_square_RL = square_J(J_RL);

cc = make_cc(Seq1, TilingC01, NN1);

Corr_LS=corrcoef(16*cc, J_square_RL);
Corr_LS(1,2)

%eigenvalue decomposition

[MM_ti1,EE_ti1] = eig(J_square_RL);
[ee_ti1,ii_ti1] = sort(real(diag(EE_ti1)),'descend');

dna_map = ['A' 'T' 'C' 'G'];

rl_1st_mode = reshape(MM_ti1(:,ii_ti1(1)),4,50);
[M,I_max] = max(rl_1st_mode);

rl_1st_sequence1 = dna_map(I_max);
rl_1st_sequence1

rl_1st_mode = reshape(-MM_ti1(:,ii_ti1(1)),4,50);
[M,I_max] = max(rl_1st_mode);

rl_1st_sequence2 = dna_map(I_max);
rl_1st_sequence2

rl_200th_mode = reshape(MM_ti1(:,ii_ti1(200)),4,50);
[M,I_max] = max(rl_200th_mode);

rl_200th_sequence1 = dna_map(I_max);
rl_200th_sequence1

rl_200th_mode = reshape(-MM_ti1(:,ii_ti1(200)),4,50);
[M,I_max] = max(rl_200th_mode);

rl_200th_sequence2 = dna_map(I_max);
rl_200th_sequence2

%Tiling Library

Tdata = readtable('/Users/rita/Downloads/Tiling-c0free-readCounts (1).txt', 'Delimiter', '\t');
Seq1 = table2array(Tdata(:, "Sequence"));
TilingC01 = table2array(Tdata(:, "C0free"));

NN1 = height(Tdata);

J = make_J_flat(Seq1, TilingC01);
J_square = square_J(J);

Corr_RL_TL=corrcoef(J_square_RL, J_square);
Corr_RL_TL(1,2)

%eigenvalue decomposition

[MM_ti1,EE_ti1] = eig(J_square);
[ee_ti1,ii_ti1] = sort(real(diag(EE_ti1)),'descend');

dna_map = ['A' 'T' 'C' 'G'];

tl_1st_mode = reshape(MM_ti1(:,ii_ti1(1)),4,50);
[M,I_max] = max(tl_1st_mode);

tl_1st_sequence1 = dna_map(I_max);
tl_1st_sequence1

tl_1st_mode = reshape(-MM_ti1(:,ii_ti1(1)),4,50);
[M,I_max] = max(tl_1st_mode);

tl_1st_sequence2 = dna_map(I_max);
tl_1st_sequence2

tl_200th_mode = reshape(MM_ti1(:,ii_ti1(200)),4,50);
[M,I_max] = max(tl_200th_mode);

tl_200th_sequence1 = dna_map(I_max);
tl_200th_sequence1

tl_200th_mode = reshape(-MM_ti1(:,ii_ti1(200)),4,50);
[M,I_max] = max(tl_200th_mode);

tl_200th_sequence2 = dna_map(I_max);
tl_200th_sequence2


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


function x_matrix = make_x_matrix(ssdata)

    x_matrix = zeros(size(ssdata, 1), 49*6);

    for sample_ind = 1:size(ssdata, 1)
        sample = ssdata(sample_ind, :);
        sample = sample{1};
        x_vector = x_matrix(sample_ind, :);
        for i = 1:49
            for j = (i+1):50
                x_vector = x_completion(sample, x_vector, i, j);
            end
        end
        x_matrix(sample_ind, :) = x_vector;
    end
end


function J = make_J_flat(Seq, TilingC0)
    XXdata = TilingC0 - mean(TilingC0);
    x_matrix = make_x_matrix(Seq);
    J = pinv((x_matrix-mean(x_matrix)))*XXdata;
end


function J_square = square_J(J)
    J_square = zeros(200, 200);
    for i = 1:199
        for j = (i+1):200  
            J_square = J_completion(J_square, J, i, j);
            
        end
    end
end

function x_vector = x_completion(sample, x_vector, i, j)

    pos_diff = j - i - 1;
    x = x_vector(6*pos_diff + 1: 6*(pos_diff+1));
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
    x_vector(6*pos_diff + 1 : 6*(pos_diff+1)) = x;
end

function J = J_completion(J, J_flat, i, j)
    if i == j
        J(i, j) = 0;
        J(j, i) = 0;
        return;
    end

    seq_i = floor((i-1)/4);
    seq_j = floor((j-1)/4);
    
    if seq_j > seq_i
        pos_diff = seq_j - seq_i - 1;

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


function cc = make_cc(Seq, TilingC0, NN)

    for n = 1:NN
        seq = char(Seq(n));    
        [~, x] = ismember(seq(1:50), 'ATCG');
    
        for k = 1:50
            SS(n, (k-1)*4+1:k*4) = x(k) == [1 2 3 4];
        end
    
    end

    XXdata = TilingC0-mean(TilingC0);
    ss = SS - ones(NN,1)*mean(SS);
    
    cc = zeros(200);
    for n=1:NN
            cc = cc + XXdata(n)*ss(n,:)'*ss(n,:)/NN;
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