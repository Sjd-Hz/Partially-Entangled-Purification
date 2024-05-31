clear
clc
close all;

sum=0;
jj=0;
for a=0:0.01:0.5
    jj=jj+1;
    alpha=sqrt(a);
    beta=sqrt(1-alpha^2);
        for N=1:1:8
            P_suc(N)=(2*alpha^(2^N)*beta^(2^N))/(alpha^(2^N)+beta^(2^N))^2;
            mult=1;
                for i=1:1:N-1
                    mult=mult*(alpha^(2^(i+1))+beta^(2^(i+1)))/(alpha^(2^i)+beta^(2^i))^2; %failure of previous rounds
                end
            sum=sum+P_suc(N)*mult; %Sum of net probabilities
            P_sum(N)=sum;
                if N==1
                    val=0.5;
                else
                val=P_sum(N)-P_sum(N-1); %check if the sum of net probabilities reaches a fixed value
                end
                if val<0.03
                    round(jj)=N-1; %Required rounds to reach a fixed value of sum of net probabilities 
                    break
                end        
        end
end
figure(2)
plot(round,':o','LineWidth', 2)
axis tight
grid on
yticks([1 2 3 4])
yticklabels({'1','2','3','4'})
xticks([1 6 11 16 21 26 31 36 41 46 51])
xticklabels({'0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5'})
xlabel('\alpha^2','FontSize',12);
ylabel('Required rounds','FontSize',12);
