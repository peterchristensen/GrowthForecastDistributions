%% Import the data
[~, ~, raw] = xlsread('~\ExpertForecastData.xlsx','Sheet1','A2:H193');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[1,2]);
raw = raw(:,[3,4,5,6,7,8]);

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Allocate imported array to column variable names
Expert = cellVectors(:,1);
Region = cellVectors(:,2);
P_10_2050 = data(:,1);
P_25_2050 = data(:,2);
P_50_2050 = data(:,3);
P_75_2050 = data(:,4);
P_90_2050 = data(:,5);
period = data(:,6);

%% Clear temporary variables
clearvars data raw cellVectors;

%%Generate output variables mu and sigma
mu_out = zeros(192,1);
sigma_out = zeros(192,1);
SumSquares_out = zeros(192,1);
output = zeros(192,3);

% Loop through 193 sets of quantile estimates.
%%% mu is starting value and is set to median estimate.  
%%% mu_ and sigma_ are ranges of possible values.
%%% P is quantile values
%%% P is set of quantile estimates

for j = 1:192
    mu = P_50_2050(j)
    mu_lb = mu-1.5
    mu_ub = mu+1.5
    mu_ = mu_lb:.01:mu_ub
    range_7525 = P_75_2050(j)-P_25_2050(j)
    sigma_lb = 0
    sigma_ub = range_7525+1
    sigma_ = sigma_lb:.01:sigma_ub
    P = [.1 .25 .5 .75 .9]
    q_hat = [P_10_2050(j) P_25_2050(j) P_50_2050(j) P_75_2050(j) P_90_2050(j)]
    
   
%% Compute sum of squared differences for each possible value of sigma with starting value for mu.  Find min value. 

    y = zeros(1,length(sigma_))
    for i = 1:length(sigma_)
        y(i) = sum((norminv(P,mu,sigma_(i))-q_hat).^2)
    end

    A = [sigma_.' y.']

    [M,I] = min(A(:,2))
    sigma_min = A(I,1)

%% Compute sum of squared differences for each possible value of mu.  Find min value. 
 
    x = zeros(1,length(mu_))
    for i = 1:length(mu_)
        x(i) = sum((norminv(P,mu_(i),sigma_min)-q_hat).^2)
    end

    A_1 = [mu_.' x.']
    [M,I] =min(A_1(:,2))
    mu_min = A_1(I,1)

%% Compute sum of squared differences for each possible value of sigma with updated mu.  Find min value. 
  
    z = zeros(1,length(sigma_))
    for i = 1:length(sigma_)
        z(i) = sum((norminv(P,mu_min,sigma_(i))-q_hat).^2)
    end

    A_2 = [sigma_.' z.']

    [M,I] =min(A_2(:,2))
    sigma_min2 = A_2(I,1)

%% Test whether estimated mu is within defined range. If so, export mu.
    
    if mu_min==mu_lb
        mu_out(j) = 111
    elseif mu_min==mu_ub
        mu_out(j) = 555
    else
        mu_out(j) = mu_min
    end  

%% Test whether both values of sigma are equal for starting and estimated values of mu. If so, export sigma.

    if sigma_min==sigma_min2
        sigma_out(j) = sigma_min
    else 
        sigma_out(j) = 999
    end

%% Test whether estimated sigma is within defined range. If so, export sigma.

    if sigma_min == sigma_lb
        sigma_out(j) = 111
    elseif sigma_min == sigma_ub
        sigma_out(j) = 555
    end

%% Export sum of squares and gather output.

    SumSquares_out(j) = A_2(I,2)

end


%% Gather and save output.

output = [period, P_10_2050, P_25_2050, P_50_2050, P_75_2050, P_90_2050, mu_out, sigma_out, SumSquares_out]
datacell = num2cell(output)
xlswrite('~\ReplicationFiles\Data\estimatedparameters.xlsx', [Expert, Region, datacell])
save('~\ReplicationFiles\Data')


%%% Robust Results

h=figure
hold off;
plot(Y,newpdf(14,1:length(Y)),'-r','LineWidth',2)
ax = gca
ax.XTick = [-.04 -.03 -.02 -.01 0 .01 .02 .03 .04 .05 .06];
axis([-.04 .06 0 40]) 
grid on
hold on;
plot(Y,newpdf(15,1:length(Y)),'--r','LineWidth',2)
plot(Y,newpdf(16,1:length(Y)),'-bl','LineWidth',1.5)
plot(Y,newpdf(17,1:length(Y)),'--bl','LineWidth',1.5)
plot(Y,newpdf(18,1:length(Y)),'-.bl','LineWidth',1.5)
legend('mean(\theta)','trmean(\theta)','mean(\mu,\sigma)','median(\mu,\sigma)','trmean(\mu,\sigma)') 
set(legend,'FontSize',12)
ylabel('probability')
xlabel('growth rate')
print(h,'C:\Users\pchrist\Dropbox\TFP_Paper\Figures\2100norm_pdf_robust','-dpdf', '-r0')
