function STATS=wilcoxon(x1,x2,varargin)
%This file execute the non parametric Wilcoxon test to evaluate the difference between paired (dependent) samples. 
% If the number of difference is less than 15, the algorithm calculate the exact ranks distribution; 
% else it uses a normal distribution approximation. 
% Now, the MatLab function SIGNRANK returns the same p-value. 
% Anyway, this Wilcoxon function gives a more detailed output (that is necessary for publications...)
% 
% Syntax: 	STATS=WILCOXON(X1,X2,PLTS)
%      
%     Inputs:
%           X1 and X2 - data vectors.
%           ALPHA - significance level (default = 0.05).
%           PLTS - Flag to set if you don't want (0) or want (1) view the plots
%     Outputs:
%           - W value and p-value when exact ranks distribution is used.
%           - W value, Z value, Standard deviation (Mean=0), p-value when normal distribution is used
%        If STATS nargout was specified the results will be stored in the STATS
%        struct.
% 
%      Example: 
% 
%         X1=[77 79 79 80 80 81 81 81 81 82 82 82 82 83 83 84 84 84 84 85 85 86 86 87 87];
% 
%         X2=[82 82 83 84 84 85 85 86 86 86 86 86 86 86 86 86 87 87 87 88 88 88 89 90 90];
% 
%           Calling on Matlab the function: wilcoxon(X1,X2)
% 
%           Answer is:
% 
% WILCOXON TEST
%  
%                                 Mean_of_differences    Confidence_interval
%                                 ___________________    ___________________
% 
%     Binomial_estimator            3                    3    4             
%     Hodges_Lehmann_estimator    3.5                    3    4             
% 
% Sample size is good enough to use the normal distribution approximation
%  
%      W     Mean      SD        Z       p_value_two_tails
%     ___    ____    ______    ______    _________________
% 
%     325    0       73.161    4.4354    9.1886e-06       
% 
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
% 
% To cite this file, this would be an appropriate format:
% Cardillo G. (2006). Wilcoxon test: non parametric Wilcoxon test for paired samples.
% http://www.mathworks.com/matlabcentral/fileexchange/12702

%Input Error handling
p = inputParser;
addRequired(p,'x1',@(x) validateattributes(x,{'numeric'},{'row','real','finite','nonnan','nonempty'}));
addRequired(p,'x2',@(x) validateattributes(x,{'numeric'},{'row','real','finite','nonnan','nonempty'}));
addOptional(p,'alpha',0.05, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',1}));
addOptional(p,'plts',0, @(x) isnumeric(x) && isreal(x) && isfinite(x) && isscalar(x) && (x==0 || x==1));
parse(p,x1,x2,varargin{:});
alpha=p.Results.alpha; plts=p.Results.plts;
clear p
assert(length(x1)==length(x2),'Warning: x1 and x2 must have the same length')

disp('WILCOXON TEST')
disp(' ')
dff=sort(x2-x1); %difference between x1 and x2
dff(dff==0)=[]; %eliminate null variations
n=length(dff); %number of ranks
if length(x1)~=n %tell me if there are null variations
    fprintf('There are %d null variations that will be deleted\n',length(x1)-n)
end
if isempty(dff) %if all variations are null variations exit function
    disp('There are not variations. Wilcoxon test can''t be performed')
    return       
end

%Ranks of absolute value of samples differences with sign
[r,t]=tiedrank(abs(dff)); %ranks and ties
W=sum(r.*sign(dff)); %Wilcoxon statics (sum of ranks with sign)
pem=median(dff); %point estimation of median of differences
m=ceil(n/2); %location of the median
if mod(n,2)==0 %If the length of the series is even
    tmp=[dff(1:m) pem dff(m:end)]; %add the median in the middle
    dff=tmp; clear tmp 
    m=m+1;
end
%find how many differences far from the median we have to choose
C=cumsum(binopdf(0:1:n,n,0.5));
T=find(C<=alpha/2,1,'last')-1;
cintpem=dff([m-T m+T]); %construct the interval
clear C T m
[I,J]=ndgrid(dff,dff); d=triu(I+J)./2; %Walsh averages triangular matrix
ld=sort(d(d~=0)); %linearization of Walsh averages matrix
clear I J 
HLe=median(ld); %Hodges-Lehmann estimator
if n>15
    A=n*(n+1)/4; B=realsqrt(n*(n+1)*(2*n+1)/24); 
    Za=-realsqrt(2).*erfcinv(2.*(1-alpha/2));
    T=fix(A-Za.*B);
else
    TC=[0 0 0 0 0 0 2 3 5 8 10 13 17 21 25];
    T=TC(n); clear TC
end
cintHLe=ld([T+1 end-T])';
clear dff d ld T%clear unnecessary variable
disp(table([pem;HLe],[cintpem;cintHLe],...
    'VariableNames',{'Mean_of_differences' 'Confidence_interval'},...
    'RowNames',{'Binomial_estimator','Hodges_Lehmann_estimator'}))
%If the number of elements N<15 calculate the exact distribution of the
%signed ranks (the number of combinations is 2^N); else use the normal
%distribution approximation.

if n<=15
    ap=ff2n(n); %the all possibilities based on two-level full-factorial design.
    ap(ap~=1)=-1; %change 0 with -1
    k=1:1:n; 
    J=ap*k'; %all possible sums of ranks for k elements
    %to compute the p-value see how many values are more extreme of the observed
    %W and then divide for the total number of combinations
    p=length(J(abs(J)>=abs(W)))/length(J); %p-value
    %display results
    disp('The exact Wilcoxon distribution was used')
    disp(' ')
    disp(table(W,p,'VariableNames',{'W' 'p_value_two_tails'}))
    if nargout
        STATS.method='Exact distribution';
        STATS.W=W;
        STATS.p=p;
    end
else
    mW=0;
    sW=sqrt((2*n^3+3*n^2+n-t)/6); %standard deviation
    zW=(abs(W)-0.5)/sW; %z-value with correction for continuity
    p=1-normcdf(zW); %p-value
    %display results
    disp('Sample size is good enough to use the normal distribution approximation')
    disp(' ')
    disp(table(W,mW,sW,zW,p,2*p,'VariableNames',{'W' 'Mean' 'SD' 'Z' 'p_value_one_tail' 'p_value_two_tails'}))
    if nargout
        STATS.method='Normal approximation';
        STATS.W=W;
        STATS.mean=0;
        STATS.std=sW;
        STATS.z=zW;
        STATS.p=p;
    end
end

if plts
    scrsz = get(groot,'ScreenSize');
    hfig1=figure; POS=scrsz; POS(3)=POS(3)/2;
    set(hfig1,'Position',POS)
    xg=repmat([1 2],length(x1),1); yg=[x1; x2]';
    plot(xg,yg,'b.',xg',yg','r-')
    axis square
    set(gca,'XLim',[0 3],'XtickMode','manual','Xtick',0:1:3,'XtickLabel',{' ','Before','After',' '})
    title('Wilcoxon''s Plot')
    hfig2=figure; POS(1)=POS(1)+POS(3);
    set(hfig2,'Position',POS)
    boxplot(yg(:),xg(:),'notch','on')
    set(gca,'XtickLabel',{'Before','After'})
end