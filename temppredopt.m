% Author: Amanda Laubmeier, 2022
% code solves optimization problems 
% to find best pred community for pest control
% using field temperatures and future temps

% For Laubmeier, Tabassum, Tenhumberg 2022

% Requires data from Curtsdotter 2020 to run

function temppredopt
close all,
clear all,

%parameter values from Wootton 2022
a0=9.56; phi=1.34; h0=0.03; b0=10.83; x0=4*exp(19.75);
Ropts=[1, 200, 200, 200, 200, 200, 50, 50, 50, 50];

R=aphgrowthfun; % aphid growth from Asin 2001
%pred masses, abundances from Curtsdotter 2020
[I,W,N0,names]=genFields;
Nu=habitatfun;
solveresults(a0, phi, h0, b0, x0, Ropts, R, I, W, N0, Nu);
showresults(a0, phi, h0, b0, x0, Ropts, I, W, N0, Nu, names);

end

function solveresults(a0, phi, h0, b0, x0, Ropts, R, I, W, N0, Nu)

[a,ah,x]=atnParams(I,W,a0,phi,h0,Ropts,x0); %compute allometric params
tBd=[0,34]; %time to evaluate (is highest possible based on data)
fields={'JC','JO','KC','KO','MC','MO','OC','OO','SC','SO'}; %field names

for k=1:length(fields) %for all fields
    A=load(['Templog numDate comma ',fields{k},'.csv']); %load raw data
    Ts=A(:,2); %store temperatures
    t=A(:,1)-2; %store times
    
    %Find optimal community with these temps
    Opt=solveProb(a,ah,b0,R,x,Nu,tempdat2fun(t,Ts),tBd,N0,W);
    
    %compute future temps
    [Ti, Tv, Tiv]=tempchange(t,Ts,2.5,1.1);
    
    %find optimal community using future temps
    Bi=solveProb(a,ah,b0,R,x,Nu,tempdat2fun(t,Ti),tBd,N0,W); %inc avg
    Bv=solveProb(a,ah,b0,R,x,Nu,tempdat2fun(t,Tv),tBd,N0,W); %inc diff
    Biv=solveProb(a,ah,b0,R,x,Nu,tempdat2fun(t,Tiv),tBd,N0,W); %Inc both
    
    
    Bs=[Opt, Bi, Bv, Biv]; %store everything
    save(['optimization',fields{k}],'Bs'); %save
end

end

function yOpt=solveProb(a,ah,b0,R,x,Nu,Tf,tBd,N0,W)
%function obtains optimal predator community (yOpt) to minimize
%aphid abundance at temperatures Tf

%options for the matlab minimizer
opts = optimoptions('fmincon','Algorithm','interior-point');
th0=N0(2:end); %true predator community

%construct constraints
Acon=W(2:end); %for Ax=b, use A "sum up over predators times mass"
Bcon=Acon*th0; %for Ax=b, use b "original community biomass"

thLB=zeros(size(th0)); %lower bound is 'none of a predator'
thUB=ones(size(th0))*Bcon; %upper bound is 'total size of true community'

%initialize problem using fmincon and constraints above
problem=createOptimProblem('fmincon','objective',...
    @(th) controlCost(a,ah,b0,R,x,Nu,Tf,[N0(1);th],tBd),...
    'x0',th0,'Aeq', Acon, 'beq',Bcon, 'lb',thLB, 'ub',thUB,'options',opts);
ms=MultiStart; %create multistart problem
[yOpt,~]=run(ms,problem,20); %evaluate fmincon 20 times

yOpt=round(100*(yOpt./sum(yOpt))); %report 'percent of optimal predator community'

end

function [J,t,ys]=controlCost(a,ah,b0,R,x,Nu,Tf,N0,tBd)
%function returns the cost (average daily aphid abundance)
%for given predator community in N0 with temperatures Tf
[t,ys]=ode45(@atnODE,tBd(1):1:tBd(2),N0',[],a,ah,b0,R,x,Nu,Tf); %solve the ODE
J=sum(ys(:,1))/length(t); %compute average daily prey abundance
end

function score=overlapexplore(T,Nu)
%function computes average activity level of predators
%with activity Nu and at temperatures T
scores=zeros(9,1); %initialize
for j=1:length(T) %for all temps
    Tc=T(j); %take current temp
    Nuc=Nu(Tc); %evaluate activity
    score=scores+Nuc(1,2:end)'; %add activity with aphids
end
score=scores/length(T); %turn sum into an average
end

function dy=atnODE(t,y,a,ah,b0,R,x,Nu,Tf)
%function computes the change in populations y at time t
%for the allometric model (schneider) with temp-dependences
%in predator activity (Nu), death rates (x), and growth (R)
%and using temperature function Tf

T=Tf(t); %find current temperature
nu=Nu(T); %compute predator overlaps
r=zeros(size(y)); %make growth parameter of pop size
r(1)=R(T); %only aphids have growth
x=x*exp(-.69/(8.917e-5*(T+273.15))); %compute metabolic death 
a=a.*nu; %attack rates depend on overlap
ah=ah.*nu; %attack rates depend on overlap
F=ones(size(y))+(ah'+b0*a)*y; %compute functional response
effy=y./F; %scale predators by response
dy=y.*(r-a*effy-x'); %compute ODE
end

function [a,ah,x]=atnParams(I,W,a0,phi,h0,Ropt,x0)
%function returns matrices of attack rates (a), attack*handling (ah)
%and death rates (x) from allometric model (Schneider)

W2=W.^(1/2); 
W4=W.^(1/4); %powers appearing in model
A=(W4'*W4).*I; %encounter term for aij. I makes this zero if no consumption occurs
H=repmat(W2',1,length(I)).*I; %for handling time
R=((W.^(-1))'*W).*I; %success term for aij
X=W.^(3/4); %metabolic death term

Rmat=repmat(Ropt,length(Ropt),1); %stretch out Ropt values
Rm=R./Rmat; %compute Wj/Wi/Ropt
rat=(Rm.*exp(1-Rm)).^phi; %compute success curve in aij
a=a0.*A.*rat; %compute aij
ah=h0.*a0.*H.*rat; %compute aij*hij
x=x0*X; x(1)=0; %compute death rate xj, but ignore for prey
end

function tempfun=tempdat2fun(xipl, yipl)
%function turns vectors of temperature data (yipl) at times (xipl)
%into matlab function that draws line between values

 tempstring=['@(t) ',num2str(yipl(1)),'.*(t<=',num2str(xipl(1)),')']; %start a function
for jj=2:length(yipl) %for all data
    %make fun that takes this data between prev/curr times
    newfun=[num2str(yipl(jj)),'.*','((',num2str(xipl(jj-1)),'<t).*(t<=',num2str(xipl(jj)),'))'];
    tempstring=[tempstring,' + ',newfun]; %add it to the total function
end
%take on final value
tempstring=[tempstring,' + ',num2str(yipl(jj)),'.*(t>',num2str(xipl(jj)),')'];
tempfun=str2func(tempstring); %make it into a function for matlab
end

function [Ti, Tv, Tiv]=tempchange(t,T,m,v)
% function outputs modified temps from current T at time t
% increases by average m to create Ti
% increases variability by v to create Tv
% does both to create Tiv

tcurr=t(1); %start at first entry
tmax=ceil(t(end)); %continue until past final data
ts=tmax-tcurr; %how far to go in time
Dev=[]; %initialize
Avg=[];
count=1;
for j=1:ts %until you cover all times 
    Tsum=0; %start at 0
    day=floor(tcurr); %identify what day you are on (whole number)
    Ts=[]; %initialize
    Tnum=1;
    while tcurr<day+1 %until it's tomorrow
        Tsum=Tsum+T(count); %add current time to total
        Ts=[Ts, T(count)]; %store current time
        count=count+1; %move forward
        Tnum=Tnum+1; %add to count of temps
        tcurr=t(count); %take the new time
        if count==length(t) %if it is the end
            Tsum+T(count); %don't contribute to total
            Ts=[Ts, T(count)]; %store the temp
            Tnum=Tnum+1; %add counts
            tcurr=day+1; %add counts
        end
    end
    Tsum=Tsum/(Tnum-1); %turn total into average
    Ts=Ts-Tsum; %compute difference from average
    Ta=Tsum*ones(size(Ts)); %create vec of averages
    Dev=[Dev;Ts']; %store differences
    Avg=[Avg; Ta']; %store averages
end

Ti=(Avg+m)+Dev; %create temp with increase avg and original diff
Tv=Avg+v*Dev; %create temp with original avg and increase diff
Tiv=(Avg+m)+v*Dev; %create temp with increase avg and diff
figure, hold on,
plot(t,T, 'k-'), plot(t,Ti, 'r-'), plot(t,Tv, 'k--'), plot(t,Tiv, 'r--'), 
xlabel('time (days)'), ylabel('temperature (C)'),
legend('current','inc. avg','inc. var', 'inc. both'),
xlim([t(1), floor(t(1))+1]),
end

function Nu=habitatfun
%function outputs temperature-dependent Nu,
%   matrix of nu_i*nu_j for all species in system

f1=@(T) 1;%aphid
f2=@(T) max((4.2.*(1./(1+3.*exp(-.2.*(T-15)))).*(1-exp(.15.*(T-24)))),0);%pterostichus
f3=@(T) max((2.6.*(1./(1+3.*exp(-.2.*(T-13)))).*(1-exp(.15.*(T-27)))),0);%harpalus
f4=@(T) max((1.8.*(1./(1+3.*exp(-.2.*(T-9)))).*(1-exp(.15.*(T-28)))),0);%poecilus
f5=@(T) max((2.6.*(1./(1+3.*exp(-.2.*(T-8)))).*(1-exp(.15.*(T-22)))),0);%carabid
f6=@(T) max((2.*(1./(1+3.*exp(-.2.*(T-12)))).*(1-exp(.15.*(T-30)))),0);%bembidion
f7=@(T) max((2.1.*(1./(1+5.*exp(-.2.*(T-20)))).*(1-exp(.15.*(T-38)))),0);%other spider
f8=@(T) max((2.1.*(1./(1+5.*exp(-.2.*(T-20)))).*(1-exp(.15.*(T-38)))),0);%lyco
f9=@(T) max((1.5.*(1./(1+5.*exp(-.3.*(T-15)))).*(1-exp(.15.*(T-25)))),0); %tetra 
f10=@(T) max((.5.*(1./(1+5.*exp(-.3.*(T-10)))).*(1-exp((.15.*(T-30))))),0); %Liny

% %foraging similarity between prey and predators
Nu=@(T) [ f1(T).*f1(T) f1(T).*f2(T) f1(T).*f3(T) f1(T).*f4(T) f1(T).*f5(T) f1(T).*f6(T) f1(T).*f7(T) f1(T).*f8(T) f1(T).*f9(T) f1(T).*f10(T);...
    f2(T).*f1(T) f2(T).*f2(T) f2(T).*f3(T) f2(T).*f4(T) f2(T).*f5(T) f2(T).*f6(T) f2(T).*f7(T) f2(T).*f8(T) f2(T).*f9(T) f2(T).*f10(T);...
    f3(T).*f1(T) f3(T).*f2(T) f3(T).*f3(T) f3(T).*f4(T) f3(T).*f5(T) f3(T).*f6(T) f3(T).*f7(T) f3(T).*f8(T) f3(T).*f9(T) f3(T).*f10(T);...
    f4(T).*f1(T) f4(T).*f2(T) f4(T).*f3(T) f4(T).*f4(T) f4(T).*f5(T) f4(T).*f6(T) f4(T).*f7(T) f4(T).*f8(T) f4(T).*f9(T) f4(T).*f10(T);...
    f5(T).*f1(T) f5(T).*f2(T) f5(T).*f3(T) f5(T).*f4(T) f5(T).*f5(T) f5(T).*f6(T) f5(T).*f7(T) f5(T).*f8(T) f5(T).*f9(T) f5(T).*f10(T);...
    f6(T).*f1(T) f6(T).*f2(T) f6(T).*f3(T) f6(T).*f4(T) f6(T).*f5(T) f6(T).*f6(T) f6(T).*f7(T) f6(T).*f8(T) f6(T).*f9(T) f6(T).*f10(T);...
    f7(T).*f1(T) f7(T).*f2(T) f7(T).*f3(T) f7(T).*f4(T) f7(T).*f5(T) f7(T).*f6(T) f7(T).*f7(T) f7(T).*f8(T) f7(T).*f9(T) f7(T).*f10(T);...
    f8(T).*f1(T) f8(T).*f2(T) f8(T).*f3(T) f8(T).*f4(T) f8(T).*f5(T) f8(T).*f6(T) f8(T).*f7(T) f8(T).*f8(T) f8(T).*f9(T) f8(T).*f10(T);...
    f9(T).*f1(T) f9(T).*f2(T) f9(T).*f3(T) f9(T).*f4(T) f9(T).*f5(T) f9(T).*f6(T) f9(T).*f7(T) f9(T).*f8(T) f9(T).*f9(T) f9(T).*f10(T);...
    f10(T).*f1(T) f10(T).*f2(T) f10(T).*f3(T) f10(T).*f4(T) f10(T).*f5(T) f10(T).*f6(T) f10(T).*f7(T) f10(T).*f8(T) f10(T).*f9(T) f10(T).*f10(T)];
end

function [IM,W,y0,spec]=genFields
%function outputs quantities needed for ATN model solutions: 
% interaction matrix (IM), vector of masses (W), 
% initial abundancecs (y0), and names for all species (spec)
%using data from Curtstdotter 2020 (must be available in directory)

fields={'JC','JO','KC','KO','MC','MO','OC','OO','SC','SO'}; %field refs
specs={'Aphid','Bembidion','Harpalus','Poecilus','Pterostichus','Other_Carabid',...
    'Linyphiidae','Lycosidae','Tetragnathidae','Other_Spider'}; %pred refs

Wtot=zeros(1,10); ytot=zeros(10,1); %initialize vectors

for k=1:length(fields) %loop over fields
    Wfull=load(['./BM FoodWeb Fulltime ',fields{k},'.csv']); %load the masses
    W=[Wfull(1),Wfull(6:14)']; %save masses for prey, beetles, and spiders
    y0=zeros(length(W),1); %initialize vector
    for j=2:length(specs) %loop over predators
        dats=load(['./Density numDate ',specs{j},' Fulltime ',fields{k},'.csv']); %load abundances
        pop=dats(:,2); %discard times
        averagePop=mean(pop); %compute average abundance
        y0(j)=averagePop; %save this predator's average abundance
    end
    aphs=load(['./Density numDate Aphid Shifted Truncated ',fields{k},'.csv']); %load abundances
    y0(1)=aphs(1,2); %save initial abundance
    Wtot=Wtot+W; %add up the masses we've found
    ytot=ytot+y0; %add up the initial abundances we've found
end
Wtot=Wtot./length(fields); %compute average mass over all fields

[WB,Bind]=sort(-1*Wtot(2:6)); %arrange beetles by size (easier to view)
[WS,Sind]=sort(-1*Wtot(7:end)); %arrange spiders by size (easier to view)
W=[Wtot(1),-1*WB,-1*WS]; %re-pack all species

ytot=ytot./length(fields); %compute average init. abund. over all fields
yB=ytot(2:6); yS=ytot(7:end); %separate beetles/spiders
y0=[ytot(1);yB(Bind);yS(Sind)]; %pack using ordering from above (by mass)

specs={'Aphid','Bembidion','Harpalus','Poecilus','Pterostichus','Other Carabid',...
    'Linyphiidae','Lycosidae','Tetragnathidae','Other Spider'}; %readable species names

spec={specs{Bind(1)+1},specs{Bind(2)+1},specs{Bind(3)+1},specs{Bind(4)+1},specs{Bind(5)+1},...
    specs{Sind(1)+6},specs{Sind(2)+6},specs{Sind(3)+6},specs{Sind(4)+6}}; 
%pack using ordering from above (by mass)

IMfull=load('./AlvaDat/FoodWeb.csv'); %load interaction matrix
IMhalf=[IMfull(1,:);IMfull(6:14,:)]; %save the part for prey, beetles, spiders
IM=[IMhalf(:,1),IMhalf(:,6:14)];
IM(9,2)=1; %real foodweb introduced inconsistencies
%the only non-consumption is Bembidion won't eat tetragnathidae, 
%allow this link to make comparison possible between groups
end

function R=aphgrowthfun
% function fits model to data from Asin, Pons 2001
% and outputs temperature-dependent aphid growth

Ts=[18, 22, 25, 27.5, 30]; %Asin, Pons, Env Ent 2001
rs=[0.26, 0.37, 0.45, 0.52, 0.05]; %Asin, Pons, Env Ent 2001
dat=[Ts; rs]; %pack data
th0=[.55; 29; .5; .5]; %initial guess
thL=[.3; 25; .0001; .0001]; %min param
thU=[.7; 33; 10; 10]; %max param
th=fmincon(@(th) aphgrowthfit(dat, th(1), th(2), th(3), th(4)),...
    th0, [], [], [], [], thL, thU); %fit to data
R=@(T) th(1)*(exp(-.5*((T-th(2))/th(3)).^2).*(T<=th(2))+...
    exp(-.5*((T-th(2))/th(4)).^2).*(T>th(2))); %create function
end

function J=aphgrowthfit(dat, rm, tu, c1, c2)
% function outputs difference between data (dat) and model
% from Asin/Pons 2001, given model params rm, tu, c1, c2

Rtf=@(T) rm*(exp(-.5*((T-tu)/c1).^2).*(T<=tu)+...
    exp(-.5*((T-tu)/c2).^2).*(T>tu)); %create function
Rdat=Rtf(dat(1,:)); %evaluate at times
J=sum((Rdat-dat(2,:)).^2); %compute distance from data

end

function showresults(a0, phi, h0, b0, x0, Ropts, I, W, N0, Nu, names)
%creates plots used in manuscript, requires results saved from solveresults
%note: inefficient/repetitive so that we can modify individual figures
%uncomment relevant figure to edit/save that one

[a,ah,~]=atnParams(I,W,a0,phi,h0,Ropts,x0); %create allometric parameters
fields={'JC','JO','KC','KO','MC','MO','OC','OO','SC','SO'}; %field names
pcol=[72, 0, 145; 0, 109, 218; 182, 109, 255; 109, 181, 255; 182, 219, 255;...
    146, 0, 0; 218, 108, 0; 226, 1, 52; 255, 255, 109]/255; %colors for plots

ts=0:.1:34; %times to show results
Tspan=5:.1:50; %temperatures to show results

%% Fig 1 and supplemental temp
% 
% figure, hold on,
% Temps=[];
% for k=1:length(fields) %for all fields
%     A=load(['Templog numDate comma ',fields{k},'.csv']); %load raw data
%     Ts=A(:,2); %store temps
%     t=A(:,1)-2; %store times
%     Tf=tempdat2fun(t,Ts); %create temperature function
%     Ts=Tf(ts); %evaluate function at times to show
%     Temps=[Temps, Ts]; %store all temperatures for later
%     plot(ts,Ts) %plot temperatures
% end
% xlim([0,34]), xlabel('time (days)'), ylabel('temperature (C)')
% 
% As=zeros(length(Tspan),9); %initialize
% for j=1:length(Tspan) %for all temperatures to show
%    Tc=Tspan(j); %look at current temperature
%    nu=Nu(Tc); %evaluate current activity
%    an=a.*nu; %evaluate current attack rate based on activity
%    As(j,:)=an(1,2:end); %store attacks on aphid
% end
% figure, hold on,
% subplot(3,1,1), hold on,
% for p=1:5 %for beetle predators
%     %plot the attack rate on aphids
%     plot(Tspan,As(:,p),'-','Color',pcol(p,:),'LineWidth',3),
% end
% ylabel('attack rate on aphid'),
% legend(names{1:5})
% ylim([0 25]), xlim([5 50])
% subplot(3,1,2), hold on,
% for p=6:9 %for spider predators
%     %plot the attack rate on spiders
%     plot(Tspan,As(:,p),'-','Color',pcol(p,:),'LineWidth',3),
% end
% legend(names{6:9})
% ylabel('attack rate on aphid'),
% ylim([0 25]), xlim([5 50])
% subplot(3,1,3), hold on,
% histogram(Temps)
% xlabel('temperature (C)'),
% ylabel('frequency')
% xlim([5,50])

%% supplemental pred activity
% 
% As=zeros(length(Tspan),9); %initialize
% for j=1:length(Tspan) %for all temps to show
%    Tc=Tspan(j); %take current temp
%    nu=Nu(Tc); %compute activity
%    As(j,:)=nu(1,2:end); %store activity with aphids (it's just predator activity)
% end
% figure, hold on,
% subplot(2,1,1), hold on,
% for p=1:5 %For beetle predators
%     %plot activiity
%     plot(Tspan,As(:,p),'-','Color',pcol(p,:),'LineWidth',3),
% end
% ylabel('activity level'),
% legend(names{1:5})
% ylim([0 1]), xlim([5 50])
% subplot(2,1,2), hold on,
% for p=6:9 %for spider predators
%     %plot activity
%     plot(Tspan,As(:,p),'-','Color',pcol(p,:),'LineWidth',3),
% end
% legend(names{6:9})
% ylabel('activity level'),
% ylim([0 1]), xlim([5 50])

%% Fig 2
% 
% tss=3:.01:5; %take a finer timespan
% A=load(['Templog numDate comma ',fields{10},'.csv']); %load data
% Ts=A(:,2); %store temps
% t=A(:,1)-2; %store times
% Tf=tempdat2fun(t,Ts); %create temperature function
% A1s=zeros(length(tss),5); %evaluate at fine times
% N0(1)=0; %ignore aphid abundance
% for j=1:length(tss) %for the fine timespan
%     Tc=Tf(tss(j)); %look at current temp
%     nu=Nu(Tc); %compute activity
%     an=a.*nu; %Compute attacks
%     ahn=ah.*nu; %compute handling times
%     F=ones(size(N0))+(ahn'+b0*an)*N0; %compute functional response
%     A1s(j,5)=Tc; %store the temp
%     A1s(j,1)=an(1,4); %store attack rates on aphids without FR by Poecilus
%     A1s(j,2)=an(1,4)/F(4); %store attack rates on aphids with FR by poecilus
%     A1s(j,3)=an(1,7); %store attack rates on aphids without FR by Other Spider
%     A1s(j,4)=an(1,7)/F(7); %store attack rates on aphids with FR by Other Spider
% end
% figure, hold on,
% subplot(3,1,1), hold on,
% plot(tss, A1s(:,1),'-','Color',pcol(3,:),'LineWidth',3),
% plot(tss, A1s(:,2),'--','Color',pcol(3,:),'LineWidth',3),
% ylabel('Poecilus attack rate'),
% xlim([3,5]), xticks([3, 3.5, 4, 4.5, 5]), xticklabels({'0','12','24','36','48'})
% subplot(3,1,2), hold on,
% plot(tss,A1s(:,3),'-','Color',pcol(6,:),'LineWidth',3),
% plot(tss,A1s(:,4),'--','Color',pcol(6,:),'LineWidth',3),
% ylabel('Other Spider attack rate'),
% xlim([3,5]), xticks([3, 3.5, 4, 4.5, 5]), xticklabels({'0','12','24','36','48'})
% subplot(3,1,3),
% plot(tss,A1s(:,5),'k-','LineWidth',3),
% ylabel('temperature (C)'),
% xlabel('time (hours)'),
% xlim([3,5]), xticks([3, 3.5, 4, 4.5, 5]), xticklabels({'0','12','24','36','48'})

%% Figure 3

% A=zeros(9,10); %initialize
% B=A;
% Biv=A;
% S=A;
% 
% old=N0(2:end); %original predator abundances
% 
% for k=1:length(fields) %for all fields
%     load(['optimization',fields{k}],'Bs'); %load optimization solution
%     est=Bs(:,1); %first entry is the optimal community percentages
%     abun=(W(2:end)*old)/((W(2:end)*est)/100); %conversion factor to abundance
%           %original biomass / 100*percentages biomass
%     new1=abun*est/100; %convert percentages to abundance
%     B(:,k)=new1; %store the abundance of optimal community
%     
%     A=load(['Templog numDate comma ',fields{k},'.csv']); %load temperature data
%     Ts=A(:,2); %store temperatures
%     S(:,k)=overlapexplore(Ts,Nu); %report predator activity at temperatures
% end
% S=mean(S'): %compute average predator activity
% figure, hold on,
% subplot(2,1,1),
% yyaxis left
% bar(.875:1:8.875, S, .2,'LineStyle','-','FaceColor',[0 0.4470 0.7410])
% ylabel('average activity'),
% yyaxis right
% bar(1.125:1:9.125, W(2:end),.2,'LineStyle','-','FaceColor',[109/255, 182/255, 1]);
% ylabel('body size (mg)'),
% ax=gca;
% ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = 'k';
% xlim([0.5,9.5])
% 
% xticklabels({})
% specs=names;
% subplot(2,1,2),
% boxchart(B')
% ylim([0 .8])
% ylabel('optimal density')
% xticklabels(specs)
% xtickangle(45)

%% Figure 4

% A=zeros(9,10); %initialize
% B=A;
% Biv=A;
% S=A;
% 
% for k=1:length(fields) %for all fields
%     load(['optimization',fields{k}],'Bs'); %load optimization
%     B(:,k)=Bs(:,1); %store optimal community
% end
% 
% old=N0(2:end);%original predator community
% est=mean(B')'; %compute average optimal community
% abun=(W(2:end)*old)/((W(2:end)*est)/100); %conversion factor to abundance
% new=abun*est/100; %turn percents into abundance
% N0=[0; new]; %store optimal predators with no aphid
% 
% AnoT=zeros(length(Tspan),1); %initialize
% AwT=AnoT;
% AFnoT=AnoT;
% AFwT=AnoT;
% 
% noT1=a*N0; %attack rates without temperatures or FR
% noT1=noT1(1); %attack rates on aphid without temperatures or FR
% F=ones(size(N0))+(ah+b0*a)*N0; %functional response without temperatures
% effy=N0./F; %scale predators by functional response
% noT2=(a*effy); %attack rates without temperatures but with FR
% noT2=noT2(1); %attack rates on aphids without temperature but with FR
% 
% for j=1:length(Tspan) %for temps to show
%        Tc=Tspan(j); %look at current temp
%        nu=Nu(Tc); %compute activity
%        an=a.*nu; %compute attacks with activity
%        ahn=ah.*nu; %compute handling time with activity
%        F=ones(size(N0))+(ahn'+b0*an)*N0; %compute functional response
%        effy=N0./F; %scale predators by response
%        dr1=(an*N0); %attack rates with temperature and without FR
%        dr2=(an*effy); %attack rates with temperatures and FR
%        AwT(j)=dr1(1); %store attacks on aphids with temp / without FR
%        AFwT(j)=dr2(1); %store attacks on aphids with temp/ with FR
%        AnoT(j)=noT1; %store attacks on aphids without temp / without FR (never changes)
%        AFnoT(j)=noT2; %store attacks on aphids without temp / without FR (never changes)
% end
% 
% figure, hold on, plot(Tspan,AnoT,'k-','LineWidth',3),plot(Tspan,AFnoT,'k--','LineWidth',3),
% plot(Tspan,AwT,'b-','LineWidth',3),plot(Tspan,AFwT,'b--','LineWidth',3),
% xlim([5 50]), 
% ylim([0 7.5])
% xlabel('temperature (C)'), ylabel('attacks on aphids'),
% legend({'no temp., no F.R.','no temp., with F.R.','with temp., no F.R.','with temp., with F.R.'});

%% Figure 5

% for k=1:length(fields) %for all fields
%     load(['optimization',fields{k}],'Bs'); %load optimization result
%     est=Bs(:,1); %grab optimal community
%     abun=(W(2:end)*old)/((W(2:end)*est)/100); %conversion factor to abundance
%     new1=abun*est/100; %Convert optimal community to abundance
%     est=Bs(:,4); %grab optimal community under future temps (inc avg and inc diff)
%     abun=(W(2:end)*old)/((W(2:end)*est)/100); %conversion factor to abundance
%     new2=abun*est/100;  %Convert optimal community to abundance
%     
%     Biv(:,k)=-(new1-new2); %report change in optimal communities
% end
% specs=names;
% figure, hold on,
% boxchart(Biv'),
% yline(0,'k--')
% ylabel('change in optimal density')
% xticklabels(specs)
% xtickangle(45)    

end
