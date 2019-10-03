%Plot Results

%%Analyze a Simulation
Sim = 'Biaxial'; %'Uniaxial'; %
varName = 'alphaPS';
n = 1;

%load(['Simulation Results\' Sim 'Simulation-' varName '(' num2str(n) ').mat']);

b = SimArteryVessel;

%Update material parameters that changed
b.UpdateParameters(a);

%Update vectors
b.riVec = a.V.ri;
b.ufsVec = a.V.ufsN;
b.nAMpVec = a.V.nAMp;
b.nAMVec = a.V.nAM;
b.timeVec = a.V.time;

%%
%Plot Do and FT as a function of time for the last loaded simulation
Do = zeros(length(b.timeVec),1);
FT = zeros(length(b.timeVec),1);
for i=1:length(b.timeVec)
    b.ri = b.riVec(i);
    b.ufs = b.ufsVec(i,:);
    b.nAMp = b.nAMpVec(i);
    b.nAM = b.nAMVec(i);
    
    Do(i) = b.ro*2e3; %um
    FT(i) = b.FTCalc*1e3; %mN
end

figure();
%hold on;
subplot(1,2,1);
plot(b.timeVec,Do);
ylabel('Do (\mum)'); xlabel('time (min)');
ylim(limCalc(Do,[0.5 1.2],0));
subplot(1,2,2);
plot(b.timeVec,FT);
ylabel('F_T (mN)'); xlabel('time (min)');
ylim(limCalc(FT,[0.7 1.3],0));
hold off;
%%

%Plot Do and FT as a function of time for multiple simulations
N = 4; %number of files to open
Sim = 'Biaxial';
varName = 'LS20';

%legendOptions
%['\alpha_{PS} = ' num2str(b.alphaPS*180/pi) '^o']

b = SimArteryVessel;
Do = zeros(length(b.timeVec),N);
FT = zeros(length(b.timeVec),N);
strLegend = strings(1,N);
for n=1:N
    load(['Simulation Results\' Sim 'Simulation-' varName '(' num2str(n) ').mat']);
    
    %Update material parameters that changed
    b.UpdateParameters(a);
    strLegend{n} = [varName ' = ' num2str(b.(varName)) ' (?)'];
    
    %Update vectors
    b.riVec = a.V.ri;
    b.ufsVec = a.V.ufsN;
    b.nAMpVec = a.V.nAMp;
    b.nAMVec = a.V.nAM;
    b.timeVec = a.V.time;
    
    for i=1:length(b.timeVec)
        b.ri = b.riVec(i);
        b.ufs = b.ufsVec(i,:);
        b.nAMp = b.nAMpVec(i);
        b.nAM = b.nAMVec(i);
        
        Do(i,n) = b.ro*2e3; %um
        FT(i,n) = b.FTCalc*1e3; %mN
    end
end
figure();
subplot(1,2,1);
plot(b.timeVec,Do);
ylabel('Do (\mum)'); xlabel('time (min)');
ylim(limCalc(Do,[0.5 1.2],0));
legend(strLegend);
subplot(1,2,2);
plot(b.timeVec,FT);
ylabel('F_T (mN)'); xlabel('time (min)');
ylim(limCalc(FT,[0.7 1.3],0));
legend(strLegend);

%%
%Plot Cauchy stress as a function of normalized radii and time

%Choose N to show results as a function of time or Passive and SS values
N = [1,length(b.riVec)]; %round(linspace(1,length(b.riVec),10)); %
p = zeros(20,length(N));
S = zeros(20,length(N),3);
for i=1:length(N)
    b.ri = b.riVec(N(i));
    b.ufs = b.ufsVec(N(i),:);
    b.nAMp = b.nAMpVec(N(i));
    b.nAM = b.nAMVec(N(i));
    
    x = linspace(b.ri,b.ro,20);
    for j=1:length(x)
        [p(j,i),S(j,i,:)] = b.CauchyStress(x(j));
        if isnan(p(j,i))
            disp(['p is NaN: ' num2str([j,i])]);
        end
    end
end
yLabel = {'\sigma_r','\sigma_\theta','\sigma_z'};
x = (x-b.ri)./(b.ro-b.ri);

if length(N)>2
    for k=1:3
        figure();
        mesh(b.timeVec(N),x,S(:,:,k)*1e3);
        
        zlabel([yLabel{k} ' (kPa)']);
        xlabel('time (min)');
        ylabel('Normalized position');

    end
else
    lShape = {'--','-'};
    figure();
    for k=1:3
    subplot(1,3,k); hold on;
    for i=1:length(N)
        plot(x,S(:,i,k)*1e3,lShape{i});
        lim = [floor(0.95*min(S(:,:,k)*1e3,[],'all')) ceil(1.01*max(S(:,:,k)*1e3,[],'all'))];
        if lim(2) > lim(1)
            ylim(lim);
        end
    end
    hold off;
    ylabel([yLabel{k} ' (kPa)']);
    xlabel('Normalized radial position');
    legend('Relaxed','Steady-state Active');
    end
end

%%
%Plot SMC parameters as a function of normalized radii (passive and
%steady-state active

N = [1,length(b.riVec)];

dMArMA0 = zeros(20,length(N));
LMr = zeros(20,length(N));
Lfor = zeros(20,length(N));

for i=1:length(N)
    b.ri = b.riVec(N(i));
    b.ufs = b.ufsVec(N(i),:);
    b.nAMp = b.nAMpVec(N(i));
    b.nAM = b.nAMVec(N(i));

    x = linspace(b.ri,b.ro,20);
    b.Lfoi(x);
    dMArMA0(:,i) = b.dMArMA0;
    LMr(:,i) = b.LMr;
    Lfor(:,i) = b.Lfor;
end
curves(:,:,1) = dMArMA0;
curves(:,:,2) = LMr*1e3;
curves(:,:,3) = Lfor;

yLabel = {'dMAr/dMA0','L_{M,r} (\mum)','L_{fo,r}'};
lShape = {'--','-'};
x = (x-b.ri)./(b.ro-b.ri);
figure();
for k=1:3
    subplot(1,3,k); hold on;
    for i=1:length(N)
        plot(x,curves(:,i,k),lShape{i});
    end
    hold off;
    ylabel([yLabel{k}]);
    xlabel('Normalized radial position');
    legend('Relaxed','Steady-state Active');
    lim = [round(95*min(curves(:,:,k),[],'all'))/100 round(105*max(curves(:,:,k),[],'all'))/100];
    if lim(2) > lim(1)
        ylim(lim);
    end
end

%%
%Plot SMC parameters changing through time for a specific radii
N = 1; %Number of radii coordinates to calculate

strdMAiMA0 = {'dMAr/dMA0','dMAz/dMA0'};
strLMi = {'L_{M,r}','L_{M,z}'};
strLfoi = {'L_{fo,r}','L_{fo,z}'};
streS2i = {'e_{S2,r}','e_{S2,z}'};

for k=1:2 %1-r and 2-z direction
    dMAiMA0 = zeros(length(b.timeVec),N);
    LMi = zeros(length(b.timeVec),N);
    Lfoi = zeros(length(b.timeVec),N);
    eS2i = zeros(length(b.timeVec),N,2);
    
    for i=1:length(b.timeVec)
        b.ri = b.riVec(i);
        b.ufs = b.ufsVec(i,:);
        b.nAMp = b.nAMpVec(i);
        b.nAM = b.nAMVec(i);
        
        r = (b.ri+b.ro)/2;
        
        for j=1:length(r)
            b.Lfoi(r(j));
            b.eS2;
            
            if k==1
                dMAiMA0(i,j) = b.dMArMA0;
                LMi(i,j) = b.LMr;
                Lfoi(i,j) = b.Lfor;
                eS2i(i,j,:) = [b.eS2xr, b.eS2yr];
            elseif k==2
                dMAiMA0(i,j) = b.dMAzMA0;
                LMi(i,j) = b.LMz;
                Lfoi(i,j) = b.Lfoz;
                eS2i(i,j,:) = [b.eS2xz, b.eS2yz];
            end
        end
    end
    
    figure();
    subplot(2,2,1);
    plot(b.timeVec,dMAiMA0); title([strdMAiMA0{k} ' vs time']);
    ylabel(strdMAiMA0{k}); xlabel('time (min)');
    ylim(limCalc(dMAiMA0,[0.9 1.05],2));
    subplot(2,2,2);
    plot(b.timeVec,LMi*1e3);  title([strLMi{k} ' vs time']);
    ylabel([strLMi{k} ' (\mum)']); xlabel('time (min)');
    ylim(limCalc(LMi,[0.9 1.02],2));
    subplot(2,2,3);
    plot(b.timeVec,Lfoi); title([strLfoi{k} ' vs time']);
    ylabel(strLfoi{k}); xlabel('time (min)');
    ylim(limCalc(Lfoi,[0.9 1.02],2));
    subplot(2,2,4);
    strDir = ['r','z'];
    plot(b.timeVec,[eS2i(:,1,1), eS2i(:,1,2)]); title([streS2i{k} ' vs time']);
    legend(['\epsilon_{S2x,' strDir(k) '}'],['\epsilon_{S2y,' strDir(k) '}']);
    ylabel(streS2i{k}); xlabel('time (min)');
    ylim(limCalc(eS2i,[0.9 1.1],2));
end


%%
%Plot SS Do and FT in parameter scaling

Sim = 'Biaxial';
%varName = {'alphaPS','AS2','beta','deltam','dMA0','LLA0','lMD','LMmax','LS20'};
varName = {'LLA0','lMD','LMmax','LS20','dMA0','deltam','alphaPS','AS2'};
N = length(varName); %number of files to open

%strLegend = {'\alpha_{PS}','A_{S2}','\beta','\delta_m','dMA0','L_{LA0}','l_{MD}','L_{Mmax}','L_{S20}'};
strLegend = {'L_{LA0}','l_{MD}','L_{Mmax}','L_{S20}','dMA0','\delta_m','\alpha_{PS}','A_{S2}'};
marker = {'+','o','*','.','s','d','^','v','<','>','p','h'};
lineSpec = {'-','--','-.',':','-','--','-.',':','-'};

%legendOptions
%['\alpha_{PS} = ' num2str(b.alphaPS*180/pi) '^o']

b = SimArteryVessel;
Do = zeros(25,N);
FT = Do;
%strLegend = strings(1,N);
for n=1:N
    for k=1:size(Do,1)
        load(['Simulation Results2\' varName{n} ' Simulations\' Sim 'Simulation-' varName{n} '(' num2str(k) ').mat']);

        %Update material parameters that changed
        b.UpdateParameters(a);
        %strLegend{n} = [varName ' = ' num2str(b.(varName)) ' (?)'];

        %Update vectors
        b.riVec = a.V.ri;
        b.ufsVec = a.V.ufsN;
        b.nAMpVec = a.V.nAMp;
        b.nAMVec = a.V.nAM;
        b.timeVec = a.V.time;
        
        i = length(b.timeVec);
        
        b.ri = b.riVec(i);
        b.ufs = b.ufsVec(i,:);
        b.nAMp = b.nAMpVec(i);
        b.nAM = b.nAMVec(i);
        
        Do(k,n) = b.ro*2e3; %um
        FT(k,n) = b.FTCalc*1e3; %mN
    end
end

DoCom = 834.6234;
FTCom  = 4.5765;

DoErr = (Do-DoCom)/DoCom*100;
FTErr = (FT-FTCom)/FTCom*100;
scaling = linspace(0.4,1.6,size(DoErr,1));

figure();
subplot(1,2,1);
for n=1:4
    plot(scaling,DoErr(:,n),lineSpec{n});
    hold on;
end
hold off;
legend(strLegend{1:n});
ylabel('Do Difference (%)'); xlabel('Parameter Scaling (-)');
xlim([0.4 1.6]); grid on;

%figure();
subplot(1,2,2);

lineSpec = {'-','--','-.',':','-','--','-.',':','-'};

for n=5:length(strLegend)
    plot(scaling,DoErr(:,n),lineSpec{n});
    hold on;
end
hold off;
legend(strLegend{5:end});
ylabel('Do Difference (%)'); xlabel('Parameter Scaling (-)');
xlim([0.4 1.6]); grid on;

figure();
subplot(1,2,1);
lineSpec = {'-','--','-.',':','-','--','-.',':','-'};

for n=1:4
    plot(scaling,FTErr(:,n),lineSpec{n});
    hold on;
end
hold off;
legend(strLegend{1:n});
ylabel('F_T Difference (%)'); xlabel('Parameter Scaling (-)');
xlim([0.4 1.6]); grid on;

%figure();
subplot(1,2,2);
lineSpec = {'-','--','-.',':','-','--','-.',':','-'};

for n=5:length(strLegend)
    plot(scaling,FTErr(:,n),lineSpec{n});
    hold on;
end

hold off;
legend(strLegend{5:end});
ylabel('F_T Difference (%)'); xlabel('Parameter Scaling (-)');
xlim([0.4 1.6]); grid on;

%%
%Plot SS Do and FT for different Pin and lz

Sim = 'Biaxial';

n = [16,8];

N = n(1)*n(2); %number of files to open

%strLegend = {};
%marker = {'+','o','*','.','s','d','^','v','<','>','p','h'};
%lineSpec = {'-','--','-.',':','-','--','-.',':','-'};

b = SimArteryVessel;
Do = zeros(n); FT = Do;

lz = zeros(n(1),1); %linspace(0.5,2,16);
Pin = zeros(n(2),1); %linspace(40,110,8);

for k1=1:n(1)
    for k2=1:n(2)
        load(['Simulation lzandPin\' Sim ' Simulation' '(' num2str(n(2)*(k1-1)+k2) ').mat']);
        
        %Update material parameters that changed
        b.UpdateParameters(a);
        %strLegend{n} = [varName ' = ' num2str(b.(varName)) ' (?)'];

        %Update vectors
        b.riVec = a.V.ri;
        b.ufsVec = a.V.ufsN;
        b.nAMpVec = a.V.nAMp;
        b.nAMVec = a.V.nAM;
        b.timeVec = a.V.time;
        
        i = length(b.timeVec);
        
        b.ri = b.riVec(i);
        b.ufs = b.ufsVec(i,:);
        b.nAMp = b.nAMpVec(i);
        b.nAM = b.nAMVec(i);
        
        Do(k1,k2) = b.ro*2e3; %um
        FT(k1,k2) = b.FTCalc*1e3; %mN
        
        lz(k1) = b.lz;
        if k1==1
            Pin(k2) = b.Pin/(133.322387415*1e-6); %Convert to mmHg
        end
    end
end

figure();
for k2=1:n(2)
    plot(lz,Do(:,k2));
    hold on;
    strLegend{k2} = [num2str(Pin(k2)) ' mmHg'];
end
hold off;
legend(strLegend);
xlabel('\lambda_z'); ylabel('Do (\mum)');



