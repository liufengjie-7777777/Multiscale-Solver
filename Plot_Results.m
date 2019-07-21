%Plot Results

%%Analyze biaxial Simulation
b = SimArteryVessel;

%Update here material parameters that changed

%Update vectors
b.riVec = a.V.ri;
b.ufsVec = a.V.ufsN;
b.nAMpVec = a.V.nAMp;
b.nAMVec = a.V.nAM;
b.timeVec = a.V.time;

%Plot Cauchy stress as a function of normalized radii (passive and
%steady-state active
N = [1,length(b.riVec)];
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
    end
end

yLabel = {'\sigma_r','\sigma_\theta','\sigma_z'};
lShape = {'--','-'};
x = (x-b.ri)./(b.ro-b.ri);
figure();
for k=1:3
    subplot(1,3,k); hold on;
    for i=1:length(N)
        plot(x,S(:,i,k)*1e3,lShape{i});
    end
    hold off;
    ylabel([yLabel{k} ' (kPa)']);
    xlabel('Normalized radial position');
    legend('Relaxed','SS Active');
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

yLabel = {'dMA_r/dMA0','L_{M,r} (\mum)','L_{fo,r}'};
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
    legend('Relaxed','SS Active');
end

%%
%Plot SMC parameters changing through time for a specific radii
N = 1; %Number of radii coordinates to calculate

dMArMA0 = zeros(length(b.timeVec),N);
LMr = zeros(length(b.timeVec),N);
Lfor = zeros(length(b.timeVec),N);
eS2xr = zeros(length(b.timeVec),N);
eS2yr = zeros(length(b.timeVec),N);

for i=1:length(b.timeVec)
    b.ri = b.riVec(i);
    b.ufs = b.ufsVec(i,:);
    b.nAMp = b.nAMpVec(i);
    b.nAM = b.nAMVec(i);
    
    r = (b.ri+b.ro)/2;
    
    for j=1:length(r)
        b.Lfoi(r(j));
        b.eS2;
        
        dMArMA0(i,j) = b.dMArMA0;
        LMr(i,j) = b.LMr;
        Lfor(i,j) = b.Lfor;
        
        eS2xr(i,j) = b.eS2xr;
        eS2yr(i,j) = b.eS2yr;
    end
end

figure();
plot(b.timeVec,dMArMA0); title('dMAr/dMA0 vs time');
ylabel('dMAr/dMA0'); xlabel('time (min)');

figure();
plot(b.timeVec,LMr*1e3);  title('L_{M,r} vs time');
ylabel('LMr (\mum)'); xlabel('time (min)');


figure();
plot(b.timeVec,Lfor); title('L_{fo,r} vs time');
ylabel('L_{fo,r}'); xlabel('time (min)');    

figure();
plot(b.timeVec,eS2xr,b.timeVec,eS2yr); title('\epsilon_{S2,r} vs time');
legend('\epsilon_{S2xr}','\epsilon_{S2yr}');
ylabel('\epsilon_{S2,r}'); xlabel('time (min)');
%ylim([0.06 0.15]);
