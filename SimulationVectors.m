classdef SimulationVectors < handle
    %Simulationtors contains the variables calculated during simulation
    properties
        %Both Simulations
        ufs
        stretch
        sECM
        sSMC
        sMMy
        PMMCU
        
        nAMp
        nAM
        time
        
        %Uniaxial Simulation
        Pisom
        lr
        
        dMAiMA0
        LMi
        Lfoi
        eS2
        I4f
        x2c
        y2c
        
        %Biaxial Simulation
        ufsN
        ri
        Do
        FT
    end
    
    methods
        function obj = InitialVectors(obj,SamplePoints,Sim)
            %Both
            obj.ufs = zeros(SamplePoints,1);
            obj.stretch = zeros(SamplePoints,4);
            
            obj.sECM = zeros(SamplePoints,3);
            obj.sSMC = zeros(SamplePoints,3);
            obj.sMMy = zeros(SamplePoints,3);
            
            obj.PMMCU = zeros(SamplePoints,2);
            
            obj.dMAiMA0 = zeros(SamplePoints,2);
            obj.LMi = zeros(SamplePoints,2);
            obj.Lfoi = zeros(SamplePoints,2);
            obj.eS2 = zeros(SamplePoints,4);
            obj.x2c = zeros(SamplePoints,2);
            obj.y2c = zeros(SamplePoints,2);
            
            if Sim %Uni
                obj.ufs = zeros(SamplePoints,1);
                obj.Pisom = zeros(SamplePoints,1);
                obj.lr = zeros(SamplePoints,1);       
            else %Bi
                obj.ufsN = zeros(SamplePoints,31);
                obj.ufs = zeros(SamplePoints,1);
                obj.ri = zeros(SamplePoints,1);
                obj.Do = zeros(SamplePoints,1);
                obj.FT = zeros(SamplePoints,1);
            end
        end
        
        function obj = UpdateVectors(obj,i,cs)
            if length(cs.lrNum)==1
                j = 1;
                symVar = cs.lrG;
                numVar = cs.lrNum;
            else
                j = 2;
                symVar = cs.riG;
                numVar = cs.riNum;
            end
            
            lrt = cs.lrNum(j);
            lt = cs.ltNum(j);
            lz = cs.lzNum;
            
            obj.ufs(i,:) = cs.ufs(j);
            obj.stretch(i,:) = [lrt,lt,lz, (lrt*lt*lz)];
            obj.sECM(i,:) = subs(cs.sECM(:,j),symVar,numVar);
            obj.sSMC(i,:) = subs(cs.sSMC(:,j),symVar,numVar);
            obj.sMMy(i,:) = subs(cs.sMMy(:,j),symVar,numVar);
            obj.PMMCU(i,:) = subs([cs.PMM(j),cs.PCU(j)],symVar,numVar);

            if ~isempty(obj.Pisom)
                obj.ufs(i,:) = cs.ufs;
                obj.lr(i,:) = cs.lrNum;  
                
                obj.Pisom(i) = cs.Pisom;
                obj.dMAiMA0(i,:) = [cs.dMArMA0,cs.dMAzMA0];
                obj.LMi(i,:) = [cs.LMr,cs.LMz];
                obj.Lfoi(i,:) = [cs.Lfor,cs.Lfoz];
                obj.eS2(i,:) = [cs.eS2xr,cs.eS2yr,cs.eS2xz,cs.eS2yz];
                obj.x2c(i,:) = [cs.x2cr,cs.x2cz];
                obj.y2c(i,:) = [cs.y2cr,cs.y2cz];
            end
            
            if ~isempty(obj.Do)
                obj.ri(i) = cs.riNum;
            	obj.Do(i) = cs.roNum*2e3;
            	obj.FT(i) = cs.FT*1e3;
                
                obj.dMAiMA0(i,:) = [cs.dMArMA0(2),cs.dMAzMA0];
                obj.LMi(i,:) = [cs.LMr(2),cs.LMz];
                obj.Lfoi(i,:) = [cs.Lfor(2),cs.Lfoz(2)];
                obj.eS2(i,:) = [cs.eS2xr(2),cs.eS2yr(2),cs.eS2xz,cs.eS2yz];
                obj.x2c(i,:) = [cs.x2cr(2),cs.x2cz];
                obj.y2c(i,:) = [cs.y2cr(2),cs.y2cz];
            end
        end
    end
end

