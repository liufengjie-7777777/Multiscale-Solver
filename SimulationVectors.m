classdef SimulationVectors < handle
    %SimulationVectors contains the variables calculated during simulation
    properties
        %Both Simulations
        ufsVec
        stretchVec
        sECMVec
        sSMCVec
        sMMyVec
        PMMCUVec
        
        nAMpVec
        nAMVec
        timeVec
        
        %Uniaxial Simulation
        PisomVec
        
        %Biaxial Simulation
        DoVec
        FTVec
    end
    
    methods
        function obj = InitialVectors(obj,SamplePoints,Sim)
            %Both
            obj.ufsVec = zeros(SamplePoints,1);
            obj.stretchVec = zeros(SamplePoints,4);
            
            obj.sECMVec = zeros(SamplePoints,3);
            obj.sSMCVec = zeros(SamplePoints,3);
            obj.sMMyVec = zeros(SamplePoints,3);
            
            obj.PMMCUVec = zeros(SamplePoints,2);
            
            if Sim %Uni
                obj.PisomVec = zeros(SamplePoints,1);
            else %Bi
                obj.DoVec = zeros(SamplePoints,1);
                obj.FTVec = zeros(SamplePoints,1);
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
            obj.ufsVec(i,:) = cs.ufs(j);
            lr = cs.lrNum(j);
            lt = cs.ltNum(j); 
            lz = cs.lzNum;
            sECM = cs.sECM(:,j);
            sSMC = cs.sSMC(:,j);
            sMMy = cs.sMMy(:,j);
            PMM = cs.PMM(j);
            PCU = cs.PCU(j);
            
            obj.stretchVec(i,:) = [lr,lt,lz, (lr*lt*lz)];
            obj.sECMVec(i,:) = subs(sECM,symVar,numVar);
            obj.sSMCVec(i,:) = subs(sSMC,symVar,numVar);
            obj.sMMyVec(i,:) = subs(sMMy,symVar,numVar);
            obj.PMMCUVec(i,:) = subs([PMM,PCU],symVar,numVar);
            
            if ~isempty(obj.PisomVec)
                obj.PisomVec(i) = cs.Pisom;
            end
            
            if ~isempty(obj.DoVec)
            	obj.DoVec(i) = cs.roNum*2e3;
            	obj.FTVec(i) = cs.FT*1e3;
            end
        end
    end
end

