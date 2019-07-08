classdef SimulationVectors < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %Both Simulations
        ufsVec
        stretchVec
        sECMVec
        sSMCVec
        sMMyVec
        PMMCUVec
        
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
        
        function obj = UpdateVectors(obj,i,ufs,stretch,sECM,sSMC,sMMy,PMMCU,Pisom,Do,FT)
            obj.ufsVec(i,:) = ufs;
            obj.stretchVec(i,:) = [stretch, stretch(1)*stretch(2)*stretch(3)];
            obj.sECMVec(i,:) = sECM;
            obj.sSMCVec(i,:) = sSMC;
            obj.sMMyVec(i,:) = sMMy;
            obj.PMMCUVec(i,:) = PMMCU;
            
            if ~isempty(obj.PisomVec)
                obj.PisomVec(i) = Pisom;
            end
            
            if exist('Do','var')
                if ~isempty(obj.DoVec)
                    obj.DoVec(i) = Do;
                    obj.FTVec(i) = FT;
                end
            end
        end
    end
end

