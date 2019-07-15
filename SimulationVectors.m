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
            obj.ufsVec(i,:) = cs.ufs;
            obj.stretchVec(i,:) = [cs.lrNum,cs.ltNum,cs.lzNum, (cs.lrNum*cs.ltNum*cs.lzNum)];
            obj.sECMVec(i,:) = subs(cs.sECM,[cs.riG,cs.lrG],[cs.riNum,cs.lrNum]);
            obj.sSMCVec(i,:) = subs(cs.sSMC,[cs.riG,cs.lrG],[cs.riNum,cs.lrNum]);
            obj.sMMyVec(i,:) = subs(cs.sMMy,[cs.riG,cs.lrG],[cs.riNum,cs.lrNum]);
            obj.PMMCUVec(i,:) = subs([cs.PMM,cs.PCU],[cs.riG,cs.lrG],[cs.riNum,cs.lrNum]);
            
            if ~isempty(obj.PisomVec)
                obj.PisomVec(i) = cs.Pisom;
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

