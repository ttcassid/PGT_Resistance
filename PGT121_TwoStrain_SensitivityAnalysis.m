close all
clear all
%% Model generic parameters
PA.dT = 0.01;
PA.f = 0.05;
PA.c = 23;
%% Load PK parameters
PopulationPKParametersStruct = load('PatientPKParameters_PopulationFit_03Jul23.mat'); %  load('PatientPKParameters12022020V2.mat'); % load('PatientPKParameters12022020.mat');
PatientPKParameters = cell2mat(struct2cell(PopulationPKParametersStruct));
%% Load the Dose size
PopulationDoseSize = load('PatientDoseSize12022020V1.mat');
VectorDoseSize = cell2mat(struct2cell(PopulationDoseSize));

%% Load the VL parameters for Two Strain w/mutation
% disp('Two strain model with Mutation')
% PopulationVLParametersStructTwoStrainMutation = load('TwoStrainMutationFittingResults4Jul2023'); % load('TwoStrainMutationFittingResults9Jan2023.mat');
% PatientVLParametersTwoStrain = cell2mat(struct2cell(PopulationVLParametersStructTwoStrainMutation));
%    PA.mu =  2.16e-5 ;

disp('Two strain model without mutation')
PopulationVLParametersStructTwoStrainMutation = load('TwoStrainFittingResults4Jul2023'); load('TwoStrainFittingResults9Jan2023.mat');
PatientVLParametersTwoStrain = cell2mat(struct2cell(PopulationVLParametersStructTwoStrainMutation));
PA.mu =  0 ;

InitialViralLoad = [55700,0,21040,23930,9660,0,2670,0,17570,190,180,750,350];
PatientIndexVec =   [1,3:5,7,9:13];
%% Initialize storage
ParameterNames = {'omega' 'p'  'alphaS' 'tau' 'rho' 'alphaR'};   %The name of parameters to be tested
NumberOfParameters = length(ParameterNames);

VLNadirLowerParameter = zeros(length(PatientIndexVec),NumberOfParameters);
VLNadirUpperParameter = zeros(length(PatientIndexVec),NumberOfParameters);
ReboundTimeLowerParameter = zeros(length(PatientIndexVec),NumberOfParameters);
ReboundTimeUpperParameter = zeros(length(PatientIndexVec),NumberOfParameters);
VLNadirBaseline  = zeros(1,length(PatientIndexVec));
TimeToReboundBaseline = zeros(1, length(PatientIndexVec));
ReSensitiveTimeBaseline = zeros(1, length(PatientIndexVec));
ReSensitiveTimeLowerParameter = zeros(length(PatientIndexVec),NumberOfParameters);
ReSensitiveTimeUpperParameter = zeros(length(PatientIndexVec),NumberOfParameters);

TimeNadirBaseline = zeros(1, length(PatientIndexVec));
TargetBaseline = zeros(1, length(PatientIndexVec));
AntibodyBaseline = zeros(1, length(PatientIndexVec));
AntibodyRatio = zeros(1, length(PatientIndexVec));
AntibodyRatioResistant = zeros(1, length(PatientIndexVec));

for ii =   1: length(PatientIndexVec);
    PatientIndex = PatientIndexVec(ii) 
    InitialVL = InitialViralLoad(PatientIndex);
    
    Y =  PatientPKParameters(:,PatientIndex);
    W = VectorDoseSize(PatientIndex);
    %% PK Parameters from fitting
    PA.k12 = Y(1); %[1/day] Blood to tissue rate
    PA.k0 = Y(2); %[1/day] Clearance rate
    PA.k21 = Y(3); %[1/day] Tissue to blood
    %% Set the fixed PK model parameters
    PA.Vol1 = 3000; %[mL of plasma distribution]
    PA.Vol2 = PA.k12.*PA.Vol1./PA.k21; % [Liters, imposing k_12 x Vol1 = k_21 x Vol2]
    PA.Tinf = 1/24; % [day] Infusion time
    % Load the Dose size
    PA.Dose = W;
    
    %% Load patient specific VL data
    X = PatientVLParametersTwoStrain(:,PatientIndex) ;
    if PA.mu == 0
        if PatientIndex == 10;
            PA.a = 0.8;
        elseif PatientIndex == 11;
            PA.a = 0.65;
        else
            PA.a = 0;
        end
    else
        if PatientIndex == 10;
            PA.a = 0.4;
        elseif PatientIndex == 11;
            PA.a = 0.7;
        else
            PA.a = 0;
        end
    end
    %Parameters for infected cells
    PA.omega =   X(1);
    %Parameters for virus
    PA.p = 10^(X(2));
    PA.alphaS = 10^(X(3)) ;
    PA.tau = X(4); % delay for antibody mediated loss of infectivity
    PA.rho = X(5);
    PA.alphaR = 10^(X(6)) ;
    
    % Initial conditions
    TIC = 649500;
    VIC =(1-max(PA.rho,0)).*InitialVL;
    ResistantVIC =max(PA.rho,0).*InitialVL;
    
    IIC = PA.c.*VIC./PA.p;
    RIC = PA.c.*ResistantVIC./PA.p;
    A1IC = 0;
    A2IC = 0;
    TotalTime = [0 700];
    IC = [TIC,IIC,VIC,A1IC,A2IC,RIC,ResistantVIC];
    
    %Parameters dependent on IC
    PA.delta = 1.5/(RIC+IIC)^(PA.omega-1);
    if X(5) == 1
        PA.beta = 0;
    else
        PA.beta = PA.delta.*(IIC)^(PA.omega)./(PA.f.*(VIC).*TIC);
    end
    
    if X(5) == 0
        PA.betaR = 0;
    else
        PA.betaR = PA.delta.*(RIC)^(PA.omega)./(PA.f.*(ResistantVIC).*TIC);
    end
    
    PA.lambda = (PA.dT+PA.beta.*VIC+PA.betaR.*ResistantVIC ).*TIC;
    %% Set up sampling procedure
    ParameterNames = {'omega' 'p'  'alphaS' 'tau'  'rho' 'alphaR'};
    ParameterHomeostasis = [PA.omega ;PA.p ; PA.alphaS; PA.tau; PA.rho; PA.alphaR ];
    NumberOfParameters = length(ParameterNames);
    PA.ReboundVL = 0.75.*InitialVL;
    
    % Create lower and upper paraemter values
    P = 0.10; %Percentage of parameter change
    LowerBoundVector = (1-P).*ParameterHomeostasis;
    LowerBoundVector(1) = max( (1-P).*PA.omega,1);
    LowerBoundVector(3) =  PA.alphaS.*(1-P);
    LowerBoundVector(6) =  PA.alphaR.*(1-P);
    UpperBoundVector = (1+P).*ParameterHomeostasis;
    UpperBoundVector(3) =   PA.alphaS.*(1+P);
    UpperBoundVector(5) = min( (1+P).*PA.rho,0.5);
    UpperBoundVector(6) =  PA.alphaR.*(1+P);
    
    %% Simulate for baseline
    [solBnAbVLTwoStrain] = BnAbVLDynamicsSolverTwoStrain(TotalTime,IC,PA);
    % % Measurements
    [VLNadirBaseline(ii),TimeNadirBaseline(ii)] = min(solBnAbVLTwoStrain.y(3,:)+solBnAbVLTwoStrain.y(7,:));
    TargetBaseline(ii) = solBnAbVLTwoStrain.y(1,TimeNadirBaseline(ii));
    AntibodyBaseline(ii) = solBnAbVLTwoStrain.y(4,TimeNadirBaseline(ii));
    
    VLResistantNadirBaseline(ii) = min(solBnAbVLTwoStrain.y(7,:));
    ResistantRatio(ii) = -log( VLResistantNadirBaseline(ii)./solBnAbVLTwoStrain.y(7,1) )/log(10) ;
    
    AntibodyRatio(ii) = AntibodyBaseline(ii)./PA.alphaS;
    AntibodyRatioResistant(ii) = AntibodyBaseline(ii)./PA.alphaR;
    
    TimeToReboundBaseline(ii) = solBnAbVLTwoStrain.xe(1,1);
    if length(solBnAbVLTwoStrain.xe(1,:) ) > 1
        ReSensitiveTimeBaseline(ii) =  solBnAbVLTwoStrain.xe(1,end);
    else
        ReSensitiveTimeBaseline(ii) = NaN;
    end
    
    % Set up parameter variation
    for kk = 1:NumberOfParameters;
        %Reset parameters to healthy levels
        for jj = 1:kk
            PA.(ParameterNames{jj}) = ParameterHomeostasis(jj);
        end
        %% Simulate model for lower parameter value
        PA.(ParameterNames{kk}) = LowerBoundVector(kk);
        
        % calculate parameters and initial conditions
        VIC =(1-max(PA.rho,0)).*InitialVL(1);
        ResistantVIC =max(PA.rho,0).*InitialVL(1);
        IIC = PA.c.*VIC./PA.p;
        RIC = PA.c.*ResistantVIC./PA.p;
        A1IC = 0;
        A2IC = 0;
        
        IC = [TIC,IIC,VIC,A1IC,A2IC,RIC,ResistantVIC];
        %Parameters dependent on IC
        PA.delta = 1.5/(RIC+IIC)^(PA.omega-1);
        if PA.rho == 1
            PA.beta = 0;
        else
            PA.beta = PA.delta.*(IIC)^(PA.omega)./(PA.f.*(VIC).*TIC);
        end
        if PA.rho == 0
            PA.betaR = 0;
        else
            PA.betaR = PA.delta.*(RIC)^(PA.omega)./(PA.f.*(ResistantVIC).*TIC);
        end
        
        PA.lambda = (PA.dT+PA.beta.*VIC+PA.betaR.*ResistantVIC ).*TIC;
        %Simulate for lower parameter
        [solBnAbVLTwoStrainLower] = BnAbVLDynamicsSolverTwoStrain(TotalTime,IC,PA); 
        % % Measurements
        VLNadirLowerParameter(ii,kk) = min(solBnAbVLTwoStrainLower.y(3,:)+solBnAbVLTwoStrainLower.y(7,:));
        if ~isempty(solBnAbVLTwoStrainLower.xe(:,1))
            ReboundTimeLowerParameter(ii,kk) = solBnAbVLTwoStrainLower.xe(1,1);
        else
            ReboundTimeLowerParameter(ii,kk) = inf;
        end
        
        if length(solBnAbVLTwoStrainLower.xe(1,:)) > 1  
            ReSensitiveTimeLowerParameter(ii,kk) = solBnAbVLTwoStrainLower.xe(1,end);
        else
            ReSensitiveTimeLowerParameter(ii,kk) = NaN;
        end
        
        %% Simulate model for upper parameter value
        PA.(ParameterNames{kk}) = UpperBoundVector(kk);
        
        % calculate parameters and initial conditions
        VIC =(1-max(PA.rho,0)).*InitialVL(1);
        ResistantVIC =max(PA.rho,0).*InitialVL(1);
        IIC = PA.c.*VIC./PA.p;
        RIC = PA.c.*ResistantVIC./PA.p;
        A1IC = 0;
        A2IC = 0;
        IC = [TIC,IIC,VIC,A1IC,A2IC,RIC,ResistantVIC];
        %Parameters dependent on IC
        PA.delta = 1.5/(RIC+IIC)^(PA.omega-1);
        if PA.rho == 1
            PA.beta = 0;
        else
            PA.beta = PA.delta.*(IIC)^(PA.omega)./(PA.f.*(VIC).*TIC);
        end
        if PA.rho == 0
            PA.betaR = 0;
        else
            PA.betaR = PA.delta.*(RIC)^(PA.omega)./(PA.f.*(ResistantVIC).*TIC);
        end
        
        PA.lambda = (PA.dT+PA.beta.*VIC+PA.betaR.*ResistantVIC ).*TIC;
        
        %Simulate for upper parameters
        [solBnAbVLTwoStrainUpper] = BnAbVLDynamicsSolverTwoStrain(TotalTime,IC,PA);
        % % Measurements
        VLNadirUpperParameter(ii,kk) = min(solBnAbVLTwoStrainUpper.y(3,:)+solBnAbVLTwoStrainUpper.y(7,:));
        
        if ~isempty(solBnAbVLTwoStrainUpper.xe)
            ReboundTimeUpperParameter(ii,kk) = solBnAbVLTwoStrainUpper.xe(1,1);
        else
            ReboundTimeUpperParameter(ii,kk) = inf;
        end
        
        if length(solBnAbVLTwoStrainUpper.xe) > 1 
            ReSensitiveTimeUpperParameter(ii,kk) = solBnAbVLTwoStrainUpper.xe(1,end);
        else
            ReSensitiveTimeUpperParameter(ii,kk) = NaN;
        end
        
    end
    
end
%% Aggregate data plotting
AggregateReboundTimeUpperParameterSensitivity = zeros(1,NumberOfParameters );
AggregateReboundTimeLowerParameterSensitivity = zeros(1,NumberOfParameters );
AggregateVLNadirUpperParameterSensitivity = zeros(1,NumberOfParameters );
AggregateVLNadirLowerParameterSensitivity = zeros(1,NumberOfParameters );
AggregateReSensitiseTimeUpperParameterSensitivity = zeros(1,NumberOfParameters );
AggregateReSensitiseTimeLowerParameterSensitivity = zeros(1,NumberOfParameters );
for kk = 1:NumberOfParameters % for each parameter
    AggregateVLNadirLowerParameterSensitivity(kk) = nanmedian( [VLNadirLowerParameter(:,kk)./VLNadirBaseline(:)]);
    AggregateVLNadirUpperParameterSensitivity(kk) = nanmedian( [VLNadirUpperParameter(:,kk)./VLNadirBaseline(:)]);
    AggregateReboundTimeLowerParameterSensitivity(kk) = nanmedian( [ReboundTimeLowerParameter(:,kk)./TimeToReboundBaseline(:)]);
    AggregateReboundTimeUpperParameterSensitivity(kk) = nanmedian( [ReboundTimeUpperParameter(:,kk)./TimeToReboundBaseline(:)]);
    AggregateReSensitiseTimeUpperParameterSensitivity(kk) =  nanmedian( [ReSensitiveTimeUpperParameter(:,kk)./ReSensitiveTimeBaseline(:)]);
    AggregateReSensitiseTimeLowerParameterSensitivity(kk) =  nanmedian( [ReSensitiveTimeLowerParameter(:,kk)./ReSensitiveTimeBaseline(:)]);
end

disp('Correlation between alphaS and resensitize time')
[BaselineAlphaMatrix,BaselineAlphaPVal] = corr([PatientVLParametersTwoStrain(3,PatientIndexVec)',ReSensitiveTimeBaseline'],'Type','Pearson','tail','right')

MedianTimeToRebound = median(TimeToReboundBaseline)
MedianTimeToResensitize = median(ReSensitiveTimeBaseline)

Fig1 = figure(1);
xvalues = {'-10%','+10%'};
yvalues = ParameterNames;  
VLNadirHM = heatmap(xvalues,yvalues,100.*[AggregateVLNadirLowerParameterSensitivity',AggregateVLNadirUpperParameterSensitivity'],'ColorbarVisible','off') ; 

% %Time to rebound duration
Fig2 = figure(2);
xvalues = {'-10%','+10%'};
yvalues = ParameterNames;  
ReboundDelayHM = heatmap(xvalues,yvalues,100.*[AggregateReboundTimeLowerParameterSensitivity',AggregateReboundTimeUpperParameterSensitivity'],'ColorbarVisible','off') ;

Fig3 = figure(3);
xvalues = {'-10%','+10%'};
yvalues = ParameterNames; 
ReSensitiseDelayHM = heatmap(xvalues,yvalues,100.*[AggregateReSensitiseTimeLowerParameterSensitivity',AggregateReSensitiseTimeUpperParameterSensitivity'],'ColorbarVisible','off') ; %,'Title','Tumour Doubling Time (% of normal)');



function DoseVal = Dose(t,PA)
if 0< t && t < PA.Tinf
    DoseVal = PA.Dose;
else
    DoseVal = 0;
end

end

function Loss = InfectivityLoss(Ab,t,PA)
if t < PA.tau
    Loss = 1;
else
    Loss = 1./(1+PA.alphaS.*Ab);
end
end


function Loss = InfectivityResistantLoss(Ab,t,PA)
if t < PA.tau
    Loss =1;
else
    Loss = 1./(1+PA.alphaR.*Ab);
end
end


function [sol] = BnAbVLDynamicsSolverTwoStrain(TotalTime,IC,PA) %ODE model without therapy
opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2,'Events',@TimeToRebound);
sol = ode45(@BnAbVLDynamics,TotalTime,IC,opts);
    function dydt = BnAbVLDynamics(t,y) ;
        dydt(1) = PA.lambda-PA.dT.*y(1)-InfectivityLoss(y(4),t,PA).*PA.beta.*y(1).*y(3)-InfectivityResistantLoss(y(4),t,PA).*PA.betaR.*y(1).*y(7); %DE for target
        dydt(2) = PA.a + InfectivityLoss(y(4),t,PA).*PA.f.*PA.beta.*y(1).*y(3).*(1-PA.mu)-PA.delta.*y(2)^(PA.omega); %Differential equation for infected
        dydt(3) = PA.p.*y(2)-PA.c.*y(3); %Differential equation for virus
        dydt(4) =  Dose(t,PA) - PA.k12.*y(4)+ PA.k21.*(PA.Vol2./PA.Vol1).*y(5); %Differential equation for A1
        dydt(5) =  PA.k12.*(PA.Vol1./PA.Vol2).*y(4) - PA.k21.*y(5)-PA.k0.*y(5) ; %Differential equation for A2
        dydt(6) =  InfectivityLoss(y(4),t,PA).*PA.f.*PA.beta.*y(1).*y(3).*PA.mu+ InfectivityResistantLoss(y(4),t,PA).*PA.f.*PA.betaR.*y(1).*y(7)-PA.delta.*y(6)^(PA.omega);  %Differential equation for resistant infected
        dydt(7) =   PA.p.*y(6)-PA.c.*y(7); %Differential equation for resistant virus
        dydt = dydt';
    end
    function [value,isterminal,direction] = TimeToRebound(t,y) 
        
        Rebound = y(3)+y(7) - PA.ReboundVL ;
        isterminal(1) = 0;   % Continue the integration
        direction(1) = 1;   % increasing direction only
        
        % value(2) = y(7)./( y(3)+y(7) ) - PA.rho;
        Sensitive =  y(7)./( y(3)+y(7) ) - 0.5; % PA.rho;
        isterminal(2) = 1;   % End the integration
        direction(2) = -1;   % decreasing direction only 
        
        value = [Rebound,Sensitive]; %,Nadir];
        
    end

end
