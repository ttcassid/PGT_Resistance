% BnABSensitivityAnalysis script
close all
clear all
%% Load PK parameters
PopulationPKParametersStruct = load('PatientPKParameters_PopulationFit_03Jul23.mat');
PatientPKParameters = cell2mat(struct2cell(PopulationPKParametersStruct));
%% Load the Dose size
PopulationDoseSize = load('PatientDoseSize12022020V1.mat');
VectorDoseSize = cell2mat(struct2cell(PopulationDoseSize));

%% Load the VL parameters for single strain model
PopulationVLParametersStructOneStrainMutation = load('OneStrainFittingResults4Jul2023.mat'); 
PatientVLParametersOneStrain = cell2mat(struct2cell(PopulationVLParametersStructOneStrainMutation)); 

%% Model generic parameters
PA.dT = 0.01;
PA.f = 0.05;
PA.dV = 23;
InitialViralLoad = [55700,0,21040,23930,9660,0,2670,0,17570,190,180,750,350];

%% Initialize storage
ParameterNames = {'omega' 'p'  'alpha' 'tau'}; %    %The name of parameters to be tested
NumberOfParameters = length(ParameterNames);
PatientIndexVec =   [1,3:5,7,9:13]; 

VLNadirLowerParameter = zeros(length(PatientIndexVec),NumberOfParameters);
VLNadirUpperParameter = zeros(length(PatientIndexVec),NumberOfParameters);
ReboundTimeLowerParameter = zeros(length(PatientIndexVec),NumberOfParameters);
ReboundTimeUpperParameter = zeros(length(PatientIndexVec),NumberOfParameters);
VLNadirBaseline  = zeros(1,length(PatientIndexVec));
TimeToReboundBaseline = zeros(1, length(PatientIndexVec));
BetaVec = zeros(1, length(PatientIndexVec));
for ii =  1: length(PatientIndexVec); 
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
X = PatientVLParametersOneStrain(:,PatientIndex) ;
PA.a = 0; 
%Parameters for infected cells
PA.omega = X(1);
%Parameters for virus
PA.p = 10^(X(2));
PA.alphaS = 10^(X(3)) ;  
PA.tau = X(4); % delay for antibody mediated loss of infectivity

% Initial conditions
TIC = 649500;  
VIC = InitialVL(1);
IIC = PA.dV.*VIC./PA.p; 
A1IC = 0;
A2IC = 0;
TotalTime = [0 500];
IC = [TIC,IIC,VIC,A1IC,A2IC];

%Parameters dependent on IC
PA.dI = 1.5/IIC^(PA.omega-1);
PA.beta = PA.dI.*IIC^(PA.omega)./(PA.f.*VIC.*TIC);
PA.lambda = (PA.dT+PA.beta.*VIC).*TIC;
BetaVec(ii) = PA.beta;

%% Set up sampling procedure
ParameterNames = {'omega' 'p'  'alphaS' 'tau'};   %The name of parameters to be tested
ParameterHomeostasis = [PA.omega ;PA.p ; PA.alphaS; PA.tau]; 
NumberOfParameters = length(ParameterNames);
PA.ReboundVL = 0.75.*InitialVL(1);

% Create lower and upper paraemter values
P = 0.1; %Percentage of parameter change
LowerBoundVector = (1-P).*ParameterHomeostasis;
LowerBoundVector(1) = max( (1-P).*PA.omega,1); 
UpperBoundVector = (1+P).*ParameterHomeostasis; 

if PatientIndex == 10;
    PA.a = 0.8;
elseif PatientIndex == 11;
    PA.a = 0.5;
else
    PA.a = 0;
end
%% Simulate for baseline 
[solBnAbVLOneStrain] = BnAbVLDynamicsSolverOneStrain(TotalTime,IC,PA);
% % Measurements
VLNadirBaseline(ii) = min(solBnAbVLOneStrain.y(3,:));
TimeToReboundBaseline(ii) = solBnAbVLOneStrain.xe(1);

% Set up parameter variation

for kk = 1:NumberOfParameters;
    %Reset parameters to healthy levels
    for jj = 1:kk
        PA.(ParameterNames{jj}) = ParameterHomeostasis(jj);
    end
    %% Simulate model for lower parameter value
    PA.(ParameterNames{kk}) = LowerBoundVector(kk);
    
    % calculate parameters and initial conditions
   % Initial conditions
   TIC = 649500;  % PA.lambda./PA.dT;
   VIC = InitialVL;
   IIC = PA.dV.*VIC./PA.p;
   A1IC = 0;
   A2IC = 0;
   IC = [TIC,IIC,VIC,A1IC,A2IC];
   
   %Parameters dependent on IC
   PA.dI = 1.5/IIC^(PA.omega-1);
   PA.beta = PA.dI.*IIC^(PA.omega)./(PA.f.*VIC.*TIC);
   PA.lambda = (PA.dT+PA.beta.*VIC).*TIC; 
   
    %Simulate for lower parameter 
    [solBnAbVLOneStrainLower] =  BnAbVLDynamicsSolverOneStrain(TotalTime,IC,PA);
    % % Measurements
    VLNadirLowerParameter(ii,kk) = min( solBnAbVLOneStrainLower.y(3,:) );
    if ~isempty(solBnAbVLOneStrainLower.xe)
        ReboundTimeLowerParameter(ii,kk) = solBnAbVLOneStrainLower.xe(end);
    else
        ReboundTimeLowerParameter(ii,kk) = inf;
    end
    
    %% Simulate model for upper parameter value
    PA.(ParameterNames{kk}) = UpperBoundVector(kk);
    
       % Initial conditions
   TIC = 649500;  % PA.lambda./PA.dT;
   VIC = InitialVL;
   IIC = PA.dV.*VIC./PA.p;
   A1IC = 0;
   A2IC = 0;
   IC = [TIC,IIC,VIC,A1IC,A2IC];
   
   %Parameters dependent on IC
   PA.dI = 1.5/IIC^(PA.omega-1);
   PA.beta = PA.dI.*IIC^(PA.omega)./(PA.f.*VIC.*TIC);
   PA.lambda = (PA.dT+PA.beta.*VIC).*TIC; 
   
    
    %Simulate for upper parameters 
    [solBnAbVLOneStrainUpper] = BnAbVLDynamicsSolverOneStrain(TotalTime,IC,PA); 
    % % Measurements
    VLNadirUpperParameter(ii,kk) = min( solBnAbVLOneStrainUpper.y(3,:) );
    
    if ~isempty(solBnAbVLOneStrainUpper.xe)
        ReboundTimeUpperParameter(ii,kk) = solBnAbVLOneStrainUpper.xe(1);
    else
        ReboundTimeUpperParameter(ii,kk) = inf;
    end

end

end
%% Aggregate data plotting
AggregateReboundTimeUpperParameterSensitivity = zeros(1,NumberOfParameters );
AggregateReboundTimeLowerParameterSensitivity = zeros(1,NumberOfParameters );
AggregateVLNadirUpperParameterSensitivity = zeros(1,NumberOfParameters );
AggregateVLNadirLowerParameterSensitivity = zeros(1,NumberOfParameters );
for kk = 1:NumberOfParameters % for each parameter
        AggregateVLNadirLowerParameterSensitivity(kk) = median( [VLNadirLowerParameter(:,kk)./VLNadirBaseline(:)]);
    AggregateVLNadirUpperParameterSensitivity(kk) = median( [VLNadirUpperParameter(:,kk)./VLNadirBaseline(:)]);
    AggregateReboundTimeLowerParameterSensitivity(kk) = median( [ReboundTimeLowerParameter(:,kk)./TimeToReboundBaseline(:)]);
    AggregateReboundTimeUpperParameterSensitivity(kk) = median( [ReboundTimeUpperParameter(:,kk)./TimeToReboundBaseline(:)]);
end

Fig1 = figure(1);
xvalues = {'-10%','+10%'};
yvalues = ParameterNames;  
VLNadirHM = heatmap(xvalues,yvalues,100.*[AggregateVLNadirLowerParameterSensitivity',AggregateVLNadirUpperParameterSensitivity'],'ColorbarVisible','off') ; 
 
%Time to rebound duration
Fig2 = figure(2 );
xvalues = {'-10%','+10%'};
yvalues = ParameterNames;  
ReboundDelayHM = heatmap(xvalues,yvalues,100.*[AggregateReboundTimeLowerParameterSensitivity',AggregateReboundTimeUpperParameterSensitivity'],'ColorbarVisible','off') ; 

[a,b] = corr(PatientVLParametersOneStrain(:,PatientIndexVec)')

[BetaR,BetaPval] = corr(BetaVec',PatientVLParametersOneStrain(3,PatientIndexVec)')

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


function [sol] = BnAbVLDynamicsSolverOneStrain(TotalTime,IC,PA) 
opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2,'Events',@TimeToRebound);
sol = ode45(@BnAbVLDynamics,TotalTime,IC,opts);
function dydt = BnAbVLDynamics(t,y)  
    dydt(1) = PA.lambda-PA.dT.*y(1)-PA.beta.*InfectivityLoss(y(4),t,PA).*y(1).*y(3); %DE for target
    dydt(2) = PA.a + InfectivityLoss(y(4),t,PA).*PA.f.*PA.beta.*y(1).*y(3)-PA.dI.*y(2)^(PA.omega); %Differential equation for infected
    dydt(3) = PA.p.*y(2)-PA.dV.*y(3); %Differential equation for virus
    dydt(4) =  Dose(t,PA) - PA.k12.*y(4)+ PA.k21.*(PA.Vol2./PA.Vol1).*y(5); %Differential equation for A1
    dydt(5) =  PA.k12.*(PA.Vol1./PA.Vol2).*y(4) - PA.k21.*y(5)-PA.k0.*y(5) ; %Differential equation for A2
    dydt = dydt';
end

function [value,isterminal,direction] = TimeToRebound(t,y)  

value(1)= y(3) - PA.ReboundVL ;  %What we are setting to 0: 
isterminal(1) = 1;   % Continue the integration
direction(1) = 1;   % increasing direction only

end
end







