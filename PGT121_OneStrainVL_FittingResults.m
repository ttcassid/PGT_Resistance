% Visualize the parameter fitting from the BnAb fitting
close all
clear all
%% Choose patient
PatientIndex = 13;

%% Patient specific data
PatientID    = ["6292","9372","6775","2990","2305","5257",...
    "7190","8074","2319","6113","1536","2936",...
    "4236"];
PatientNumber = length( PatientID);  
InitialViralLoad = [55700,0,21040,23930,9660,0,2670,0,44800,190,180,750,350];
%% Load PK parameters
PopulationPKParametersStruct = load('PatientPKParameters_PopulationFit_03Jul23.mat'); % load('PatientPKParameters12022020V2.mat'); % load('PatientPKParameters12022020.mat');
PatientPKParameters = cell2mat(struct2cell(PopulationPKParametersStruct));
%% Load the Dose size
PopulationDoseSize = load('PatientDoseSize12022020V1.mat');
VectorDoseSize = cell2mat(struct2cell(PopulationDoseSize));
%% Load the VL parameters for One Strain without mutation 
disp('One strain model') %
PopulationVLParametersStructOneStrainMutation = load('OneStrainFittingResults4Jul2023.mat'); %  load('OneStrainFittingResults9Jan2023_V2.mat'); %
PatientVLParametersOneStrain = cell2mat(struct2cell(PopulationVLParametersStructOneStrainMutation));
VectorOfMonolixID = PatientVLParametersOneStrain(:,1);

PA.mu =  0;
InitialViralLoadVec = zeros(1,length(PatientNumber));

BICTable = zeros(3,length(PatientNumber));
PatientIndexVec = [1,3:5,7,9:13];
for ii =  1: length(PatientIndexVec);  
    
    PatientIndex = PatientIndexVec(ii);
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
    
    %% Set up solver for viral dynamics model
    if PatientIndex == 9
        X = [1.023738653439307,1.276791464757337,-1.090982543311391,1.169620572768574];
        PA.a = 0;
    elseif PatientIndex == 10
        TotalTime = [0 275];
        X = [1.05530843167287;1.20516254759505;2.70891893663471;1.45481605858392;0.424404677052051]; 
        PA.a = 8/10;  
    elseif PatientIndex == 11
        TotalTime = [0 283];
        X =  [1.00058863648203;1.33678974350726;2.459757842934110;1.22347501265493;0.269059035649851];  
        PA.a = 5/10; 
    else
        PA.a = 0;
        X = PatientVLParametersOneStrain(:,PatientIndex); %
        TotalTime = [0 56];
    end
    
    %% Model specific parameters
    PA.dT = 0.01; 
    PA.f = 0.05;
    PA.c = 23;
    
    %Patient specific parameters
    PA.omega =   X(1);  
    %Parameters for virus
    PA.p = 10^(X(2));
    PA.alphaS = 10^(X(3)) ; 
    PA.tau = X(4);  
    
    % Initial conditions
    TIC = 649500;
    VIC = InitialVL;
    IIC = PA.c.*VIC./PA.p;
    A1IC = 0;
    A2IC = 0;

    IC = [TIC,IIC,VIC,A1IC,A2IC];
    
    %Parameters dependent on IC
    PA.delta = 1.5/(IIC)^(PA.omega-1);
    PA.beta = PA.delta.*(IIC)^(PA.omega)./(PA.f.*(VIC).*TIC);
    
    PA.lambda = (PA.dT+PA.beta.*VIC ).*TIC;
    
    [solPGT121OneStrainMutation] = PGT121_OneStrainMutation_Solver(TotalTime,IC,PA); 
    
    %% Plotting solPGT121TwoStrainMutation with PK
    Fig = figure(PatientIndex);  
    set(Fig,'defaultAxesColorOrder',[ [33,102,172]./255 ; [227,26,28]./255]);
    yyaxis left
    hold on
    g5 = plot(TotalTime,[40 40], 'LineWidth',0.75,'Color',[67,147,195]/255,'LineStyle',':');
    hold on
    g3 =  plot(solPGT121OneStrainMutation.x(:),max(solPGT121OneStrainMutation.y(3,:),0),'LineWidth',2.25,'Color', [171,217,233]/255,'LineStyle','-'); %grey
    hold on
    set(gca,'yscale','log')
    ylim([10^(-2),10^6])
    xlim([0,inf])
    
    ylabel('Viral Load (copies/mL) ','FontSize',15); 
    yyaxis right
    set(gca,'yscale','log')
    g31 = plot(solPGT121OneStrainMutation.x(:),solPGT121OneStrainMutation.y(4,:),'LineWidth',2,'Color',[227,26,28]./256,'LineStyle','-');
    hold on
    ylim([10^(-2), inf])
    xlim([0,inf])
    xlabel('Time (days)','FontSize',15)
    ylabel(' PGT-121 (mg/mL) ','FontSize',15);  
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'r';
    title(PatientID(PatientIndex))

    
end

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



%% Single strain Model in
function [sol] = PGT121_OneStrainMutation_Solver(TotalTime,IC,PA) %DDE model without therapy
opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2);
sol = ode45(@BnAbVLDynamics,TotalTime,IC,opts);
    function dydt = BnAbVLDynamics(t,y)
        dydt(1) = PA.lambda-PA.dT.*y(1)-InfectivityLoss(y(4),t,PA).*PA.beta.*y(1).*y(3); %DE for target
        dydt(2) = PA.a + InfectivityLoss(y(4),t,PA).*PA.f.*PA.beta.*y(1).*y(3)-PA.delta.*y(2)^(PA.omega); %Differential equation for infected
        dydt(3) = PA.p.*y(2)-PA.c.*y(3); %Differential equation for virus
        dydt(4) =  Dose(t,PA) - PA.k12.*y(4)+ PA.k21.*(PA.Vol2./PA.Vol1).*y(5); %Differential equation for A1
        dydt(5) =  PA.k12.*(PA.Vol1./PA.Vol2).*y(4) - PA.k21.*y(5)-PA.k0.*y(5) ; %Differential equation for A2
        dydt = dydt';
    end

end
