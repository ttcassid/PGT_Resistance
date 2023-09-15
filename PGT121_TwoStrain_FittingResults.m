% Visualize the parameter fitting from the BnAb fitting
close all
clear all
%% Choose patient 
 PatientID    = ["6292","9372","6775","2990","2305","5257",...
                  "7190","8074","2319","6113","1536","2936",...
                  "4236"];
 PatientNumber = length( PatientID);     
 InitialViralLoad = [55700,0,21040,23930,9660,0,2670,0,44800,190,180,750,350];
%% Load PK parameters PatientPKParameters_PopulationFit_03Jul23
PopulationPKParametersStruct = load('PatientPKParameters_PopulationFit_03Jul23.mat'); %  load('PatientPKParameters12022020V2.mat'); % load('PatientPKParameters12022020.mat');
PatientPKParameters = cell2mat(struct2cell(PopulationPKParametersStruct));
%% Load the Dose size
PopulationDoseSize = load('PatientDoseSize12022020V1.mat');
VectorDoseSize = cell2mat(struct2cell(PopulationDoseSize));
%% Load the VL parameters for Two Strain w/mutation 
disp('Two strain model with Mutation')
PopulationVLParametersStructTwoStrainMutation = load('TwoStrainMutationFittingResults4Jul2023'); % load('TwoStrainMutationFittingResults9Jan2023.mat');
PatientVLParametersTwoStrain = cell2mat(struct2cell(PopulationVLParametersStructTwoStrainMutation));
   PA.mu =  2.16e-5 ;

% disp('Two strain model without mutation')
% PopulationVLParametersStructTwoStrainMutation = load('TwoStrainFittingResults4Jul2023'); load('TwoStrainFittingResults9Jan2023.mat');
% PatientVLParametersTwoStrain = cell2mat(struct2cell(PopulationVLParametersStructTwoStrainMutation));
% PA.mu =  0 ;

CalculatedParameters = zeros(length(PatientNumber),4);
PatientParameters = zeros(length(PatientNumber),6);


PatientIndexVec =  [1,3:5,7,9:13]; 
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
    %% Set up solver for viral load
    if PatientIndex == 9
        %% Load patient specific VL data
        X =  PatientVLParametersTwoStrain(:,PatientIndex);
        PA.a = 0;
    elseif PatientIndex == 10
        TotalTime = [0 275];
        %% Load patient specific VL data
        if PA.mu == 0
            X =  [1.0355019128306;1.71651810054430;2.06138505268831023;1.39403999930103;0.0413635065483600;0.4423030771614954;0.864828824313235];
            PA.a = 0.8;
        else
            %% w/  mutation
            X = [1.0345019128306;1.71451810054430;1.980124205268831023;1.42703999930103;0.040710635065483600;0.913030771614954;0.634828824313235];
            PA.a =   0.4;
            
        end
    elseif PatientIndex == 11
        TotalTime = [0 283];
        %% No mutation
        if PA.mu == 0
            X = [1.0082428208151;1.46382009362400907;2.3175499856270271655;1.166302696461533;0.000718253585773;0.445001230825949411;0.608416712494084];
            PA.a = 0.65;
        else
            %% Mutation 
            X = [1.0082828208151;1.485009362400907;2.258599856270271655;1.176302696461533;0.00093789310665773;1.77230825949411;0.608416712494084];
            PA.a = 0.7;
        end
    else
        TotalTime = [0 56];
        PA.a = 0;
        X = PatientVLParametersTwoStrain(:,PatientIndex) ;
        
    end
    %% Model specific parameters
    PA.dT = 0.01; 
    PA.f = 0.05;
    PA.c = 23;
    
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
    VIC =(1-max(PA.rho ,0)).*InitialVL(1);
    ResistantVIC = max(PA.rho,0).*InitialVL(1);
    IIC = PA.c.*VIC./PA.p;
    RIC = PA.c.*ResistantVIC./PA.p;
    A1IC = 0;
    A2IC = 0; 
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
    
    PA.ReactivationBound = [1e-2*PA.delta.*(PA.c.*40./PA.p)^(PA.omega),PA.delta.*(PA.c.*40./PA.p)^(PA.omega)];
    PatientParameters(PatientIndex,:) = [PA.omega;PA.p;PA.alphaS;PA.tau;PA.rho;PA.alphaR]';
    CalculatedParameters(PatientIndex,:) = [PA.beta;PA.betaR;log10(VIC+ResistantVIC);log10(RIC+IIC)]';

    [solPGT121TwoStrainMutation] = PGT121_TwoStrainMutation_Solver(TotalTime,IC,PA);
    
    %% Plotting solPGT121TwoStrainMutation with PK
    
    Fig = figure(PatientIndex); 
    set(Fig,'defaultAxesColorOrder',[ [33,102,172]./255 ; [227,26,28]./255]);
    yyaxis left
    hold on
    g5 = plot(TotalTime,[40 40], 'LineWidth',0.75,'Color',[67,147,195]/255,'LineStyle',':');
    hold on
    g3 =  plot(solPGT121TwoStrainMutation.x(:),max(solPGT121TwoStrainMutation.y(3,:)+solPGT121TwoStrainMutation.y(7,:),0),'LineWidth',2.25,'Color', [171,217,233]/255,'LineStyle','-'); %grey
    hold on
    g4 = plot(solPGT121TwoStrainMutation.x(:),max(solPGT121TwoStrainMutation.y(7,:),0),'LineWidth',1.75,'Color',[239,138,98]/256,'LineStyle','--'); %Grey
    hold on
    g4a = plot(solPGT121TwoStrainMutation.x(:),max(solPGT121TwoStrainMutation.y(3,:),0),'LineWidth',1.75,'Color',[118,42,131]/256,'LineStyle','--'); %grey [118,42,131][26,152,80]/255
    hold on
    set(gca,'yscale','log')
    ylim([10^(-2),10^6])
    xlim([0,inf])
    ylabel('Viral Load (copies/mL) ','FontSize',15);  
    yyaxis right
    set(gca,'yscale','log')
    g31 = plot(solPGT121TwoStrainMutation.x(:),solPGT121TwoStrainMutation.y(4,:),'LineWidth',2,'Color',[227,26,28]./256,'LineStyle','-');
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

function Loss = InfectivityResistantLoss(Ab,t,PA)
if t < PA.tau
    Loss =1;
else
    Loss = 1./(1+PA.alphaR.*Ab);
end
end

function [sol] = PGT121_TwoStrainMutation_Solver(TotalTime,IC,PA) %ODE model without therapy
opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',1e-2);
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

end


