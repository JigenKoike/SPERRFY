% sandbox_for_randomtest_visualization

% memo: Load data before runnning


%% 

%% Compare　Correlation

DimShow = 5;
r_data = mouseConnectomeDataAnalysis.OverallModelAndResults.HoldoutSummary.TestAllCorrelations(:,1:DimShow);
r_null1 = nullBatchRunner_1.HoldoutTestCorrMeans(:,1:DimShow);
r_null2 = nullBatchRunner_2.HoldoutTestCorrMeans(:,1:DimShow);
r_null3 = nullBatchRunner_3.HoldoutTestCorrMeans(:,1:DimShow);
r_null4 = nullBatchRunner_4.HoldoutTestCorrMeans(:,1:DimShow);
r_null5 = nullBatchRunner_5.HoldoutTestCorrMeans(:,1:DimShow);
r_null6 = nullBatchRunner_6.HoldoutTestCorrMeans(:,1:DimShow);

%%
intervalWidth = 5;
x_data = ones(height(r_data),1) * (1:intervalWidth:1+intervalWidth*(DimShow-1));
x_gen = ones(height(r_null1),1) * (2:intervalWidth:2+intervalWidth*(DimShow-1));
x_dis = ones(height(r_null1),1) * (3:intervalWidth:3+intervalWidth*(DimShow-1));
x_net = ones(height(r_null1),1) * (4:intervalWidth:4+intervalWidth*(DimShow-1));


%%
legendsStr = ["Original Connectome","Random Generation","Distance Preserved","Network Preserved"]; 
xticksVec = [intervalWidth/2:intervalWidth:2.5+intervalWidth*(DimShow-1)];
xticklabelsVec = [1:5];
xlimVec = [0,intervalWidth*DimShow];
ylimVec = [-0.1:1.01];




%%
sz = 3;
szData = 6;
colorListTmp = colorList13([1,3,7,4],:);

%% Global 

figure;
h1 = swarmchart(x_data,r_data,szData,colorListTmp(1,:),'filled');
hold on
h2 = swarmchart(x_gen,r_null1,sz,colorListTmp(2,:),'filled');
hold on
h3 = swarmchart(x_dis,r_null3,sz,colorListTmp(3,:),'filled');
hold on
h4 = swarmchart(x_net,r_null5,sz,colorListTmp(4,:),'filled');
legend([h1(1),h2(1),h3(1),h4(1)],legendsStr);
xlim(xlimVec)
xlabel("Correlation Component (Rank)")
xticks(xticksVec)
xticklabels(xticklabelsVec)
ylim(ylimVec)
ylabel("Correlation Coefficient")
title("Globally Randomized Model")


%% Local
figure;
h1 = swarmchart(x_data,r_data,szData,colorListTmp(1,:),'filled');
hold on
h2 = swarmchart(x_gen,r_null2,sz,colorListTmp(2,:),'filled');
hold on
h3 = swarmchart(x_dis,r_null4,sz,colorListTmp(3,:),'filled');
hold on
h4 = swarmchart(x_net,r_null6,sz,colorListTmp(4,:),'filled');
legend([h1(1),h2(1),h3(1),h4(1)],legendsStr);
xlim(xlimVec)
xlabel("Correlation Component (Rank)")
xticks(xticksVec)
xticklabels(xticklabelsVec)
ylim(ylimVec)
ylabel("Correlation Coefficient")
title("Locally Randomized Model")


%% Global & Local
%%
intervalWidth = 8;
x_data = ones(height(r_data),1) * (1:intervalWidth:1+intervalWidth*(DimShow-1));
x_genL = ones(height(r_null1),1) * (2:intervalWidth:2+intervalWidth*(DimShow-1));
x_disL = ones(height(r_null1),1) * (3:intervalWidth:3+intervalWidth*(DimShow-1));
x_netL = ones(height(r_null1),1) * (4:intervalWidth:4+intervalWidth*(DimShow-1));
x_genG = ones(height(r_null1),1) * (5:intervalWidth:5+intervalWidth*(DimShow-1));
x_disG = ones(height(r_null1),1) * (6:intervalWidth:6+intervalWidth*(DimShow-1));
x_netG = ones(height(r_null1),1) * (7:intervalWidth:7+intervalWidth*(DimShow-1));


%%
legendsStr = ["Original Connectome","Random Generation","Distance Preserved","Network Preserved"]; 
xticksVec = [intervalWidth/2:intervalWidth:intervalWidth/2+intervalWidth*(DimShow-1)];
xticklabelsVec = [1:5];
xlimVec = [0,intervalWidth*DimShow];
ylimVec = [-0.15,1.01];

%%
figure;
h1 = swarmchart(x_data,r_data,szData,colorListTmp(1,:),'filled');
hold on
h2 = swarmchart(x_genL,r_null2,sz,colorListTmp(2,:),'filled');
hold on
h3 = swarmchart(x_disL,r_null4,sz,colorListTmp(3,:),'filled');
hold on
h4 = swarmchart(x_netL,r_null6,sz,colorListTmp(4,:),'filled');
hold on
h5 = swarmchart(x_genG,r_null1,sz,colorListTmp(2,:),'filled');
hold on
h6 = swarmchart(x_disG,r_null3,sz,colorListTmp(3,:),'filled');
hold on
h7 = swarmchart(x_netG,r_null5,sz,colorListTmp(4,:),'filled');
legend([h1(1),h2(1),h3(1),h4(1)],legendsStr);
xlim(xlimVec)
xlabel("Correlation Component (Rank)")
xticks(xticksVec)
xticklabels(xticklabelsVec)
ylim(ylimVec)
ylabel("Correlation Coefficient")
title("Globally Randomized Model")

%% AUC
AUC_data = mouseConnectomeDataAnalysis.OverallModelAndResults.HoldoutSummary.TestAllAUCs_ROC.';
AUC_null1 = nullBatchRunner_1.HoldoutTestAUCMeans;
AUC_null2 = nullBatchRunner_2.HoldoutTestAUCMeans;
AUC_null3 = nullBatchRunner_3.HoldoutTestAUCMeans;
AUC_null4 = nullBatchRunner_4.HoldoutTestAUCMeans;
AUC_null5 = nullBatchRunner_5.HoldoutTestAUCMeans;
AUC_null6 = nullBatchRunner_6.HoldoutTestAUCMeans;

%%
x_data = ones(height(AUC_data),1) * 1;
x_gen = ones(height(AUC_null1),1) * 2;
x_dis = ones(height(AUC_null1),1) * 3;
x_net = ones(height(AUC_null1),1) * 4;

%%
legendsStr = ["Original Connectome","Random Generation","Distance Preserved","Network Preserved"];
intervalWidth = 5; 
xlimVec = [0,5];
ylimVec = [0.4,1];


%% Global 

figure;
h1 = swarmchart(x_data,AUC_data,szData,colorListTmp(1,:),'filled');
hold on
h2 = swarmchart(x_gen,AUC_null1,sz,colorListTmp(2,:),'filled');
hold on
h3 = swarmchart(x_dis,AUC_null3,sz,colorListTmp(3,:),'filled');
hold on
h4 = swarmchart(x_net,AUC_null5,sz,colorListTmp(4,:),'filled');
legend([h1,h2,h3,h4],legendsStr);
xlim(xlimVec)
xlabel("Correlation Component (Rank)")
ylim(ylimVec)
yticks([0.5:0.1:1])
ylabel("AUC")
title("Globally Randomized Model")


%% Local 
figure;
h1 = swarmchart(x_data,AUC_data,szData,colorListTmp(1,:),'filled');
hold on
h2 = swarmchart(x_gen,AUC_null2,sz,colorListTmp(2,:),'filled');
hold on
h3 = swarmchart(x_dis,AUC_null4,sz,colorListTmp(3,:),'filled');
hold on
h4 = swarmchart(x_net,AUC_null6,sz,colorListTmp(4,:),'filled');
legend([h1,h2,h3,h4],legendsStr);
xlim(xlimVec)
xlabel("Correlation Component (Rank)")
ylim(ylimVec)
yticks([0.5:0.1:1])
ylabel("AUC")
title("Locally Randomized Model")

%% PI Similarity

DimShow = 5;

PISim_null1 = nullBatchRunner_1.WirigPISimilarities(:,1:DimShow);
PISim_null2 = nullBatchRunner_2.WirigPISimilarities(:,1:DimShow);
PISim_null3 = nullBatchRunner_3.WirigPISimilarities(:,1:DimShow);
PISim_null4 = nullBatchRunner_4.WirigPISimilarities(:,1:DimShow);
PISim_null5 = nullBatchRunner_5.WirigPISimilarities(:,1:DimShow);
PISim_null6 = nullBatchRunner_6.WirigPISimilarities(:,1:DimShow);


%%
intervalWidth = 4;
x_gen = ones(height(PISim_null1),1) * (1:intervalWidth:2+intervalWidth*(DimShow-1));
x_dis = ones(height(PISim_null1),1) * (2:intervalWidth:3+intervalWidth*(DimShow-1));
x_net = ones(height(PISim_null1),1) * (3:intervalWidth:4+intervalWidth*(DimShow-1));
xticklabelsVec = [1:5];
xticksVec = [intervalWidth/2:intervalWidth:2.5+intervalWidth*(DimShow-1)];
xlimVec = [0,intervalWidth*DimShow];
ylimVec = [0,1];

%% Global
figure;
h2 = swarmchart(x_gen,PISim_null1,sz,colorListTmp(2,:),'filled');
hold on
h3 = swarmchart(x_dis,PISim_null3,sz,colorListTmp(3,:),'filled');
hold on
h4 = swarmchart(x_net,PISim_null5,sz,colorListTmp(4,:),'filled');
legend([h2(1),h3(1),h4(1)],legendsStr(2:4));
xlim(xlimVec)
xlabel("Correlation Component (Rank)")
xticks(xticksVec)
xticklabels(xticklabelsVec)
ylim(ylimVec)
ylabel("PI Similarity Score")
title("Globally Randomized Model")

%% Local
figure;
h2 = swarmchart(x_gen,PISim_null2,sz,colorListTmp(2,:),'filled');
hold on
h3 = swarmchart(x_dis,PISim_null4,sz,colorListTmp(3,:),'filled');
hold on
h4 = swarmchart(x_net,PISim_null6,sz,colorListTmp(4,:),'filled');
legend([h2(1),h3(1),h4(1)],legendsStr(2:4));
xlim(xlimVec)
xlabel("Correlation Component (Rank)")
xticks(xticksVec)
xticklabels(xticklabelsVec)
ylim(ylimVec)
ylabel("PI Similarity Score")
title("Locally Randomized Model")