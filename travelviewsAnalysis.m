clear all
close all

%% Academic Travel Survey - Online November 2019
% Survey by Rachel Bedder and Steve Fleming (UCL)
% Code by Rachel Bedder (rachel.bedder.15@ucl.ac.uk)

% Updates posted on https://github.com/RachelBedder/AcademicTravelSurvey


%% ...extract data from excel files
cd
[~,rawExcel,~] = xlsread('travelviewsData.xlsx');
[numExcel,~,~] = xlsread('travelviewsNumerical.xlsx');



%% get the excel labels

rawExcel(:,[1:9]) = [];
numExcel(:,[1:9]) = []; %numExcel skips the questions

quesExcel   =   rawExcel(1:2,:)';
textExcel   =   rawExcel(3:end,:);

%% remove all incomplete data, or those that did not consent for their data to be used

consent     =   cellfun('isempty',textExcel(:,2));
textExcel   =   textExcel(consent,:);
numExcel    =   numExcel(consent,:);

incomplete  =   sum(cellfun('isempty',textExcel),2)>50;
textExcel   =   textExcel(~incomplete,:);
numExcel    =   numExcel(~incomplete,:);


%% For each of the questions, match the text answers to the numerical response

quests = [3:59]; %...this refers to the columns in the excel files

for iQ = quests 
    
tmpNum   =   unique(numExcel(:,iQ),'stable'); 
tmpText  =   unique(textExcel(:,iQ),'stable');   

if cellfun(@isempty,tmpText)

tmpText  =  cellstr(num2str(tmpNum));
    
end

tmpText(cellfun(@isempty, tmpText)) = [];

if sum(isnan(tmpNum)~=0)
tmpNum(isnan(tmpNum)) = [];     
tmpNum(end+1) = nan;
tmpText(isnan(tmpNum)) = {'NaN'};
end 
    
[tmpNum,sortNum] = sort(tmpNum);

respKey{iQ} = [num2cell(tmpNum) tmpText(sortNum)];

end
    




%% Create plots and tables for different splits of the data.

%... first identify the questions (columns in the excel data)
quests          =       [3 5 6 16 17 31 40]; 
%...chose an index key from below, this will just run one split at a time
%so do not add a vector. These will select from the index keys below
splitGroup      =       [2];

if splitGroup == [1];
tableName           =       {'CarbonNeutralAgree'};
i.lower  = 1;
i.higher  = 2;
i.all = 3
groupUnique         =        [i.lower i.higher];
groupLabels         =        {'Lower','Higher'};
groupColors         =        {[0.75 0.75 0.75],[0.4 0.4 0.4],[1 1 1]};

groups              =       numExcel(:,3);
groups(groups<=1)    =       i.lower;
groups(groups>1)   =       i.higher;

elseif splitGroup == [2];
    
tableName           =       {'JuniorSenior'};
i.junior  = 1;
i.senior  = 2;
groupUnique         =        [i.junior i.senior];
groupLabels         =        {'Junior','Senior'};
groupColors         =        {[0.75 0.75 0.75],[0.4 0.4 0.4]}

groups              =       numExcel(:,47);
groups(groups<5)    =       i.junior;
groups(groups>=5)   =       i.senior;

elseif splitGroup == [3];

tableName                     =       {'OriginalSplits'}
groupUnique                  =       [4 5 6 7]; %...here refer to respKey{47} for the numbers in the first column
groupLabels                  =       respKey{47}(groupUnique,2);
groupUnique                  =       cell2mat(respKey{47}(groupUnique,1))';
groupColors(groupUnique)     =       {[1 1 1],[0.75 0.75 0.75],[0.4 0.4 0.4],[0 0 0]};

groups                              =       numExcel(:,47);
groups(~ismember(groups,groupUnique))  =       nan;

end


for iQ = quests; %...for each different question (is it's own figure)
    
    figure; f = 1;getTicks = [];
       
yResp               =     respKey{iQ}(:,2); %...get the different responses for this question
yNums               =     cell2mat(respKey{iQ}(:,1)); %...get the number keys for each answer      
yResp(isnan(yNums)) =     []; %...remove nans caused by people skipping the question
yNums(isnan(yNums)) =     []; %...remove nans caused by people skipping the question
yData               =     numExcel(:,iQ);%...get the response data for each participant


%%build a table of proportions and other details
tmpTable = array2table([num2cell(NaN(length(yResp)+10,length(groupLabels)));[{'_'},{'empty'}]]);
tmpTable.Properties.VariableNames = groupLabels;
tmpTable.Properties.RowNames = [yResp;'Total';'Median';'Mean';'Std';'zStat';'zPvalue';'tStat';'tPvalue';...
                                        'zStat Compare';'zP Compare';'Question'];
tmpTable(length(yResp)+11,2) = quesExcel(iQ,2);                                  

for cat = 1:length(yResp) %...for each category of answer
   
for sub = groupUnique %...for each of the participant groups
       
groupidx  =     (groups==sub);
  
plotCat   =     sum(groupidx & yData==yNums(cat))/sum(groupidx); %...get the data for that response

bar(f,plotCat,'FaceColor',groupColors{sub},'EdgeColor','k'); hold on

f = f+1;
if sub == groupUnique(end)   %everytime we get to a new group, skip an extra xtick     
f = f+2;    
end

%get statistics!
tmpTable{cat,sub}                  =      {round(plotCat*100)/100};
tmpTable{length(yResp)+2,sub}      =      {nanmedian(yData(groupidx))};
tmpTable{length(yResp)+3,sub}      =      {nanmean(yData(groupidx))};
tmpTable{length(yResp)+4,sub}      =      {nanstd(yData(groupidx))};

[pZ,~,stats]           =    signrank(yData(groupidx));

try
    tmpTable{length(yResp)+5,sub}      =      {stats.zval};
    catch
    tmpTable{length(yResp)+5,sub}      =       {nan};
end
tmpTable{length(yResp)+6,sub}      =      {pZ};

[h,pT,~,stats]        =    ttest(yData(groupidx));

tmpTable{length(yResp)+7,sub}      =      {stats.tstat};
tmpTable{length(yResp)+8,sub}      =      {pT};

end

% if juniot and senior do some statistics here
if strcmp(tableName,'JuniorSenior')
%     
%     if Y==4;
%         abc = 1;
%     end
    
    %compare the average score 
    [p,~,stats]   =   ranksum(yData(groups==1),yData(groups==2));
    tmpTable{length(yResp)+9,1}  =   {stats.zval};  
    tmpTable{length(yResp)+10,1}  =  {p};     
    
end

end

%This builds where the question responses are on the x axis
xLabels = 1+ceil(length(groupUnique)/2):length(groupUnique)+2:(length(groupUnique)+2)*length(yResp);

title(quesExcel(iQ,2),'FontSize',20)
ylabel('Proportion of Responses')
set(gca,'xtick',[xLabels],'xticklabels',yResp,'FontSize',20)
xtickangle(45)
axis square
ylim([0 0.8])
xlim([0 max(xLabels)+ceil(length(groupUnique)/2)+1])
legend(groupLabels)

for sub = groupUnique
tmpTable{length(yResp)+1,sub}     =      {sum(cell2mat(tmpTable{1:length(yResp),sub}))};
end;questTables.(tableName{1}).(strcat('quest_',num2str(iQ))) = tmpTable; clear tmpTable;

end




%% Is there support for changing the way we travel?
% Here we provide all the statistics given for this paragraph in a table

%take an hour to build a table, you'll thank yourself...to give all the
%percentages for different groups

JUNIOR      =       numExcel(:,47)<5 & numExcel(:,47)~=0;
SENIOR      =       numExcel(:,47)>=5 & numExcel(:,47)~=0;
ALL         =       numExcel(:,47)==0;

groupidx    =       [JUNIOR SENIOR ALL];  

quests   =       [3 40];

for group = 1:3;
    
idx = groupidx(:,group);
    
for q = 1:length(quests)


noSubj(q,group)              =   sum(idx);

end
end
%%
% headers = {quesExcel{3 40}) 
% table(headers,



% What is the difference in their flying habits?
nanmedian(numExcel(SENIOR,40))
nanmedian(numExcel(JUNIOR,40))
[p h stats] = ranksum(numExcel(JUNIOR,3),numExcel(SENIOR,3))

% Academics should aim to be carbon neutral - [3] 

%...everyone
nanmedian(numExcel(ALL,3))
[p h stats] = signrank(numExcel(ALL,3))

%...junior academics
nanmedian(numExcel(JUNIOR,3))
[p h stats] = signrank(numExcel(JUNIOR,3))

%...senior academics
nanmedian(numExcel(SENIOR,3))
[p h stats] = signrank(numExcel(SENIOR,3))

%...junior vs senior
[p h stats] = ranksum(numExcel(JUNIOR,3),numExcel(SENIOR,3))





