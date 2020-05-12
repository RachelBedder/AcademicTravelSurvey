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

quests = [1:56]; %...this refers to the columns in the excel files

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
quests          =       [23:32]; 
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
groupLabels         =        {'Junior','Senior'}; %stats just between 1 and 2
groupColors         =        {[0.75 0.75 0.75],[0.4 0.4 0.4]}

groups              =       numExcel(:,47);
groups(groups<5)    =       i.junior;
groups(groups>=5)   =       i.senior;

elseif splitGroup == [3];

tableName                    =       {'OriginalSplits'}
groupUnique                  =       [4 5 6 7]; %...here refer to respKey{47} for the numbers in the first column
groupLabels                  =       respKey{47}(groupUnique,2);
groupUnique                  =       cell2mat(respKey{47}(groupUnique,1))';
groupColors(groupUnique)     =       {[1 1 1],[0.75 0.75 0.75],[0.4 0.4 0.4],[0 0 0]};

groups                              =       numExcel(:,47);
groups(~ismember(groups,groupUnique))  =       nan;

elseif splitGroup == [4];

tableName           =       {'DisabilityChronic'};
i.dischron  = 1;
i.non  = 2;
groupUnique         =        [i.dischron i.non];
groupLabels         =        {'Dishchron','Non'}; %stats just between 1 and 2
groupColors         =        {[0.75 0.75 0.75],[0.4 0.4 0.4]};

groups              =       numExcel(:,55);
groups(groups==1)    =       i.dischron;
groups(groups~=1)   =       i.non;

elseif splitGroup == [5];

tableName           =       {'Gender'};
i.male  = 1;
i.female  = 2;
groupUnique         =        [i.male i.female];
groupLabels         =        {'Male','Female'}; %stats just between 1 and 2
groupColors         =        {[0.75 0.75 0.75],[0.4 0.4 0.4]};

groups              =       numExcel(:,49);
groups(groups==1)    =       i.male;
groups(groups~=1)   =       i.female;


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
    [p,~,stats]   =   ranksum(yData(groups==1 & ~isnan(yData)),yData(groups==2 & ~isnan(yData)));
   try %doesn't return a zvalue if the sample is too small (I think)
    tmpTable{length(yResp)+9,1}  =   {stats.zval};  
   catch
       tmpTable{length(yResp)+9,1} = {NaN};
   end
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


%% Look at the disability questions seperately

disIdx          =   numExcel(:,55)==1;

quests          =       [16:22]; 

for iQ = quests;
countResp(:,1)      =   num2cell(unique(numExcel(disIdx,iQ)));
countResp(:,2)      =   arrayfun(@(x)length(find(numExcel(disIdx,iQ) == x)), unique(numExcel(disIdx,iQ)), 'Uniform', false);
countResp(:,[3 4])  =   respKey{iQ};

quesExcel{iQ,:}
countResp

clear countResp
end



%% In lieu of the appropriate statistical test, look plot average differences between groups
clear groupidx
groupidx(:,1)              =       numExcel(:,47);      groupPairs{1} = {'Junior','Senior'};
groupidx(groupidx(:,1)<5) = 1;groupidx(groupidx(:,1)>=5) = 2; %SAVE INDEXES MADE ABOVE INSTEAD
groupidx(:,2)              =       numExcel(:,49);      groupPairs{2} = {'Male','Female'};

groupidx(groupidx>2) = NaN;

groupidx(:,3) = nan(length(groupidx),1);
groupidx(groupidx(:,1)==1 & groupidx(:,2)==1,3) = [1]; %junior and male
groupidx(groupidx(:,1)==1 & groupidx(:,2)==2,3) = [2]; %junior and female
groupidx(groupidx(:,1)==2 & groupidx(:,2)==1,3) = [3]; %senior and male
groupidx(groupidx(:,1)==2 & groupidx(:,2)==2,3) = [4]; %senior and female;  
groupPairs{3} = {'JunMal','JunFem','SenMal','SenFem'};

groupidx(:,4)              =       numExcel(:,53);      groupPairs{4} = {'University','Grant'};
groupidx(groupidx(:,4)==1,4) = 1;  groupidx(groupidx(:,4)~=1,4) = 2; %we assume other is a grant here, be more specific later

questPairs = {[10 14],[8 12],[10 11],[8 9]};
questTitles = {'Should support carbon offsetting','Support Travel by train were possible','Universities support of carbon offsetting','Universities support train travel'}
questLabels = {{'Universitory','Funders'},{'Universitory','Funders'},{'Voluntary','Mandatory'},{'Voluntary','Mandatory'}}

for G = 1:size(groupidx,2)
figure;
for iQ = 1:size(questPairs,2)
quests = questPairs{iQ};
groups    =   groupidx(:,G);

%FIX ERRORBARS FOR THE NANS
for j = unique(groups)'
subplot(2,2,iQ)
errorbar([1 2]+j*0.1,nanmedian(numExcel(groups==j,quests)),nanstd(numExcel(groups==j,quests))./sqrt(sum(groups==j & ~isnan(numExcel(:,quests)))),'LineWidth',3); hold on
end
set(gca,'xtick',[1 2],'xticklabels',questLabels{iQ})
title(questTitles{iQ})
xlim([0 3]);
ylim([-2 4])
if iQ == 1;
    legend(groupPairs{G})
end

end
end

%% Plotting the incentives

quests = [33:36];

%plot average ranking for each incentive for each group

groups = groupidx(:,1);

for G = 1:size(groupidx,2)
figure;

groups    =   groupidx(:,G);

for j = unique(groups)'
errorbar(nanmean(numExcel(groups==j,quests)),nanstd(numExcel(groups==j,quests))./sqrt(sum(groups==j & ~isnan(numExcel(:,quests)))),'LineWidth',3); hold on
end

xlim([0 5])
set(gca, 'YDir','reverse')
legend(groupPairs{G})
set(gca,'xtick',[1:4],'xticklabels',{'First Class','Overnight Stay','Extra Annual Leave','Green in Grants'});xtickangle(45)
ylabel('Average Ordered Ranking (1 is most prefered)')
end


%% Ranking concerns on increased train by group

quests = [23:31];

groups = groupidx(:,3);

for G = 1:size(groupidx,2)
figure;

groups    =   groupidx(:,G);

for j = unique(groups)'
errorbar(nanmean(numExcel(groups==j,quests)),nanstd(numExcel(groups==j,quests))./sqrt(sum(groups==j & ~isnan(numExcel(:,quests)))),'LineWidth',3); hold on
end

xlim([0 10]);ylim([0 2])
legend(groupPairs{G})
set(gca,'ytick',[0 1 2],'yticklabels',{'Not a concern','Minor Concern','Major Concern'})
set(gca,'xtick',[1:9],'xticklabels',{'Too expensive','Too long','Time away from home','Time away from science'...
                            'Lack of support from seniors','Take AL','Seen less serious','Harm Junior','Harm Senior'});xtickangle(45)
% ylabel('Average Ordered Ranking (1 is most prefered)')
end

%% Ranking concerns on reduced flying by group

quests = [16:22];

%plot average ranking for each incentive for each group

groups = groupidx(:,1);

for G = 1:size(groupidx,2)
figure;

groups    =   groupidx(:,G);

for j = unique(groups)'
errorbar(nanmean(numExcel(groups==j,quests)),nanstd(numExcel(groups==j,quests))./sqrt(sum(groups==j & ~isnan(numExcel(:,quests)))),'LineWidth',3); hold on
end

xlim([0 8]);ylim([0 2])
legend(groupPairs{G})
set(gca,'ytick',[0 1 2],'yticklabels',{'Not a concern','Minor Concern','Major Concern'})
set(gca,'xtick',[1:7],'xticklabels',{'Giving a video talk','Attending video talk','Reduced network','Lack of support from senior'...
                            'Seen less serious','Harm Junior','Harm Senior'});xtickangle(45)
% ylabel('Average Ordered Ranking (1 is most prefered)')
end







