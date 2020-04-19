clear all
close all

[~,rawExcel,~] = xlsread('C:\Users\rbedder\Dropbox\MATLAB\PhD_Matlab\Academic_Travel_Survey\travelviewsData.xlsx');
[numExcel,~,~] = xlsread('C:\Users\rbedder\Dropbox\MATLAB\PhD_Matlab\Academic_Travel_Survey\travelviewsNumerical.xlsx');

%edits to excel

%1 means junior (research assistant, lab manager, transition to postdoc
%7 means senior (midcareer faculty,

%When they had to give number of flights etc, the 0 exported as numerical
%and 1-2 etc as text. Changed this row to text in the data extraction and
%changed 0 to 0-0

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

%% make an index for some FAQ

i.position    =   47;
i.gender      =   49;
i.location    =   50;
i.funder      =   53;
i.disability  =   56;

%% match the text to the numbers

cols = [3:59];

for groups = cols %...for each column of the excel (i.e. for each individual response)
    
tmpNum   =   unique(numExcel(:,groups),'stable'); 
tmpText  =   unique(textExcel(:,groups),'stable');   

if cellfun(@isempty,tmpText)

tmpText  =  cellstr(num2str(tmpNum));
    
end

tmpText(cellfun(@isempty, tmpText)) = [];

if sum(isnan(tmpNum)~=0)
% tmpText(isnan(tmpNum)) = []; 
tmpNum(isnan(tmpNum)) = [];     tmpNum(end+1) = nan;
tmpText(isnan(tmpNum)) = {'NaN'};
end %it seems that the first one will be empty, not the last? what is empty and not a nan
    
[tmpNum,sortNum] = sort(tmpNum);

respKey{groups} = [num2cell(tmpNum) tmpText(sortNum)];

end
    




%% new numerical analysis

% yQuest      =     [3 5 6 16 17 40]; %the corresponds to rows in the quesExcel (so individual questions)
yQuest      =     [31]; %the corresponds to rows in the quesExcel (so individual questions)

%build an index key
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


% %build an index key
% tableName           =       {'JuniorSenior'};
% i.junior  = 1;
% i.senior  = 2;
% groupUnique         =        [i.junior i.senior];
% groupLabels         =        {'Junior','Senior'};
% groupColors         =        {[0.75 0.75 0.75],[0.4 0.4 0.4]}
% 
% groups              =       numExcel(:,47);
% groups(groups<5)    =       i.junior;
% groups(groups>=5)   =       i.senior;

%build an index key
%tableName                     =       {'OriginalSplits'}
% groupUnique                  =       [4 5 6 7]; %...here refer to respKey{47} for the numbers in the first column
% groupLabels                  =       respKey{47}(groupUnique,2);
% groupUnique                  =       cell2mat(respKey{47}(groupUnique,1))';
% groupColors(groupUnique)     =       {[1 1 1],[0.75 0.75 0.75],[0.4 0.4 0.4],[0 0 0]};
% 
% groups                              =       numExcel(:,47);
% groups(~ismember(groups,groupUnique))  =       nan;




for Y = 1:length(yQuest); %...for each different question (is it's own figure)
    
    figure; f = 1;getTicks = [];
       
yResp               =     respKey{yQuest(Y)}(:,2); %...get the different responses for this question
yNums               =     cell2mat(respKey{yQuest(Y)}(:,1)); %...get the number keys for each answer      
yResp(isnan(yNums)) =     []; %...remove nans caused by people skipping the question
yNums(isnan(yNums)) =     []; %...remove nans caused by people skipping the question
yData               =     numExcel(:,yQuest(Y));%...get the response data for each participant


%%build a table of proportions and other details
tmpTable = array2table([num2cell(NaN(length(yResp)+10,length(groupLabels)));[{'_'},{'empty'}]]);
% tmpTable.Properties.VariableTypes = {repmat('double',length(yResp)+10,1);'string'};
tmpTable.Properties.VariableNames = groupLabels;
tmpTable.Properties.RowNames = [yResp;'Total';'Median';'Mean';'Std';'zStat';'zPvalue';'tStat';'tPvalue';...
                                        'zStat Compare';'zP Compare';'Question'];
tmpTable(length(yResp)+11,2) = quesExcel(yQuest(Y),2);                                  

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

title(quesExcel(yQuest(Y),2),'FontSize',20)
ylabel('Proportion of Responses')
set(gca,'xtick',[xLabels],'xticklabels',yResp,'FontSize',20)
xtickangle(45)
axis square
ylim([0 0.8])
xlim([0 max(xLabels)+ceil(length(groupUnique)/2)+1])
legend(groupLabels)

for sub = groupUnique
tmpTable{length(yResp)+1,sub}     =      {sum(cell2mat(tmpTable{1:length(yResp),sub}))};
end;questTables.(tableName{1}).(strcat('quest_',num2str(yQuest(Y)))) = tmpTable; clear tmpTable;

end




%% Is there support for changing the way we travel?
% Here we provide all the statistics given for this paragraph in a table

%take an hour to build a table, you'll thank yourself...to give all the
%percentages for different groups

JUNIOR      =       numExcel(:,47)<5 & numExcel(:,47)~=0;
SENIOR      =       numExcel(:,47)>=5 & numExcel(:,47)~=0;
ALL         =       numExcel(:,47)==0;

groupidx    =       [JUNIOR SENIOR ALL];  

questions   =       [3 40];

for group = 1:3;
    
idx = groupidx(:,group);
    
for q = 1:length(questions)


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





