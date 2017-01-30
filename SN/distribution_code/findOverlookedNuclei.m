function nucleiSet=findOverlookedNuclei(X,nucleiSet,diskSet,anisotropy,celldiameter,numcells,parameters)

nucleiSet.centers=cell(length(nucleiSet.centerindicies),1);
nucleiSet.logodds=cell(length(nucleiSet.centerindicies),1);
nucleiSet.range=cell(length(nucleiSet.centerindicies),1);
for i=1:length(nucleiSet.centerindicies)
    nucleiSet.centers{i}=assign_planes(diskSet.centeredxymax(nucleiSet.centerindicies(i),:),diskSet.xydetdiameters(nucleiSet.centerindicies(i)),X,diskSet.centeredxymax,diskSet.xydetdiameters,anisotropy);
    nucleiSet.logodds{i}=calculateLogodds(nucleiSet.centers{i},diskSet.centeredxymax(nucleiSet.centerindicies(i),:),diskSet.xymaximavals(nucleiSet.centerindicies(i),:),diskSet.xydetdiameters(nucleiSet.centerindicies(i)),diskSet.xymaximavals,diskSet.xydetdiameters,anisotropy,diskSet.xycoverage);
    nucleiSet.range{i}=vcalculateMaximalRange(nucleiSet.centers{i},diskSet.centeredxymax(nucleiSet.centerindicies(i),:),nucleiSet.logodds{i},diskSet.xycoverage,length(nucleiSet.centers),parameters);
end

 %[centers2,logodds2,range2,unclaimedcentered,unclaimeddiams_out,unclaimedmaximavals_out]=iterate2ndround(X,oldindicies,centers2,logodds2,range2,unclaimed, unclaimedmaximavals,unclaimedcentered,unclaimeddiams,unclaimedcoverage,xymax,centeredxymax,xymaximavals,xydetdiameters,xycoverage,anisotropy,celldiameter,previous,numcells)
%indicies are unclaimed in previous 
 
%second round iterated
%unclaimed subset
[unclaimed,indicies]=removeClaimed(nucleiSet.centers,nucleiSet.range,diskSet.xymax);


%need to and indicies with previous indicies and keep subset
unclaimed=(diskSet.xymax(indicies,:));
unclaimedmaximavals=diskSet.xymaximavals(indicies,:);
unclaimedcentered3=diskSet.centeredxymax(indicies,:);
unclaimeddiams=diskSet.xydetdiameters(indicies);
unclaimedcoverage=diskSet.xycoverage(indicies);

%filtered to be isolated
if(~isempty(unclaimed))
[newfilteredpoints,indicies2]=closenessFilterRemovingIsolated(unclaimed,unclaimedmaximavals,.75*celldiameter,anisotropy,unclaimedcoverage);
else
    indicies2=[];
end
while (~isempty(find(indicies2))) 


%vals to go with filtered
unclaimedmaximavals=unclaimedmaximavals(indicies2);
unclaimedcentered3=unclaimedcentered3(indicies2,:);
unclaimeddiams=unclaimeddiams(indicies2);
unclaimedcoverage=unclaimedcoverage(indicies2);

%allocate next round
centers3=cell(length(unclaimeddiams),1);
logodds3=cell(length(unclaimeddiams),1);
range3=cell(length(unclaimeddiams),1);

for i=1:length(unclaimeddiams) 
 centers3{i}=assign_planes(unclaimedcentered3(i,:),unclaimeddiams(i),X,diskSet.centeredxymax,diskSet.xydetdiameters,anisotropy);
 logodds3{i}=calculateLogodds(centers3{i},unclaimedcentered3(i,:),unclaimedmaximavals(i),unclaimeddiams(i),diskSet.xymaximavals,diskSet.xydetdiameters,anisotropy,diskSet.xycoverage);
 range3{i}=vcalculateMaximalRange(centers3{i},unclaimedcentered3(i,:),logodds3{i},diskSet.xycoverage,numcells,parameters);
end


test=find(indicies);
nucleiSet.centerindicies=[nucleiSet.centerindicies;test(indicies2)'];

nucleiSet.centers={nucleiSet.centers{:},centers3{:}};
nucleiSet.logodds={nucleiSet.logodds{:},logodds3{:}};
nucleiSet.range={nucleiSet.range{:},range3{:}};

%unclaimedcentered=[unclaimedcentered;unclaimedcentered3];
%unclaimeddiams_out=[unclaimeddiams_out,unclaimeddiams];
%unclaimedmaximavals_out=[unclaimedmaximavals_out;unclaimedmaximavals];

%iteration
oldindicies=indicies;
[unclaimed,indicies]=removeClaimed(nucleiSet.centers,nucleiSet.range,diskSet.xymax);
indicies=indicies&oldindicies;
%need to and indicies with previous indicies and keep subset
unclaimed=(diskSet.xymax(indicies,:));
unclaimedmaximavals=diskSet.xymaximavals(indicies,:);
unclaimedcentered3=diskSet.centeredxymax(indicies,:);
unclaimeddiams=diskSet.xydetdiameters(indicies);
unclaimedcoverage=diskSet.xycoverage(indicies);

%filtered to be isolated
[newfilteredpoints,indicies2]=closenessFilterRemovingIsolated(unclaimed,unclaimedmaximavals,.75*celldiameter,anisotropy,unclaimedcoverage);

end