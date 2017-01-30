%save greedy lineaging result in esequence data structure in  acetree readable format.
function saveGreedyNucleiFiles(esequence,time,directory,anisotropy,ROIxmin,ROIymin)
set(0,'RecursionLimit',max(500,time*2));%if need more set recursion higher
if(~exist('ROIxmin','var'))
    ROIxmin=1;
    ROIymin=1;
end
nextnum=1;%nuclear 'names'

%initialize traversed data structure
traversed=cell(time,1);
%for t=1:time
%    traversed{t}=zeros(size(esequence{t}.finalpoints,1),1,'uint8');
%end
maxpoint=false;
for t=1:time
    if(~isempty(esequence{t}.finalpoints))
        %note off by 2 correction here ROI is first included pixel so
        %pos=pos+roi-1  Acetree coordinate system is 0 origin rather than 1
        %which is the origin of second subtraction
        esequence{t}.finalpoints(:,1)=  esequence{t}.finalpoints(:,1)+ROIxmin-2;
        esequence{t}.finalpoints(:,2)=  esequence{t}.finalpoints(:,2)+ROIymin-2;
    end
end

%traverse points outputting info to array
%and then outputing to file once traversed

currentcounters=ones(1,time);
assembledData=cell(1,time);
traversed=cell(1,time);
for i=1:length(traversed)
    traversed{i}=zeros(size(esequence{i}.finalpoints,1),1,'uint8');
end


for t=1:time %we stop at time -1 because any initiation at endtime time is an isolated point
    for i=1:size(esequence{t}.finalpoints,1)
    %    if(esequence{t}.suc(i,1)~=-1&~esequence{t}.delete(i)) %not an isolated unlinked point or slated for deletion
       if(~esequence{t}.delete(i)) %not an isolated unlinked point or slated for deletion
    
            if(~traversed{t}(i)) %not traversed so its sucessors need to be entered
                %output the points info
                %not a death, or a dead successor or endtime
                if(t<time&&esequence{t}.suc(i,1)~=-1&&~esequence{esequence{t}.suc_time(i,1)}.delete(esequence{t}.suc(i,1)))
                    data=[-1,currentcounters(t+1)]; %pred,suc1
                    if(esequence{t}.suc(i,2)~=-1&&~esequence{esequence{t}.suc_time(i,2)}.delete(esequence{t}.suc(i,2)))
                        data=[data,currentcounters(t+1)+1];%suc2
                    else
                        data=[data,-1];%suc2
                    end
                else
                    data=[-1,-1,-1]; %pred,suc1
                end
                
                
                data=[data,esequence{t}.finalpoints(i,:),esequence{t}.finaldiams(i),nextnum,esequence{t}.totalGFP(i)];%x,y,z,diam,nextnum(used in name)
                
                %assembledData{t}(currentcounters(t))=data;
                assembledData{t}=[assembledData{t};data];
                
                
                currentcounters(t)=currentcounters(t)+1;
                traversed{t}(i)=1;
                nextnum=nextnum+1;
                if(t<time&esequence{t}.suc(i,1)~=-1&~esequence{esequence{t}.suc_time(i,1)}.delete(esequence{t}.suc(i,1)))
                    %recurse daughter 1
             
                    [traversed,currentcounters,assembledData,nextnum]=recursive_traverse_output(esequence,traversed,currentcounters,assembledData,t,i,1,nextnum,time);
                    if(esequence{t}.suc(i,2)~=-1&~esequence{esequence{t}.suc_time(i,2)}.delete(esequence{t}.suc(i,2)))
                        [traversed,currentcounters,assembledData,nextnum]=recursive_traverse_output(esequence,traversed,currentcounters,assembledData,t,i,2,nextnum,time);
                    end
                end
            end
        end
    end
end




%ouput assembled data

for i=1:time%end_time %per time
    
    outputfile=[directory,'/t',num2str(i,'%03d'),'-nuclei'];
    fid=fopen(outputfile,'w');
    
    
    for p=1:size(assembledData{i},1)%per point
  %note arbitrary attempt to keep gfp total in grange acetree will read
            
        fprintf(fid,'%d,  %d,  %d,  %d,  %d,  %d, %d,  %2.1f,  %d,  %s,  %d, %s\n',...
            p,1,round(assembledData{i}(p,1:7)),['Nuc',num2str(assembledData{i}(p,8))],uint16(assembledData{i}(p,9)/256),' 0, 0, 0, , 0, 0, 0, 0, 0,');

    end
    
    
    fclose(fid);

end

end %end function