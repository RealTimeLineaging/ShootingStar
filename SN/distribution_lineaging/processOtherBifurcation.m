function [esequence,count] = processOtherBifurcation( esequence,t,i,splitscores_div,trackingparameters,count,parameters)
%given that current bifurcation is 'other' ie. not a division but not a FP
%or FN as far as we can tell, process it by breaking weaker link
global FNtype;
global classround;

scores=nondivScoreModelCostFunction(esequence,i,t,...
    esequence{t}.suc(i,:)',esequence{t}.suc_time(i,:)',trackingparameters);


%if(splitscores_div(2)>splitscores_div(3)) %d1 is bigger (worse) than d2 match
%d1 is worse scoring against model
if(scores(1)>scores(2))
    otheri=esequence{t}.suc(i,1);
    esequence{esequence{t}.suc_time(i,1)}.pred(esequence{t}.suc(i,1))=-1;
    esequence{esequence{t}.suc_time(i,1)}.pred_time(esequence{t}.suc(i,1))=-1;
    esequence{t}.suc_time(i,1)=esequence{t}.suc_time(i,2);
    esequence{t}.suc(i,1)=esequence{t}.suc(i,2);
else %d2 is worse match
    otheri=esequence{t}.suc(i,2);
    esequence{esequence{t}.suc_time(i,2)}.pred(esequence{t}.suc(i,2))=-1;
    esequence{esequence{t}.suc_time(i,2)}.pred_time(esequence{t}.suc(i,2))=-1;
end
esequence{t}.suc_time(i,2)=-1;
esequence{t}.suc(i,2)=-1;

%now having broken worse off try attaching it to 4 back nn and iteratively
%testing
dis=distance_anisotropic(esequence{t}.finalpoints',esequence{t+1}.finalpoints(otheri,:)',trackingparameters.anisotropyvector);
solved=false;
cind=1;
%for nnind=1:4
while (~solved)

    [d,mind]=min(dis);
    if(mind~=i&&esequence{t}.suc(mind,2)==-1)% if not where it was attached previously
        %attach it
        if(esequence{t}.suc(mind,1)==-1)
            %if is a loose end just join it? this had little effect
            esequence{t}.suc(mind,1)=otheri;
            esequence{t}.suc_time(mind,1)=t+1;
            esequence{t+1}.pred(otheri)=mind;
            esequence{t+1}.pred_time(otheri)=t;
            solved=true;
        else
           
            %make bifurcation elsewhere
            esequence{t}.suc(mind,2)=otheri;
            esequence{t}.suc_time(mind,2)=t+1;
            esequence{t+1}.pred(otheri)=mind;
            esequence{t+1}.pred_time(otheri)=t;
            
            [d1cand,d2cand,d1candt,d2candt,d1length,d2length,...
                bestFNBackCorrect,bestmatchings,FNbackcand1lengths,FNbackcand2lengths,...
                bestFNForwardLengthD1,bestFNForwardLengthD2,....
                bestmatchingsplayerstart,bestmatchingplayerend,bestdaughter,bestIndex,...
                splitscores_div,alldaughterdata,allforwarddata,allbackdata,count] ....
                = assembleBifurcationData( esequence,t,mind,trackingparameters,count,parameters );

            %pick whether to call the single, or multiple statistical model
            %classifier based on whether it is a single or multi model
            %passed in by parameter file
            if (isfield(trackingparameters.bifurcationclassifier,'ambigious'))
            predicted_class = predictBifurcationType(...
                alldaughterdata,allforwarddata,allbackdata,d1length,d2length,...
                FNbackcand1lengths,FNbackcand2lengths,bestFNForwardLengthD1,...
                bestFNForwardLengthD2,bestFNBackCorrect,trackingparameters,bestIndex, count); 
           else
                   predicted_class = predictBifurcationTypeSinglemodel(...
                alldaughterdata,allforwarddata,allbackdata,d1length,d2length,...
                FNbackcand1lengths,FNbackcand2lengths,bestFNForwardLengthD1,...
                bestFNForwardLengthD2,bestFNBackCorrect,trackingparameters,bestIndex, count); 
           end
               minsize=min(d1length,d2length);
            classround=[classround;2];
            %process bifurcation based on predicted class
            if(predicted_class==0)
                %other
                FNtype=[FNtype;-1, -1];
                %this new location failed also break and retry in loop
                %above
                %test d1 or d2 worse
                scores=nondivScoreModelCostFunction(esequence,mind,t,...
                    esequence{t}.suc(mind,:)',esequence{t}.suc_time(mind,:)',trackingparameters);
                
                %if d1 worse restart from scratch 
                if (scores(2)<scores(1))
                    [esequence,count] = processOtherBifurcation( esequence,t,mind,splitscores_div,trackingparameters,count,parameters);
                    solved=true;%even if it failed we're done
                else
                    
                    %otherwise continue trying to move this loose end
                    esequence{t}.suc_time(mind,2)=-1;
                    esequence{t}.suc(mind,2)=-1;
                    esequence{t+1}.pred(otheri)=-1;
                    esequence{t+1}.pred_time(otheri)=-1;
                end
                %{
            if(predicted_class==0)
                %other
                FNtype=[FNtype;-1, -1];
                %this new location failed also break and retry in loop
                %above
                %test d1 or d2 worse
%                scores=nondivScoreModelCostFunction(esequence,mind,t,...
%    esequence{t}.suc(mind,:)',esequence{t}.suc_time(mind,:)',trackingparameters);

                %if d1 worse restart from scratch
                %if (scores(2)<scores(1))
                %esequence = processOtherBifurcation( esequence,t,mind,splitscores_div,trackingparameters );
            %else
            %end
                %otherwise continue trying
                esequence{t}.suc_time(mind,2)=-1;
                esequence{t}.suc(mind,2)=-1;
                esequence{t+1}.pred(otheri)=-1;
                esequence{t+1}.pred_time(otheri)=-1;
                %}
            else if(predicted_class==3)
                    %FP
                    FNtype=[FNtype;-1, -1];
                    esequence = processFPBifurcation( esequence,t,mind,minsize,d1length );
                    solved=true;
                    ' find solution fp'
                else if (predicted_class==2)
                        %FN
                        [esequence,numplayers,endplayerssameasdiv ]= processFNBifurcation...
                            ( esequence,t,mind,bestmatchings, bestmatchingsplayerstart,...
                            bestmatchingplayerend, bestdaughter,bestIndex,d1cand,d2cand,...
                            d1candt,d2candt,trackingparameters);
                        %  matchinfo
                        FNtype=[FNtype;numplayers, endplayerssameasdiv];
                        solved=true;
                        ' find solution FN'
                    else
                        %division
                        FNtype=[FNtype;-1, -1];
                        solved=true;
                        ' find solution div'
                    end %dont need to do anthing for class=1 division
                end
            end
            
        end
    end
    dis(mind)=inf; %remove rejected of nn from consideration
    cind=cind+1;
    %give up if dont get sucessfull within 4nn
    if(cind>4)
        solved=true;
        'failed to find solution'
    end
end



end

%{
%if time is past where we will likely edit leave bad divisons to avoid
%bogging down acetree which has a hard time with loose ends
%if (t<trackingparameters.edit_endtime)

scores=nondivScoreModelCostFunction(esequence,i,t,...
    esequence{t}.suc(i,:)',esequence{t}.suc_time(i,:)',trackingparameters);


%if(splitscores_div(2)>splitscores_div(3)) %d1 is bigger (worse) than d2 match
%d1 is worse scoring against model
if(scores(1)>scores(2))
esequence{esequence{t}.suc_time(i,1)}.pred(esequence{t}.suc(i,1))=-1;
    esequence{esequence{t}.suc_time(i,1)}.pred_time(esequence{t}.suc(i,1))=-1;
    esequence{t}.suc_time(i,1)=esequence{t}.suc_time(i,2);
    esequence{t}.suc(i,1)=esequence{t}.suc(i,2);
else %d2 is worse match
    
    esequence{esequence{t}.suc_time(i,2)}.pred(esequence{t}.suc(i,2))=-1;
    esequence{esequence{t}.suc_time(i,2)}.pred_time(esequence{t}.suc(i,2))=-1;
end
esequence{t}.suc_time(i,2)=-1;
esequence{t}.suc(i,2)=-1;
end
%end

%}