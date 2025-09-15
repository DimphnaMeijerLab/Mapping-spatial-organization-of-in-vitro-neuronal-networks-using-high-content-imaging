function [Opt_LBW, Data, threshold]=brightsizefilteralt(LBW_Neurons,Dapi)
    %% 
    Opt_LBW=LBW_Neurons;
    CC=bwconncomp(Opt_LBW);
    stats=regionprops('table',CC,Dapi,'Area','MeanIntensity','PixelValues');
    sigma=zeros(length(stats.PixelValues),1);
    for i=1:length(stats.PixelValues)
        sigma(i)=std(double(stats.PixelValues{i}));
    end
    %%
    [idx,C] = kmedoids([stats.Area sigma],2,'Distance','seuclidean','Replicates',5);
    if C(1,2)>C(2,2)
        Idx_intensity=(idx==1);
    else
        Idx_intensity=(idx==2);
    end
    %%
    Opt_LBW(cat(1,CC.PixelIdxList{Idx_intensity}))=false;
    %figure(),plot(stats.Area(idx==1),sigma(idx==1),'r.',stats.Area(idx==2),sigma(idx==2),'b.');
    
    %figure()
    %histogram(stats.Area,50);
    [N,X]=hist(stats.Area,50);
    [Minima,P]=islocalmin(N);
    threshold=X(P==max(P));
    threshold=threshold(1);
    if threshold>550
        threshold=550;
    end
    Idx_area=(stats.Area<threshold);
    Opt_LBW(cat(1,CC.PixelIdxList{Idx_area}))=false;
    Data=[stats.Area sigma idx];
    
end