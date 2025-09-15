function [Glia_stats, Neuron_stats, LBW_Glia, LBW_Neurons, threshold]=Copy_2_of_Nuclei_analysis(LBW_Dapi,LBW_AnkG,Dapi)
    %% 
    Neur_size=500; 
    overlap=0.5;
    %Remove rubbish and watershed Nuclei
    LBW_Dapi=firstfilter(LBW_Dapi,[200,3800],1);

    LBW_AnkG=firstfilter(LBW_AnkG,[Neur_size,Inf],0);
    %% 

    %filter on overlap
    LBW_Neurons=overlapfilter(LBW_Dapi, LBW_AnkG, overlap);
    %%
    %size and intensity filter by first kmeans clustering on Area size and std of
    %intensity, afterwards a threshold is determined by finding the largest
    %minimum in the area histogram
    [LBW_Neurons, ~, threshold]=brightsizefilteralt(LBW_Neurons,Dapi);
    
    %stats second filter step
    CC=bwconncomp(LBW_Neurons);
    Neurons_stats_second=regionprops('table',CC,Dapi,'Area','MeanIntensity','Centroid');
    % Glia cells
    LBW_Glia=logical(LBW_Dapi-LBW_Neurons);
    CC=bwconncomp(LBW_Glia);
    Gliapos=regionprops('table',CC,Dapi,'Centroid','Area','MeanIntensity');

    %new reference frame with mid point as (0,0) 
    Mid=size(Dapi)/2;
    Glia_Car=[Gliapos.Centroid(:,1)-Mid(1,1) Gliapos.Centroid(:,2)-Mid(1,2)]; 
    Neurons_Car=[Neurons_stats_second.Centroid(:,1)-Mid(1,1) Neurons_stats_second.Centroid(:,2)-Mid(1,2)]; 

    Glia_stats = struct('Centroid', Glia_Car, 'Area', Gliapos.Area, 'Mean_intensity', Gliapos.MeanIntensity);
    Neuron_stats = struct('Centroid', Neurons_Car, 'Area', Neurons_stats_second.Area, 'Mean_intensity', Neurons_stats_second.MeanIntensity);
    
%     Glia_stats=[Glia_Car(:,1), Glia_Car(:,2), Gliapos.Area, Gliapos.MeanIntensity] ;
%     Neuron_stats=[Neurons_Car(:,1), Neurons_Car(:,2), Neurons_stats_second.Area, Neurons_stats_second.MeanIntensity];
   
end

