clc
clear all
close all

disp('Select the data folder')
analysepath=uigetdir;
wells=["B08"];

%%
%create parallel pool
parpool('local');

for i=1:length(wells)
    well=wells(i);


    SomaMask = imread(strcat(analysepath,'\',well,'\',well,'_SomaMasks.tif'));
    AnkGImage = imread(strcat(analysepath,'\',well,'\',well,'_NeuronMask.tif'));
    [Allgood, AISMask] = segmentAIS_thresholdMerge(SomaMask, AnkGImage);

   
    %% 
    
    
    Allgood=int16(Allgood);
    AnkG=imread(strcat(analysepath,'\',well,'\',well,'_AnkG.tif'));
    
    tic
    Somas=findSomas(Allgood);
    toc
%% 

    tic
    [Vectors,AISparts,Distances_to_Soma,Max_Ais_pos,AISlength]=AISfinding(Allgood,Somas,AnkG);
    toc
   
    tic
    [Angles_neuron,Allangles,NumAis,Vectorsperneuron]=get_angles(Vectors,Somas,Distances_to_Soma,AISlength);
    toc

    %save results
    AIS_Results.Vectors=Vectors;
    AIS_Results.NumAIS=NumAis;
    AIS_Results.Angles_neuron=Angles_neuron;
    AIS_Results.Allangles=Allangles;
    AIS_Results.Somas=Somas;
    AIS_Results.Vectorsperneuron=Vectorsperneuron;
    AIS_Results.Distances=Distances_to_Soma;
    AIS_Results.AISlength=AISlength;
    save(strcat(analysepath,'\',well,'\AISResults.mat'),'AIS_Results');
end
%end parallel pool
delete(gcp('nocreate'));

%% visualization
count=1;
for i=1:length(Max_Ais_pos)
    if isempty(Max_Ais_pos{i,1})==0
        for j=1:length(Max_Ais_pos{i,1}(:,1))
            pos(count,1)=Max_Ais_pos{i,1}(j,1);
            pos(count,2)=Max_Ais_pos{i,1}(j,2);
            count=count+1;
        end
    end
end
count=1;
for i=1:length(Vectors(:,1))
    if isempty(Vectors{i})==0
        for j=1:length(Vectors(i,:))
            if isempty(Vectors{i,j})==0
                for k=1:length(Vectors{i,j}(:,1))
                    for m=1:length(Vectors{i,j}{k,1})
                        newVectors(count,1)=Vectors{i,j}{k,1}(1,m);
                        newVectors(count,2)=Vectors{i,j}{k,2}(1,m);
                        newVectors(count,3)=Vectors{i,j}{k,3}(1,m);
                        newVectors(count,4)=Vectors{i,j}{k,4}(1,m);
                        count=count+1;
                    end
                end
            end
        end     
    end
end
%Coloredsegmented=imread(strcat(analysepath,'\SegmentedImage_rgb.tif'));
figure(),imshow(AnkG), hold on;
plot(pos(:,1),pos(:,2),'*r'), hold on;
quiver(newVectors(:,1),newVectors(:,2),newVectors(:,3),newVectors(:,4),0,'b'),
hold off;