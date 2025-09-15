clc
clear all
close all

disp('Select the folder with the cellpose segmentation')
analysepath=uigetdir;

Allgood=imread(strcat(analysepath,'\C02_Cellpose_RGB_cellpose_segmentation_custom_model.tif'));
Allgood=int16(Allgood);
AnkG=imread(strcat(analysepath,'\C02_AnkG.tif'));

%create parallel pool
parpool('local');

tic
Somas=findSomas(Allgood);
toc

tic
[Vectors,AISparts,Distances,Max_Ais_pos]=AISfinding(Allgood,Somas,AnkG);
toc


tic
[Angles_neuron,Allangles,NumAis,Vectorsperneuron]=get_angles(Vectors,Somas);
toc

%end parallel pool
delete(gcp('nocreate'));

%save results
% AIS_Results.Vectors=Vectors;
% AIS_Results.AISparts=AISparts;
% AIS_Results.NumAIS=NumAIS;
% AIS_Results.Angles_neuron=Angles_neuron;
% AIS_Results.Allangles=Allangles;
% AIS_Results.Somas=Somas;
% AIS_Results.Distances=Distances_to_Soma;
% save(strcat(analysepath,'\AISResults',well,'.mat'),'AIS_Results');

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
%saveas(gcf,fullfile(analysepath, 'vectorPlot.tif'))

%% Functions

function Somas=findSomas(Allgood)
    % This function finds Soma based on the cellpose segmentation with the
    % use of the Soma extraction algorithm
    % (https://github.com/cihanbilge/SomaExtraction). Note that only the
    % first half is used, in which the directionality of objects is found.
    % Afterwards I have set the threshold for that measure to 0.60 to find
    % Soma.
    
    % get Soma
    folder=pwd;
    addpath(...
        fullfile(folder),...
        fullfile(folder, 'anigauss'));
    %parameters, standard
    filtersize=9; 
    num_direction=10;
    disp('Detecting soma...')

    %split the picture
    sz=size(Allgood);
    Somas=zeros(sz);
    x=round(0:sz(1)/9:sz(1));
    Intersomas=cell(length(x)-1,length(x)-1);
    Intersomaspos=cell(length(x)-1,length(x)-1);
    for i=1:length(x)-1
        parfor j=1:length(x)-1
            dx=200;
            xlow=max(x(i)-dx,1);
            ylow=max(x(j)-dx,1);
            xhigh=x(i+1)+dx; xhigh(xhigh>sz(1,2))=sz(1,2);
            yhigh=x(j+1)+dx; yhigh(yhigh>sz(1,1))=sz(1,1);
            xrange=xlow:xhigh; yrange=ylow:yhigh;
            Inter=Allgood(yrange,xrange);
            [dirRatio, ~]=Main_Anigauss_2d(Inter>0,filtersize,num_direction);
            Intersomas{i,j}=dirRatio>0.60;
            Intersomaspos{i,j}={xrange,yrange};
        end
        for j=1:length(x)-1
            Somas(Intersomaspos{i,j}{1,2},Intersomaspos{i,j}{1,1})=(Somas(Intersomaspos{i,j}{1,2},Intersomaspos{i,j}{1,1})+Intersomas{i,j})>0;
        end
    end          
        
%     [dirRatio, ~]=Main_Anigauss_2d(Allgood>0,filtersize,num_direction);
%     Somas=dirRatio>0.60;
    Somas=Somas.*double(Allgood);    
    Somas=imerode(Somas>0,ones(3)); Somas=imdilate(Somas>0,ones(3));
    Somas=bwareaopen(Somas,200);
    Somas=Somas.*double(Allgood);
end

function [Vectors,AllAIS,Distances_to_Soma,Max_Ais_pos]=AISfinding(Allgood,Somas,AnkG)
    %this function finds the AIS objects in the segmented image from cellpose
    %based on area and circularity. Furthermore, it finds the distance of
    %the AIS to the soma and creates vectors following the path of the AIS
    
    %=====variables used=======%
    Number=max(Allgood,[],'all');
    sz=size(Allgood);
    maxAISnumber=5;
    D=bwdistgeodesic(logical(Allgood),logical(Somas));
    D(isnan(D))=0;
    Centers=regionprops('table',Allgood,'Centroid');
    Centers=Centers.Centroid;
    %======output variables=======%
    Vectors=cell(Number,1);
    AllAIS=cell(Number,1);
    Distances_to_Soma=cell(Number,1);
    Max_Ais_pos=cell(Number,1);
    disp('Classify AIS...')
    count=0;
    for i = 1:Number
        %used variables per round
        Value=i;
        perimetervalue=150;
        xlow=round(max(Centers(i,1)-perimetervalue,1));
        ylow=round(max(Centers(i,2)-perimetervalue,1));
        xhigh=round(Centers(i,1)+perimetervalue); xhigh(xhigh>sz(2))=sz(2);
        yhigh=round(Centers(i,2)+perimetervalue); yhigh(yhigh>sz(1))=sz(1);
       
        
        %create smaller intermediary images to save memory
        SingleAISneuron=(Allgood(ylow:yhigh,xlow:xhigh)==Value);
        Soma=(Somas(ylow:yhigh,xlow:xhigh)==Value);
        Dinter=D(ylow:yhigh,xlow:xhigh);
        Ankginter=AnkG(ylow:yhigh,xlow:xhigh);
        intersz=size(Soma);

        %make sure one Soma object
        cc=bwconncomp(Soma);
        if cc.NumObjects>0
            stats=regionprops('table',cc,'Area');
            idx=find(max(stats.Area)==stats.Area);
            Soma=ismember(labelmatrix(cc),idx);
        end
        
        fakecenter=regionprops('table',bwconncomp(SingleAISneuron),'Centroid'); fakecenter=fakecenter.Centroid;
        
        %classify AIS objects based on size (not too small) and circularity
        %(should be fairly low, since an AIS is a linear structure).
        AISobjects=SingleAISneuron-Soma; 
        AISobjects=bwareaopen(AISobjects,75);
        CCAISobjects=bwconncomp(AISobjects);
        if CCAISobjects.NumObjects>0
            stats=regionprops('table',CCAISobjects,'Area','Circularity');
            idx=find(((stats.Area<150).*(stats.Circularity>0.5))+((stats.Area>500).*(stats.Circularity>0.3)));
            fakeAISobjects=ismember(labelmatrix(CCAISobjects),idx);
            realAISobjects=AISobjects-fakeAISobjects;
            CCAISobjects=bwconncomp(realAISobjects);
        end
        
        %find vectors and distance to soma
        if CCAISobjects.NumObjects>0
            Intensitystats=regionprops('table',CCAISobjects,Ankginter,'PixelValues');
            Distances_to_Soma{i,1}=zeros(CCAISobjects.NumObjects,1);
            Max_Ais_pos{i,1}=zeros(CCAISobjects.NumObjects,2);
            
            if sum(sum(Soma))>0
                for j=1:maxAISnumber
                    if j<=CCAISobjects.NumObjects
                        %Find Distance of peak AnkG, assumed to be center
                        %of AIS. Further center means more activity.
                        Toppixel=mean(maxk(Intensitystats.PixelValues{j,1},50));
                        Toppixel=Toppixel-0.1*Toppixel;
                        idx=[Intensitystats.PixelValues{j,1}>Toppixel];
                        Distances=Dinter(CCAISobjects.PixelIdxList{1,j}(idx,1));
                        AIS=zeros(intersz); AIS(CCAISobjects.PixelIdxList{1,j})=1;
                        AIS=AIS.*Dinter;
                        Min=min(Distances);
                        [row,col]=find(AIS==Min);
                        Max_Ais_pos{i,1}(j,:)=[mean(col)-fakecenter(1,1)+Centers(i,1),mean(row)-fakecenter(1,2)+Centers(i,2)];
                        Distances_to_Soma{i,1}(j,1)=Min;
                        
                        %Find direction of each AIS
                        partAIS=zeros(intersz); partAIS(CCAISobjects.PixelIdxList{1,j})=1;                        
                        AISD=partAIS.*Dinter;
                        AISD(isnan(AISD))=0;
                        startvalue=min(min(AISD(AISD>0)));
                        Endvalue=max(AISD,[],'all');
                        nsteps=ceil((Endvalue-startvalue)/10)+1;
                        intermediaries=round(linspace(startvalue,Endvalue,nsteps));
                        ypos=cell(1,nsteps); xpos=cell(1,nsteps);
                        for k=1:nsteps
                            test=zeros(intersz);
                            [x,y]=find(AISD==intermediaries(1,k));                            
                            test(x,y)=1; test=imdilate(test,ones(5));
                            cc=bwconncomp(test);
                            if cc.NumObjects==1
                                ypos{1,k}=mean(x); xpos{1,k}=mean(y);
                              
                            elseif cc.NumObjects>1  %if ais splits
                                stats=regionprops('table',partAIS,'Extrema');
                                Extrema=round(stats.Extrema{1,1});
                                test=zeros(intersz); test(Extrema(:,2),Extrema(:,1))=1; test=imdilate(test,ones(3));
                                Dtest=test.*AISD;
                                cctest=bwconncomp(logical(Dtest));
                                stats=regionprops('table',cctest,Dtest,'MeanIntensity');
                                idx=find(stats.MeanIntensity>=min(maxk(stats.MeanIntensity,cc.NumObjects)));
                                test=ismember(labelmatrix(cctest),idx);
                                cctest=bwconncomp(test);
                                ypos{1,k}=zeros(1,cctest.NumObjects); xpos{1,k}=zeros(1,cctest.NumObjects);
                                for m=1:cctest.NumObjects
                                    [Xend,Yend]=ind2sub(intersz,cctest.PixelIdxList{1,m});
                                    ypos{1,k}(1,m)=mean(Xend);
                                    xpos{1,k}(1,m)=mean(Yend);
                                end
                                break
                            end
                        end
                        xpos=xpos(~cellfun('isempty',xpos)); ypos=ypos(~cellfun('isempty',ypos));
                        if length(xpos)>1
                        Vectors{i,j}=cell(length(xpos)-1,4);
                        for k=1:length(xpos)-2
                            Vectors{i,j}{k,1}=xpos{1,k}-fakecenter(1,1)+Centers(i,1);
                            Vectors{i,j}{k,2}=ypos{1,k}-fakecenter(1,2)+Centers(i,2);
                            Vectors{i,j}{k,3}=xpos{1,k+1}-xpos{1,k};
                            Vectors{i,j}{k,4}=ypos{1,k+1}-ypos{1,k};
                        end
                        if length(xpos{1,end})>1
                            Vectors{i,j}{end,1}=zeros(1,length(xpos{1,end}));
                            Vectors{i,j}{end,2}=zeros(1,length(xpos{1,end}));
                            Vectors{i,j}{end,3}=zeros(1,length(xpos{1,end}));
                            Vectors{i,j}{end,4}=zeros(1,length(xpos{1,end}));
                            for k=1:length(xpos{1,end})
                                Vectors{i,j}{end,1}(1,k)=xpos{1,end-1}-fakecenter(1,1)+Centers(i,1);
                                Vectors{i,j}{end,2}(1,k)=ypos{1,end-1}-fakecenter(1,2)+Centers(i,2);
                                Vectors{i,j}{end,3}(1,k)=xpos{1,end}(1,k)-xpos{1,end-1};
                                Vectors{i,j}{end,4}(1,k)=ypos{1,end}(1,k)-ypos{1,end-1};
                            end
                        else
                            Vectors{i,j}{end,1}=xpos{1,end-1}-fakecenter(1,1)+Centers(i,1);
                            Vectors{i,j}{end,2}=ypos{1,end-1}-fakecenter(1,2)+Centers(i,2);
                            Vectors{i,j}{end,3}=xpos{1,end}-xpos{1,end-1};
                            Vectors{i,j}{end,4}=ypos{1,end}-ypos{1,end-1};
                        end
                        end
                     end
                 end
            end
        end
    end  
end

function [Anglesperneuron,Angles,NumAis,EndVectorsperneuron,BeginVectorsperneuron]=get_angles(Vectors,Somas)
    %This function finds the angle of all the AIS of each neuron w.r.t. a
    %vector pointing to the middle of the well. Additionally, it returns
    %the number of AIS for each neuron.
    
    %======variables used=====%
    sz=size(Somas);
    NeuronNumber=length(Vectors(:,1));
    mid=sz/2;
    Center=regionprops('table',Somas,'Centroid');
    
    %====output=====%
    NumAis=zeros(NeuronNumber,1);
    Anglesperneuron=cell(NeuronNumber,1);
    EndVectorsperneuron=cell(NeuronNumber,1);
    BeginVectorsperneuron=cell(NeuronNumber,1);
    
    %=====fuction=====%
    for i=1:NeuronNumber
        x=length(Vectors(i,:));
        AISnumber=0;
        for j=1:x
            if isempty(Vectors{i,j})==0
                AISnumber=AISnumber+1;
                if length(Vectors{i,j}{end,1})>1
                    AISnumber=AISnumber+length(Vectors{i,j}{end,1})-1;
                end
            end
        end
        NumAis(i,1)=AISnumber;
        count=1;        
        NeuronAngles=zeros(AISnumber,1);
        EndNeuronvectors=cell(AISnumber,1);
        BeginNeuronvectors=cell(AISnumber,1);
        CenterV=[mid(1,1)-Center.Centroid(i,1), mid(1,2)-Center.Centroid(i,2)];
%         Midvectors(i,1)=Center.Centroid(i,1); Midvectors(i,2)=Center.Centroid(i,2);  
%         Midvectors(i,3)=CenterV(1,1); Midvectors(i,4)=CenterV(1,2); 
        
        for j=1:length(Vectors(i,:))
            if isempty(Vectors{i,j})==0
                BeginNeuronvectors=[Vectors{i,j}{1,3}, Vectors{i,j}{1,4}];
                for k=1:length(Vectors{i,j}{end,3})
                    V=[Vectors{i,j}{end,3}(1,k), Vectors{i,j}{end,4}(1,k)];
                    EndNeuronvectors{count,1}=V;
                    NeuronAngles(count,1)=acos(dot(V,CenterV)/(norm(V)*norm(CenterV)))*(360/2/pi);
                    count=count+1;
                end
            end
        end
        Anglesperneuron{i,1}=NeuronAngles;
        EndVectorsperneuron{i,1}=EndNeuronvectors;
        BeginVectorsperneuron{i,1}=BeginNeuronvectors;
    end
    
    count=1;
    for i=1:length(Anglesperneuron(:,1))
        if isempty(Anglesperneuron{i,1})==0
            for j=1:length(Anglesperneuron{i,:})
                Angles(count,1)=Anglesperneuron{i,1}(j,1);
                count=count+1;
            end
        end
    end    
end    

