function [Vectors,AllAIS,Distances_to_Soma,Max_Ais_pos,AISlength]=AISfinding(Allgood,Somas,AnkG)
    %this function finds the AIS objects in the segmented image from cellpose
    %based on area and circularity. Furthermore, it finds the distance of
    %the AIS to the soma and creates vectors following the path of the AIS
    
    %=====variables used=======%
    %these are used to find the distance along the AIS
    Number=max(Somas,[],'all');
    sz=size(Allgood);
    maxAISnumber=5;
    D=bwdistgeodesic(logical(Allgood),logical(Somas));
    D(isnan(D))=0;
    %these are used to cut up and put back together the image
    Centers=regionprops('table',Allgood,'Centroid');
    Centers=Centers.Centroid;
    %======output variables=======%
    Vectors=cell(Number,1);
    AllAIS=cell(Number,1);
    Distances_to_Soma=cell(Number,1);
    Max_Ais_pos=cell(Number,1);
    AISlength=cell(Number,1);
    disp('Classify AIS...')
    
    for i = 1:Number   %can be parfor
        %used variables per round
        Value=i;
        perimetervalue=200;
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
        
        if sum(Soma,'all')>0
            %make sure one Soma object
            cc=bwconncomp(Soma);
            if cc.NumObjects>0
                stats=regionprops('table',cc,'Area');
                idx=find(max(stats.Area)==stats.Area);
                Soma=ismember(labelmatrix(cc),idx);
            end

            fakecenter=regionprops('table',bwconncomp(SingleAISneuron),'Centroid'); 
            fakecenter=fakecenter.Centroid;

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
                AISlength{i,1}=zeros(CCAISobjects.NumObjects,1);
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
                        AISlength{i,1}(j,1)=max(AIS,[],'all')-Min;
                        Distances_to_Soma{i,1}(j,1)=sqrt((mean(col)-fakecenter(1,1))^2+(mean(row)-fakecenter(1,2))^2);
                        
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
                                test=zeros(intersz); test(Extrema(:,2),Extrema(:,1))=1; test=imdilate(test,ones(9));
                                Dtest=test.*AISD;
                                cctest=bwconncomp(logical(Dtest));
                                stats=regionprops('table',cctest,Dtest,'MeanIntensity');
                                idx=find(stats.MeanIntensity>=min(maxk(stats.MeanIntensity,2)));
                                test=ismember(labelmatrix(cctest),idx);
                                cctest=bwconncomp(test);
                                ypos{1,k}=zeros(1,cctest.NumObjects); xpos{1,k}=zeros(1,cctest.NumObjects);
                                if cctest.NumObjects>2
                                    disp('yes')
                                end
                                for m=1:cctest.NumObjects
                                    [Xend,Yend]=ind2sub(intersz,cctest.PixelIdxList{1,m});
                                    ypos{1,k}(1,m)=mean(Xend);
                                    xpos{1,k}(1,m)=mean(Yend);
                                end
                                break
                            end
                        end
                        %store points in vectors
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

 
