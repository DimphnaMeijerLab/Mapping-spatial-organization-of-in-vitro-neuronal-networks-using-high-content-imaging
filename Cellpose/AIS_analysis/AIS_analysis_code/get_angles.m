function [Anglesperneuron,Angles,properties,NumAis,EndVectorsperneuron,BeginVectorsperneuron]=get_angles(Vectors,Somas,Distances,lengths)
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
    conversion = 443.14 / 1104;   % µm per pixel ≈ 0.4016
    count=1;
    ccount=1;
    for i=1:length(Anglesperneuron(:,1))
        if isempty(Anglesperneuron{i,1})==0
            for j=1:length(Anglesperneuron{i,:})
                Angles(count,1)=Anglesperneuron{i,1}(j,1);
                Angles(count,2)=(Center.Centroid(i,1)-mid(1,1))*conversion;
                Angles(count,3)=(Center.Centroid(i,2)-mid(1,2))*conversion;
                count=count+1;
            end
            for j=1:length(Distances{i,:})    
                properties(ccount,1)=Distances{i,1}(j,1)*conversion;
                properties(ccount,2)=lengths{i,1}(j,1)*conversion;
                properties(ccount,3)=(Center.Centroid(i,1)-mid(1,1))*conversion;
                properties(ccount,4)=(Center.Centroid(i,2)-mid(1,2))*conversion;
                ccount=ccount+1;
            end
            
            
        end
    end    
end 