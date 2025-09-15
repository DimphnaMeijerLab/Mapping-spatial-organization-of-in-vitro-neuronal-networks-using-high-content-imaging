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
        parfor j=1:length(x)-1 %can be parfor
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