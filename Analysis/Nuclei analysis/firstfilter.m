function filtered=firstfilter(blur,range,watershed)
    %% 
    filtered=blur;
    if watershed==1
        filtered=Watershed(filtered);
    end
    filtered=bwareaopen(filtered,range(1));
    CC=bwconncomp(filtered);
    stats = regionprops('table',CC,'Area');
    removeMask=(stats.Area)>range(2);
    filtered(cat(1,CC.PixelIdxList{removeMask}))=false;
end
