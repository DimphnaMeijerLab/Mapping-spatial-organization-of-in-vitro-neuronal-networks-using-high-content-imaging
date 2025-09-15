function Opt_LBW=overlapfilter(LBW,mask,percentage)
        %% 
        Opt_LBW=LBW;
        Overlap=mask.*Opt_LBW;
        CC=bwconncomp(Opt_LBW);
        Overlap=Opt_LBW+Overlap;
        stats_Overlap=regionprops('table',CC,Overlap,'MeanIntensity');
        removeMask=(stats_Overlap.MeanIntensity-1)<percentage;
        Opt_LBW(cat(1,CC.PixelIdxList{removeMask}))=false;      
end