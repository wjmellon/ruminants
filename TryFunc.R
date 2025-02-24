### Longevity model
#longevity neo
AnalysisFuncion <- function(RunName = "",KeepCol=c(5,9,10,11,42),whattree="",var1,var2) {
  cutData <- DataRum[,KeepCol,drop=FALSE] 
  cutData[cutData$var1 < 0,] <-NA
  cutData <- na.omit(cutData)
  tree <- read.tree(whattree)
  
  view(cutData)
  
  cutData$Species <- gsub(" ", "_", cutData$Species) 
  cutData$common_name<-gsub("_", "", cutData$common_name)
  includedSpecies<-cutData$Species
  pruned.tree<-drop.tip(
    tree, setdiff(
      tree$tip.label, includedSpecies))
  pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
  cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
  cutData <- cutData[!(cutData$Keep==FALSE),]
  rownames(cutData)<-cutData$Species
  SE<-setNames(cutData$SE_simple,cutData$Species)[rownames(cutData)]
  view(cutData)
  
  #pgls model
  if (nrow(cutData) > 9) {
    LocalPgls<-pglsSEyPagel(var2~var1,data=cutData,
                                tree=pruned.tree,method="ML",se=SE)
    summary(LocalPgls)
    
    #grab r squared, lambda, and p values from summary 
    
    r.v.LocalPgls <- R2(phy = pruned.tree,LocalPgls)
    r.v.LocalPgls <- format(r.v.LocalPgls[3])
    r.v.LocalPgls <-signif(as.numeric(r.v.LocalPgls), digits= 2)
    ld.v.LocalPgls<- summary(LocalPgls)$modelStruct$corStruct
    ld.v.LocalPgls <- signif(ld.v.LocalPgls[1], digits = 2)
    p.v.LocalPgls<-summary(LocalPgls)$tTable
    p.v.LocalPgls<-signif(p.v.LocalPgls[2,4], digits = 3)
    
    longneo<-ggplot(cutData, aes(y=NeoplasiaPrevalence*100, x=log10(max_longevity.months.))) + 
      scale_color_manual(values = c("Mammalia" = "#EB1D24", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff" ))+
      scale_y_continuous(
        limits = c(0,75),
        breaks = c(0, 25,50,75),
        labels = c(0, 25,50,75))+
      geom_abline(intercept = coef(longevity.neo)[1]*100, slope =  coef(longevity.neo)[2]*100,
                  color = 'grey',size = 1.2) +
      theme_cowplot(12)+
      theme(axis.title = element_text(size = 18))+
      ylab("Neoplasia Prevalence (%)") +
      xlab("(log10) Max Longevity (Mo)") +
      geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
      scale_size(name   = "Total Necropsies",
                 breaks = c(20,100,200,300,477),
                 labels =  c(20,100,200,300,477))+
      geom_text_repel(aes(label=ifelse(NeoplasiaPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
      guides(colour = guide_legend(override.aes = list(size=5))) +
      labs(title = "B")+
      theme(
        plot.title = element_text(size = 20, face = "bold")) +
      theme(legend.position = "bottom")+
      labs(colour="Clade", size="Total Necropsies")
    
    
    #ggsave(filename='longneo.pdf', width=13, height=10, limitsize=FALSE,bg="white")
    
    #rbind adds a new row to a dataframe
    SummaryStats <- rbind(SummaryStats, list("Neo vs Longevity",nrow(cutData), r.v.longneo, p.v.longneo, ld.v.longneo))
    
    #Print stats
    print(SummaryStats)
  }
}