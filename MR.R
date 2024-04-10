#1
#install.packages("devtools")
#devtools::install_github("mrcieu/ieugwasr", force = TRUE)


#???ð?
library(ieugwasr)

inputFile="exposure_data.csv"      #?????ļ?
setwd("C:\\biowolf\\proteinMR\\05.F")     #???ù???Ŀ¼

#??ȡ?????ļ?
dat=read.csv(inputFile, header=T, sep=",", check.names=F)

#????F????ֵ
dat$R2<-(2*dat$beta.exposure*dat$beta.exposure*dat$eaf.exposure*(1-dat$eaf.exposure)/(2*dat$beta.exposure*dat$beta.exposure*dat$eaf.exposure*(1-dat$eaf.exposure)+2*dat$se.exposure*dat$se.exposure*dat$samplesize.exposure*dat$eaf.exposure*(1-dat$eaf.exposure)))     #????R2
dat$F<-dat$R2*(dat$samplesize.exposure-2)/(1-dat$R2)     #????F????ֵ

#????Fֵ>10?????ݽ??й???, ɾ???????߱?��
outTab=dat[as.numeric(dat$F)>10,]
write.csv(outTab, file="exposure.F.csv", row.names=F)



#2
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
BiocManager::install("VariantAnnotation")

#install.packages("devtools")
devtools::install_github("mrcieu/gwasglue", force = TRUE)

#install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")



#3
#引用包
library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)

exposureFile="exposure.F.csv"            #暴露数据文件
outcomeID="ukb-a-114"                      #结局数据id(需修改)
outcomeName="Hypertrophic cardiomyopathy"     #图形中展示疾病的名称
setwd("C:\\biowolf\\proteinMR\\06.MR")     #设置工作目录

#读取暴露数据
exposure_dat=read_exposure_data(filename=exposureFile,
                                sep = ",",
                                snp_col = "SNP",
                                beta_col = "beta.exposure",
                                se_col = "se.exposure",
                                pval_col = "pval.exposure",
                                effect_allele_col="effect_allele.exposure",
                                other_allele_col = "other_allele.exposure",
                                eaf_col = "eaf.exposure",
                                phenotype_col = "exposure",
                                id_col = "id.exposure",
                                samplesize_col = "samplesize.exposure",
                                chr_col="chr.exposure", pos_col = "pos.exposure",
                                clump=FALSE)

#读取结局数据
outcomeData=extract_outcome_data(snps=exposure_dat$SNP, outcomes=outcomeID)
write.csv(outcomeData, file="outcome.csv", row.names=F)

#将暴露数据和结局数据合并
outcomeData$outcome=outcomeName
dat=harmonise_data(exposure_dat, outcomeData)

#输出用于孟德尔随机化的工具变量
outTab=dat[dat$mr_keep=="TRUE",]
write.csv(outTab, file="table.SNP.csv", row.names=F)

#孟德尔随机化分析
mrResult=mr(dat)

#对结果进行OR值的计算
mrTab=generate_odds_ratios(mrResult)
write.csv(mrTab, file="table.MRresult.csv", row.names=F)

#异质性分析
heterTab=mr_heterogeneity(dat)
write.csv(heterTab, file="table.heterogeneity.csv", row.names=F)

#多效性检验
pleioTab=mr_pleiotropy_test(dat)
write.csv(pleioTab, file="table.pleiotropy.csv", row.names=F)

#绘制散点图
pdf(file="pic.scatter_plot.pdf", width=7.5, height=7)
mr_scatter_plot(mrResult, dat)
dev.off()

#森林图
res_single=mr_singlesnp(dat)      #得到每个工具变量对结局的影响
pdf(file="pic.forest.pdf", width=7, height=5.5)
mr_forest_plot(res_single)
dev.off()

#漏斗图
pdf(file="pic.funnel_plot.pdf", width=7, height=6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()

#留一法敏感性分析
pdf(file="pic.leaveoneout.pdf", width=7, height=5.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
dev.off()

pvalFilter=0.05      #pvalue过滤条件
fdrFilter=0.05        #FDR过滤条件
mrFile="table.MRresult.csv"        #孟德尔随机化分析的结果文件
pleFile="table.pleiotropy.csv"     #多效性的结果文件
setwd("C:\\biowolf\\proteinMR\\07.IVWfilter")     #???ù???Ŀ¼

#读取孟德尔随机化结果文件
rt=read.csv(mrFile, header=T, sep=",", check.names=F)

#对pvalue进行校正
ivwRT=rt[rt$method=="Inverse variance weighted",]
ivwPval=ivwRT[,"pval"]
fdr=p.adjust(ivwPval, method="fdr")
ivwRT=cbind(ivwRT, fdr)

#提取5种方法OR值方向一致的蛋白
ivw=data.frame()
for(protein in unique(ivwRT$exposure)){
  proteinData=rt[rt$exposure==protein,]
  if(sum(proteinData$or>1)==nrow(proteinData) | sum(proteinData$or<1)==nrow(proteinData)){
    ivw=rbind(ivw, ivwRT[ivwRT$exposure==protein,])
  }
}

#读取多效性的结果文件
pleRT=read.csv(pleFile, header=T, sep=",", check.names=F)
#剔除pvalue小于0.05的蛋白
pleRT=pleRT[pleRT$pval>0.05,]
proteinLists=as.vector(pleRT$exposure)
outTab=ivw[ivw$exposure %in% proteinLists,]
outTab=outTab[order(outTab$pval),]
outTab=outTab[!duplicated(outTab$exposure),]
write.csv(outTab, file="IVW.result.csv", row.names=F)

#根据pvalue和fdr设置的条件对其进行过滤
outTab=outTab[((outTab$pval<pvalFilter) & (outTab$fdr<fdrFilter)),]
write.csv(outTab, file="IVW.filter.csv", row.names=F)



#4

install.packages("ggplot2")
install.packages("ggrepel")


#引用包
library(dplyr)
library(ggplot2)
library(ggrepel)

diffFile="IVW.resultshan.csv"       #IVW方法的结果文件
geneFile="IVW.filtershan.csv"       #IVW方法过滤的结果文件
setwd("C:\\biowolf\\proteinMR\\08.vol")      #???ù???Ŀ¼

#读取IVW方法的结果文件
rt=read.csv(diffFile, header=T, sep=",", check.names=F)
rt=rt[order(rt$pval),]
row.names(rt)=rt[,"exposure"]

#读取显著蛋白列表文件
geneRT=read.csv(geneFile, header=T, sep=",", check.names=F)
showGene=as.vector(geneRT[,"exposure"])
showData=rt[showGene,]

#定义显著性
Risk=ifelse(rt$exposure %in% showGene, ifelse(rt$b>0,"High","Low"), "Not")
Risk=factor(Risk, levels=c("Low", "Not", "High"))

#绘制火山图
rt = mutate(rt, Risk=Risk)
p = ggplot(rt, aes(b, -log10(pval)))+
  geom_point(aes(col=Risk))+ 
  xlab("Beta") + ylab("-log10(pvalue)") + xlim(-1,1)+
  scale_color_manual(values=c("green", "grey","red"))+
  labs(title = " ")+
  theme(plot.title = element_text(size=16, hjust=0.5, face = "bold"))

#在图形中表注蛋白名称
p1=p+geom_text_repel(data=showData,
                     box.padding=0.2, point.padding=0.15, min.segment.length=0.01,
                     size=3, aes(label=exposure)) + theme_bw()
#输出图形
pdf(file="vol.pdf", width=5.5, height=4.5)
print(p1)
dev.off()



#5
#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")


#引用包
library(TwoSampleMR)

exposureFile="exposure.F.csv"            #暴露数据文件
outcomeFile="outcome.csv"                #结局数据文件
sigExpoFile="IVW.filter.csv"             #IVW方法过滤的结果文件
outcomeName="Atrial fibrillation"     #设置图形中展示疾病的名字
setwd("C:\\biowolf\\proteinMR\\09.MRpic")     #设置工作目录

#读取暴露数据输入文件
rt=read.csv(exposureFile, header=T, sep=",", check.names=F)

#读取IVW方法过滤的结果文件
sigExpo=read.csv(sigExpoFile, header=T, sep=",", check.names=F)

#对疾病相关的蛋白进行循环
for(protein in unique(sigExpo$id.exposure)){
  #提取这个蛋白的暴露数据
  i=gsub("(.+?)\\_(.+?)\\_(.+?)\\_(.+)", "\\3", protein)
  singleExposureFile=paste0(i, ".exposure.csv")
  exposure_set=rt[rt$id.exposure==protein,]
  write.csv(exposure_set, file=singleExposureFile, row.names=F)
  #读取这个蛋白的暴露数据
  exposure_dat=read_exposure_data(filename=singleExposureFile,
                                  sep = ",",
                                  snp_col = "SNP",
                                  beta_col = "beta.exposure",
                                  se_col = "se.exposure",
                                  pval_col = "pval.exposure",
                                  effect_allele_col="effect_allele.exposure",
                                  other_allele_col = "other_allele.exposure",
                                  eaf_col = "eaf.exposure",
                                  phenotype_col = "exposure",
                                  id_col = "id.exposure",
                                  samplesize_col = "samplesize.exposure",
                                  chr_col="chr.exposure", pos_col = "pos.exposure",
                                  clump=FALSE)
  
  #读取结局数据
  outcome_data=read_outcome_data(snps=exposure_dat$SNP,
                                 filename="outcome.csv", sep = ",",
                                 snp_col = "SNP",
                                 beta_col = "beta.outcome",
                                 se_col = "se.outcome",
                                 effect_allele_col = "effect_allele.outcome",
                                 other_allele_col = "other_allele.outcome",
                                 pval_col = "pval.outcome",
                                 eaf_col = "eaf.outcome")
  
  #将暴露数据和结局数据合并
  outcome_data$outcome=outcomeName
  dat=harmonise_data(exposure_dat, outcome_data)
  
  #输出用于孟德尔随机化的工具变量
  outTab=dat[dat$mr_keep=="TRUE",]
  write.csv(outTab, file=paste0(i, ".table.SNP.csv"), row.names=F)
  
  #MR-PRESSO异常值检测(偏倚的SNP)
  presso=run_mr_presso(dat)
  write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file=paste0(i, ".table.MR-PRESSO_Global.csv"))
  write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file=paste0(i, ".table.MR-PRESSO_Outlier.csv"))
  
  #孟德尔随机化分析
  mrResult=mr(dat)
  
  #对结果进行OR值的计算
  mrTab=generate_odds_ratios(mrResult)
  write.csv(mrTab, file=paste0(i, ".table.MRresult.csv"), row.names=F)
  
  #异质性检验
  heterTab=mr_heterogeneity(dat)
  write.csv(heterTab, file=paste0(i, ".table.heterogeneity.csv"), row.names=F)
  
  #多效性检验
  pleioTab=mr_pleiotropy_test(dat)
  write.csv(pleioTab, file=paste0(i, ".table.pleiotropy.csv"), row.names=F)
  
  #绘制散点图
  pdf(file=paste0(i, ".scatter_plot.pdf"), width=7, height=6.5)
  p1=mr_scatter_plot(mrResult, dat)
  print(p1)
  dev.off()
  
  #森林图
  res_single=mr_singlesnp(dat)      #得到每个工具变量对结局的影响
  pdf(file=paste0(i, ".forest.pdf"), width=6.5, height=5)
  p2=mr_forest_plot(res_single)
  print(p2)
  dev.off()
  
  #漏斗图
  pdf(file=paste0(i, ".funnel_plot.pdf"), width=6.5, height=6)
  p3=mr_funnel_plot(singlesnp_results = res_single)
  print(p3)
  dev.off()
  
  #留一法敏感性分析
  pdf(file=paste0(i, ".leaveoneout.pdf"), width=6.5, height=5)
  p4=mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
  print(p4)
  dev.off()
}




#6
#引用包
library(grid)
library(readr)
library(forestploter)

selectMethod=c("MR Egger", "Weighted median", "Inverse variance weighted", "Simple mode", "Weighted mode")     #设置展示的方法
setwd("C:\\biowolf\\proteinMR\\10.forest")            #设置工作目录
files=dir()                           #获取目录下所有文件
files=grep("csv$", files, value=T)    #提取csv结尾的文件

#读取孟德尔随机化分析的结果
data=data.frame()
for(i in files){
  rt=read.csv(i, header=T, sep=",", check.names=F)
  data=rbind(data, rt)
}
data=data[(data$method %in% selectMethod),]
lineVec=cumsum(c(1,table(data[,"exposure"])))

#对数据进行整理
data$' ' <- paste(rep(" ", 10), collapse = " ")
data$'OR(95% CI)'=ifelse(is.na(data$or), "", sprintf("%.3f (%.3f to %.3f)", data$or, data$or_lci95, data$or_uci95))
data$pval = ifelse(data$pval<0.001, "<0.001", sprintf("%.3f", data$pval))
data$exposure = ifelse(is.na(data$exposure), "", data$exposure)
data$nsnp = ifelse(is.na(data$nsnp), "", data$nsnp)
data[duplicated(data$exposure),]$exposure=""

#准备图形参数
tm <- forest_theme(base_size = 18,   #图形整体的大小
                   #可信区间的形状、线条类型、宽度、颜色、两端竖线高度
                   ci_pch = 16, ci_lty = 1, ci_lwd = 1.5, ci_col = "black", ci_Theight = 0.2, 
                   #参考线条的形状、宽度、颜色
                   refline_lty="dashed", refline_lwd=1, refline_col="grey20",
                   #x轴刻度字体的大小
                   xaxis_cex=0.8,
                   #脚注大小、颜色
                   footnote_cex = 0.6, footnote_col = "blue")

#绘制图形
plot <- forestploter::forest(data[, c("exposure","nsnp","method","pval"," ","OR(95% CI)")],
                             est = data$or,
                             lower = data$or_lci95,
                             upper = data$or_uci95,
                             ci_column = 5,     #可信区间所在的列
                             ref_line = 1,      #参考线条的位置
                             xlim = c(0, 2),    #X轴的范围
                             theme = tm,        #图形的参数
)

#修改图形中可信区间的颜色
boxcolor = c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148")
boxcolor = boxcolor[as.numeric(as.factor(data$method))]
for(i in 1:nrow(data)){
  plot <- edit_plot(plot, col=5,row = i, which = "ci", gp = gpar(fill = boxcolor[i],fontsize=25)) # 改col，box的列
}

#设置pvalue的字体
pos_bold_pval = which(as.numeric(gsub('<',"",data$pval))<0.05)
if(length(pos_bold_pval)>0){
  for(i in pos_bold_pval){
    plot <- edit_plot(plot, col=4,row = i, which = "text", gp = gpar(fontface="bold"))  # 改col pvalue的列
  }
}

#在图形中增加线段
plot <- add_border(plot, part = "header", row =1,where = "top",gp = gpar(lwd =2))
plot <- add_border(plot, part = "header", row = lineVec, gp = gpar(lwd =1))
#设置字体大小, 并且将文字居中
plot <- edit_plot(plot, col=1:ncol(data),row = 1:nrow(data), which = "text", gp = gpar(fontsize=12))
plot <- edit_plot(plot, col = 1:ncol(data), which = "text",hjust = unit(0.5, "npc"),part="header",
                  x = unit(0.5, "npc"))
plot <- edit_plot(plot, col = 1:ncol(data), which = "text",hjust = unit(0.5, "npc"),
                  x = unit(0.5, "npc"))

#输出图形
pdf("forest.pdf", width=12, heigh=10)
print(plot)
dev.off()





#7
#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")


#???ð?
library(TwoSampleMR)

exposureID="ieu-a-8"      #??¶????ID
outcomeFile="17441_4_ALPI_Alkaline_phosphatase__intestine.txt"     #?????????ļ?
setwd("C:\\biowolf\\proteinMR\\11.reverse")      #???ù???Ŀ¼

#??ȡ??¶????(??????????)
exposure_dat <- extract_instruments(exposureID, p1=5e-08, p2=5e-08, clump=TRUE)

#??ȡ????????????(????????)
outcome_dat<-read_outcome_data(snps=exposure_dat$SNP,
                               filename=outcomeFile, sep = "\t",
                               snp_col = "rsids",
                               beta_col = "Beta",
                               se_col = "SE",
                               effect_allele_col = "effectAllele",
                               other_allele_col = "otherAllele",
                               pval_col = "Pval",
                               eaf_col = "ImpMAF")


#????¶???ݺͽ??????ݺϲ?
i=gsub(".txt", "", outcomeFile)
outcome_dat$outcome=i
dat <- harmonise_data(exposure_dat, outcome_dat)

#MR-PRESSO?쳣ֵ????(ƫ?е?SNP)
presso=run_mr_presso(dat)
write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file=paste0(i, ".table.MR-PRESSO_Global.csv"))
write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file=paste0(i, ".table.MR-PRESSO_Outlier.csv"))

#?ϵ¶???????????
mrResult=mr(dat)
mrTab=generate_odds_ratios(mrResult)

#?????????ϵ¶????????Ĺ??߱?��
outTab=dat[dat$mr_keep=="TRUE",]
write.csv(outTab, file=paste0(i, ".table.SNP.csv"), row.names=F)

#?????ϵ¶????????????Ľ???
write.csv(mrTab, file=paste0(i, ".table.MRresult.csv"), row.names=F)

#?????Լ???
heterTab=mr_heterogeneity(dat)
write.csv(heterTab, file=paste0(i, ".table.heterogeneity.csv"), row.names=F)

#??Ч?Լ???
pleioTab=mr_pleiotropy_test(dat)
write.csv(pleioTab, file=paste0(i, ".table.pleiotropy.csv"), row.names=F)

#????ɢ??ͼ
pdf(file=paste0(i, ".scatter_plot.pdf"), width=7, height=6.5)
p1=mr_scatter_plot(mrResult, dat)
print(p1)
dev.off()

#ɭ??ͼ
res_single=mr_singlesnp(dat)      #?õ?ÿ?????߱?��?Խ??ֵ?Ӱ??
pdf(file=paste0(i, ".forest.pdf"), width=6.5, height=5)
p2=mr_forest_plot(res_single)
print(p2)
dev.off()

#©??ͼ
pdf(file=paste0(i, ".funnel_plot.pdf"), width=6.5, height=6)
p3=mr_funnel_plot(singlesnp_results = res_single)
print(p3)
dev.off()

#??һ???????Է???
pdf(file=paste0(i, ".leaveoneout.pdf"), width=6.5, height=5)
p4=mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
print(p4)
dev.off()




#8
#???ð?
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)

pvalueFilter=0.05      #pֵ?Ĺ???????
adjPvalFilter=1       #??????pֵ?Ĺ???????

#??????ɫ
colorSel="p.adjust"
if(adjPvalFilter>0.05){
  colorSel="pvalue"
}

setwd("C:\\biowolf\\proteinMR\\12.GO")      #???ù???Ŀ¼
rt=read.csv("IVW.filter.csv", header=T, sep=",", check.names=F)     #??ȡ?????ļ?

#??ȡ??????????, ??????????ת??Ϊ????id
genes=unique(as.vector(rt[,"exposure"]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=entrezIDs[!is.na(entrezIDs)]
gene=as.character(entrezIDs)
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO????????
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$p.adjust<adjPvalFilter),]
#?????????????????Ľ???
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

#??״ͼ
pdf(file="barplot.pdf", width=9, height=7)
bar=barplot(kk, drop=TRUE, showCategory=10, label_format=100, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#????ͼ
pdf(file="bubble.pdf", width=9, height=7)
bub=dotplot(kk, showCategory=10, orderBy="GeneRatio", label_format=100, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()


###########????GOȦͼ###########
ontology.col=c("#00CC33FF", "#FFC20AFF", "#CC33FFFF")
data=GO[order(GO$pvalue),]
datasig=data[data$pvalue<0.05,,drop=F]
BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]
BP = head(BP,6)
CC = head(CC,6)
MF = head(MF,6)
data = rbind(BP,CC,MF)
main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]

#????Ȧͼ????
BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(GO=data$ID,start=1,end=max(BgGene))
rownames(df) = df$GO
bed2 = data.frame(GO=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(GO=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(GO=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5

#????Ȧͼ???岿??
pdf(file="GO.circlize.pdf", width=10, height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.par(track.margin=c(0.01,0.01))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                      #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                    })
circos.clear()
#????Ȧͼ?м???ͼ??
middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)
#????GO??????ͼ??
main.legend = Legend(
  labels = c("Biological Process", "Cellular Component", "Molecular Function"),  type="points",pch=15,
  legend_gp = gpar(col=ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
#???Ƹ?????????pvalue??ͼ??
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()



#9
#???ð?
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05      #pֵ????????
adjPvalFilter=1        #????????pֵ????????

#????ͼ?ε???ɫ
colorSel="p.adjust"
if(adjPvalFilter>0.05){
  colorSel="pvalue"
}

setwd("C:\\biowolf\\proteinMR\\13.KEGG")      #???ù???Ŀ¼
rt=read.csv("IVW.filter.csv", header=T, sep=",", check.names=F)     #??ȡ?????????Ľ????ļ?

#??ȡ??????????????, ??????????ת??Ϊ????id
genes=unique(as.vector(rt[,"exposure"]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=entrezIDs[!is.na(entrezIDs)]
gene=as.character(entrezIDs)
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg????????
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
#kk@result$Description=gsub(" - Homo sapiens \\(human\\)", "", kk@result$Description)

#????idת???ɻ???????
kk = setReadable(kk, #ǰ???????Ľ???
                 OrgDb = "org.Hs.eg.db", #???????ݿ?
                 keyType = "ENTREZID") #Ҫת???Ļ???????
KEGG=as.data.frame(kk)

#?????????????Ľ???
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$p.adjust<adjPvalFilter),]
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#????չʾͨ·????Ŀ
showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

#??״ͼ
pdf(file="barplot.pdf", width=8, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=100, color=colorSel)
dev.off()

#????ͼ
pdf(file="bubble.pdf", width=8, height=7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=100, color=colorSel)
dev.off()
