args <- commandArgs(trailing=T)
df <- read.table(args[1],h=T)

## df <- head(df[order(df$rab,decreasing=T),],n=10)
## df <- head(df[order(df$rab,decreasing=T),],n=2000)

# head(df,n=200)

## use ids if present
if ("ida" %in% colnames(df)){
    namea = 3
    nameb = 4
} else {
    namea = 1
    nameb = 2
}

df$names <- paste(df[,namea], df[,nameb],sep="-")

par.old <- par()
threshold <- 0.00
######################
## RELATEDNESS PLOT ##
######################
min.rab.threshold <- threshold
df.rab <- df[df$rab>=min.rab.threshold, ]

rab.height=nrow(df.rab)/30
if(rab.height<7)
    rab.height=7

rab.order <- order(df.rab$rab, decreasing=T)
pdf(paste0(args[1], ".relatedness.pdf"), useDingbats=F, height=rab.height,width=2)
par(lty=0,xaxs="i",yaxs="i")
## names.arg=df.rab$names[rab.order],cex.names=0.2
xpos <- barplot(df.rab$rab[rab.order], horiz=T, cex.axis=0.3, las=2, xlim=c(0,1),space = 0.01, main="Relatedness")
## axis(side=1,cex=0.5)
axis(3, cex.axis=0.3, las=2)
mtext(text=df.rab$names[rab.order], line=0.1,side=2,at=xpos, cex=0.2,las=2)
abline(v = c(0.0, 0.05, 0.1, 0.25, 0.5, 1),lwd=0.1,lty='dotted')
dev.off()


##############
## 2of3 IBD ##
##############

min.two.of.three.IBD.threshold <- threshold
df.two.of.three.IBD <- df[df[,"X2of3_IDB"]>=min.two.of.three.IBD.threshold, ]

two.of.three.IBD.height <- nrow(df.two.of.three.IBD)/30
if(two.of.three.IBD.height<7)
    two.of.three.IBD.height=7

two.of.three.IBD.order <- order(df.two.of.three.IBD[,"X2of3_IDB"], decreasing=T)
pdf(paste0(args[1], ".two_of_three_IBD.pdf"), useDingbats=F, height=two.of.three.IBD.height,width=2)
par(lty=0,xaxs="i",yaxs="i")
## names.arg=df.two.of.three.IBD$names[two.of.three.IBD.order],cex.names=0.2
xpos <- barplot(df.two.of.three.IBD[,"X2of3_IDB"][two.of.three.IBD.order], horiz=T, cex.axis=0.3, las=2, xlim=c(0,1),space = 0.01, main="2 of 3 alleles IBD")
axis(3, cex.axis=0.3, las=2)
mtext(text=df.two.of.three.IBD$names[two.of.three.IBD.order], line=0,side=2,at=xpos, cex=0.2,las=2)
abline(v = c(0.0, 0.05, 0.1, 0.25, 0.5, 1),lwd=0.1,lty='dotted')
dev.off()


################
## INBREEDING ##
################
df.inbreeding <- data.frame(name=as.character(c(as.character(df[,namea]), as.character(df[,nameb]))), F=c(df$Fa, df$Fb))
inbreeding.height <- length(unique(df.inbreeding$name))/10
if (inbreeding.height<7)
    inbreeding.height <- 7
l <- list()
for (name in unique(df.inbreeding$name)){
    l[[name]] <- df.inbreeding[df.inbreeding$name == name,"F"]
}
l <- l[order(sapply(l, function(x){median(x)}),decreasing=T)]
pdf(paste0(args[1], ".inbreeding.pdf"), useDingbats=F, height=inbreeding.height, width=2)
boxplot(l,horizontal=T, las=2, space=0.01, cex.axis=0.3, cex=0.1, ylim=c(0,1),space=0.01, boxlwd=1, main="Inbreeding\nCoefficient")
abline(v = c(0.0, 0.05, 0.1, 0.25, 0.5, 1),lwd=0.1,lty='dotted')
dev.off()


################
## R0 R1 KING ##
################
## df = df[df$Fa<0.1 & df$Fb < 0.1,]
pdf(paste0(args[1], ".R0_R1_KING.pdf"), useDingbats=F) ## , height=inbreeding.height, width=2)
par(mfrow=c(2,2))
plot(df$R1, df$R0, pch=19, cex=0.1, xlab="R1", ylab="R0")
plot(df$R1, df$KING, pch=19,cex=0.1, xlab="R1", ylab="KING")
plot(df$R1, df$R0, xlim=c(0,1), ylim=c(-0.2,0.5), pch=19, cex=0.1, xlab="R1", ylab="R0")
plot(df$R1, df$KING, xlim=c(0,1), ylim=c(-0.2, 0.3), pch=19,cex=0.1, xlab="R1", ylab="KING")
dev.off()
