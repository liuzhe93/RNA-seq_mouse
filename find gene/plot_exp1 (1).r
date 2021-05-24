library('DEGseq')

data = read.table("read_counts_rmdup.txt", head=TRUE, sep='\t')
cat(names(data))

if (!file.exists('DEGseq'))
{
DEGexp(geneExpMatrix1=data, geneExpMatrix2=data,
#       expCol1=c(4,5), expCol2=c(2,3), ##case1
#       expCol1=c(6,7), expCol2=c(2,3), ##case2
       expCol1=c(6,7), expCol2=c(4,5), ##case3
       method="MARS", qValue=0.01, foldChange=2,
       thresholdKind=5,
       outputDir='DEGseq'
      )
}

res = read.table("DEGseq/output_score.txt", head=TRUE, sep='\t')
cat(names(res), '\n')
gene = res[,1]
PTB_idx1 = as.character(gene) == "PTBP1"
PTB_idx2 = as.character(gene) == "PTBP2"
PTB_idx3 = as.character(gene) == "PTBP3"

pdf("plot1_MA.pdf")
R = res[,2]
G = res[,3]
y = log2(R/G)
x = 0.5*(log2(R)+log2(G))
sig = res[,9] < 0.01
fld = abs(res[,4]) > log2(2)
plot(x, y, xlab='A', ylab='M', pch='.', cex=0.2, col='grey')
abline(h= log2(2), col='grey')
abline(h=      0 , col='grey')
abline(h=-log2(2), col='grey')
text(max(x[!is.na(x)])-2, log2(2)+0.3, labels='Fold change>2', cex=1.0)
text(max(x[!is.na(x)])-2, -log2(2)-0.3, labels='Fold change<0.5', cex=1.0)
points(x[sig&(!fld)], y[sig&(!fld)], pch='.', cex=0.8, col='black')
points(x[sig&fld], y[sig&fld], pch='+', cex=1.0, col='violet')
points(x[PTB_idx1], y[PTB_idx1], col='red')
text(x[PTB_idx1]+1.4, y[PTB_idx1], labels='PTBP1', cex=1.0, col='red')
points(x[PTB_idx2], y[PTB_idx2], col='red')
text(x[PTB_idx2]+1.4, y[PTB_idx2], labels='PTBP2', cex=1.0, col='red')
points(x[PTB_idx3], y[PTB_idx3], col='red')
text(x[PTB_idx3]+1.4, y[PTB_idx3], labels='PTBP3', cex=1.0, col='red')
dev.off()

pdf("plot2_Cum.pdf")
sel = read.table("gene_selection.txt", head=TRUE, sep='\t')
cat(colnames(sel), '\n')

fold = res[,5]
fold[is.na(fold)] = 0
idx = sort(fold, index.return=TRUE)$ix
fd = fold[idx]
qv = res[,9][idx]

valid = (abs(fd)>0) & (!is.na(qv)) & (qv<=1)
f = fd[valid]*log10(2)
s = sel[idx,][valid,] == 1

plot(c(0,0), c(0,1), pch='.', col='white',
     xlab="Expression fold change (log10)",
     ylab="Cumulative probability")
#lines(f, seq(0,1,length.out=length(f)), col='black')
#sel_names = c('PTB_body', 'PTB_up500', 'PTB_up1000', 'PTB_up2000')
#sel_names = c('nPTB_body', 'nPTB_up500', 'nPTB_up1000', 'nPTB_up2000')
sel_names = c('PTB_body', 'nPTB_body', 'PTB_anti', 'nPTB_anti')
sel_color = c('red', 'brown', 'blue', 'green')
for(i in 1:length(sel_names)) {
    vals = f[s[, sel_names[i]]]
    lines(vals, seq(0,1,length.out=length(vals)), col=sel_color[i])
}
legend("topleft", legend=sel_names, lty=c(1,1,1,1), col=sel_color)

dev.off()

