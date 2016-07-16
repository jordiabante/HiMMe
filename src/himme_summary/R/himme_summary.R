# Arguments
args <- commandArgs(trailingOnly = TRUE)
x<-read.table(args[1],sep="\t")
kmer_size<-as.numeric(args[2])
outfile<-args[3]

# Model
model_mean=-1.350
model_sd=0.03824

# Names df
colnames(x)<-c("id","len","s","shat")

# Kmer size used
out=c("K-mer size used",kmer_size)

# Number of entries
n_entries=nrow(x)
out=rbind(out,c("Number of entries",n_entries))

# Median length
x_median_l=median(x[,2])
out=rbind(out,c("Median length",x_median_l))

# Length variance
x_variance=var(x[,2])
out=rbind(out,c("Length variance",x_variance))

# Mean score
x_mean_s=mean(x[,3])
out=rbind(out,c("Mean score",x_mean_s))

# Score Variance
x_variance=var(x[,3])
out=rbind(out,c("Score variance",x_variance))

# 95% CI corrected score
x_density=density(x[,3])
x_quantiles=quantile(x_density$x,c(0.05,0.95))
out=rbind(out,c("Score interval [5%,95%]",paste(x_quantiles[1],",",x_quantiles[2])))

# Mean score hat
x_mean_shat=mean(x[,4])
out=rbind(out,c("Mean corrected score",x_mean_shat))

# Score hat Variance
x_variance=var(x[,4])
out=rbind(out,c("Score hat variance",x_variance))

# 95% CI corrected score
x_density=density(x[,4])
x_quantiles=quantile(x_density$x,c(0.05,0.95))
out=rbind(out,c("Corrected score interval [5%,95%]",paste(x_quantiles[1],",",x_quantiles[2])))

# HiMMe coeefficient
zscores=(x$shat-model_mean)/model_sd
himme_coeff=mean(zscores)/sd(zscores)
out=rbind(out,c("HiMMe coefficient",himme_coeff))

# Save output
write.table(out,file=outfile,sep='\t',quote=F,row.names=F,col.names=F)
