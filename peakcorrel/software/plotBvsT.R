library(ggplot2);
B <- read.table("B_27ac.bed",sep="\t");
T <- read.table("T_27ac.bed",sep="\t");
B_length <- B$V3-B$V2;
T_length <- T$V3-T$V2;
df <- data.frame(length =
                   c(B_length,T_length),
                 cell=c(rep_len('B',length(B_length)),rep_len('T',length(T_length))));
pdf('CAST-27ac-density.pdf',width = 6, height = 4,bg = "transparent");
ggplot(df, aes(x=length)) + geom_density(aes(group=cell, color=cell, fill=cell), alpha=0.3) + 
  scale_fill_manual(values=c("#228b22", "#0000ff","#56B4E9")) +
  scale_colour_manual(values=c("#228b22", "#0000ff","#56B4E9")) + 
  scale_x_log10() + 
  theme(axis.text = element_text(colour="black"));
dev.off();
t.test(df$length ~ df$cell);