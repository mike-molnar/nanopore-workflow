library(ggplot2)
library(data.table)
args = commandArgs(trailingOnly=TRUE)
dat1 <- fread(args[1], select = c("type", "predicted_length"))

dat2 <- subset(dat1, (type=='INS' | type=='DEL') & predicted_length >= as.integer(args[3]) & predicted_length <= as.integer(args[4]))

p <- ggplot(dat2, aes(x=predicted_length, color=type, fill=type)) + 
geom_histogram(position="dodge2", binwidth=as.integer(args[5])) +
labs(title=args[6], x="Indel Length", y="Total Indels") +
theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank())
ggsave(args[2])
