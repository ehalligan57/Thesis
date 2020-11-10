#Simona Dalin
#11/10/20
#Demonstrating FC vs. significant t-test

library(ggplot2)

######## Capture the folder that the current script is in, and change the working directory to this folder so we can save plots easily ############
currentDirName <- dirname(parent.frame(2)$ofile)
print(currentDirName)
setwd(as.character(currentDirName))

########## First look at an example with large fold change, but non significant t-test ########## 
FCbutNotSig1 <- c(3,5,3,4,20)
FCbutNotSig2 <- c(3,5,3,3,2)

FC <- mean(FCbutNotSig1)/mean(FCbutNotSig2)
Pval <- t.test(FCbutNotSig1, FCbutNotSig2, var.equal = TRUE)

FCbutNotSig.df <- data.frame(name = c(rep("Group 1", 5), rep("Group 2", 5)), 
                             value = c(FCbutNotSig1, FCbutNotSig2))

plot1 <- ggplot(FCbutNotSig.df, aes(x = name, y = value, fill = name)) +
  geom_boxplot() + 
  #geom_point() +
  #geom_jitter(position = position_jitter(0.2)) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1, fill = "black") +
  labs(title = sprintf("FC = %.3s, t-test P-val = %.2e", FC, Pval$p.value),
       x = "Groups",
       y = "Mock Gene Expression") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  theme(legend.position = "none")

ggsave("Large FC, non-sig P-val.pdf", plot1, device = "pdf")


########## Second look at an example with small fold change, but significant t-test ########## 
NotFCbutSig1 <- c(5,4.5,4,5,4)
NotFCbutSig2 <- c(3,4,3,3.5,4)

FC2 <- mean(NotFCbutSig1)/mean(NotFCbutSig2)
Pval2 <- t.test(NotFCbutSig1, NotFCbutSig2, var.equal = TRUE)

NotFCbutSig.df <- data.frame(name = c(rep("Group 1", 5), rep("Group 2", 5)), 
                             value = c(NotFCbutSig1, NotFCbutSig2))

plot2 <- ggplot(NotFCbutSig.df, aes(x = name, y = value, fill = name)) +
  geom_boxplot() + 
  #geom_point() +
  #geom_jitter(position = position_jitter(0.2)) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1, fill = "black") +
  labs(title = sprintf("FC = %.3s, t-test P-val = %.2e", FC2, Pval2$p.value),
       x = "Groups",
       y = "Mock Gene Expression") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  theme(legend.position = "none")

ggsave("Small FC, sig P-val.pdf", plot2, device = "pdf")


########## Third look at an example with large fold change, and significant t-test ########## 
FCandSig1 <- c(15, 13, 12, 15, 20)
FCandSig2 <- c(3,4,3,3.5,4)

FC3 <- mean(FCandSig1)/mean(FCandSig2)
Pval3 <- t.test(FCandSig1, FCandSig2, var.equal = TRUE)

FCandSig3.df <- data.frame(name = c(rep("Group 1", 5), rep("Group 2", 5)), 
                             value = c(FCandSig1, FCandSig2))

plot3 <- ggplot(FCandSig3.df, aes(x = name, y = value, fill = name)) +
  geom_boxplot() + 
  #geom_point() +
  #geom_jitter(position = position_jitter(0.2)) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1, fill = "black") +
  labs(title = sprintf("FC = %.3s, t-test P-val = %.2e", FC3, Pval3$p.value),
       x = "Groups",
       y = "Mock Gene Expression") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  theme(legend.position = "none")

ggsave("Large FC, sig P-val.pdf", plot3, device = "pdf")

