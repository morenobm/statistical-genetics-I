((3/4)^3)*(1/4)*4
crossprob = dbinom(0:4, size = 4, prob = 1/4); crossprob
#Soma acumulada
cumuprob = cumsum(crossprob); cumuprob

library(ggplot2)
ngreen = seq(0,4,1)
aux = data.frame(ngreen, crossprob, cumuprob)
ggplot(data = aux, aes(x = factor(ngreen))) + 
  geom_point(aes(y = crossprob, shape = "PMF"), size = 2.5, alpha = .8, color = "#10342d") + 
  geom_point(aes(y = cumuprob, shape = "CDF"), size = 2.5, alpha = .8, color = "red") +
  theme_bw() + theme(legend.position = 'top', legend.title = element_blank(),
                     text = element_text(size = 18)) +
  labs(x = "Number of green seeds", y = "Probability of occurrence") +
  scale_y_continuous(breaks = seq(0,1,0.1))

#esperança com dados de dado

dado <- c(1,2,3,4,5,6)
prob <- c(1/6,1/6,1/6,1/6,1/6,1/6)
expec = sum(dado * prob)
expec

# a esperança de um dado é 3,5

#rbinom - distribuição binomial para simulações

#var no R é sempre a amostral

set.seed(8)
sam = rnorm(n = 1000, mean = 1000, sd = 100)
sam[1:10]

rho = .75
s21 = 1.5
s22 = 4
cov12 = rho*sqrt(s21*s22)
Sigma = matrix(c(
  s21, cov12, 
  cov12, s22
), nrow = 2)

set.seed(8970011)
dat = MASS::mvrnorm(n = 100000, mu = c(5.5, 25), Sigma = Sigma)
colnames(dat) = c("hgw", "npods")

E_x = sum(ngreen * crossprob)
V_x = sum(ngreen^2*crossprob) - expec^2
n = 5
iter = 10000
zn = NULL
for (i in 1:iter){
  y = sum(rbinom(n = n, size = 4, prob = 1/4))
  E_y = n * E_x
  V_y = n * V_x
  zn[i] = (y - E_y) / sqrt(V_y)
}

ggplot() +
  geom_histogram(aes(x = zn, y = after_stat(density)), color = 'black') + 
  stat_function(fun = function(x) dnorm(x, mean = 0, sd = 1), xlim = c(-3,3),
                geom = "polygon", alpha = .5, color = 'black') + 
  theme_bw() + 
  theme(text = element_text(size = 18)) + 
  labs(x = "Standardized variable", y = "Density", title = paste("Sample size: ", n))


## exercício

 
#Using your USP number as seed, simulate phenotypic values (don’t forget, P=G+E) using the normal distribution (rnorm). 
#Then, compute the mean and the variance of P, G and E, and build a histogram, a density plot, a box plot and a scatter 
#plot depicting the relationship between P #and G . Do it for three different sample sizes: 10, 100 and 1000.

#Fiz tudo errado!! Tem que fazer uma simulação para cada um, uma para G e outra para E e soma os dois. Essa parte agora 
#está certo.

?rnorm
#Amostra com 10 indivíduos
set.seed(8970011)
gen_10 <- rnorm(n=10, mean = 2000, sd = 276)
env_10 <- rnorm(n=10, mean = 859, sd = 32)

pheno_10 <- gen_10 + env_10

hist_10 <- ggplot(data = data.frame(pheno_10)) + 
  geom_histogram(aes(x = pheno_10), bins = 20, color = 'white', fill = "#10342d") + 
  labs(x = "Pheno 10", y = "Count") + 
  theme_minimal() + theme(text = element_text(size = 18))
hist_10

density_10 <- ggplot(data = data.frame(pheno_10)) + 
  geom_density(aes(x = pheno_10),color = 'white', fill = '#10342d', alpha = .75) + 
  labs(x = "Pheno 10", y = "Density") + 
  theme_minimal() + theme(text = element_text(size = 18))
density_10

box_10 <- ggplot(data = data.frame(pheno_100)) +
  geom_boxplot(aes(x = sam), fill = "#b2bc63") + 
  annotate(geom = 'point', x = mean(sam), y = 0, shape = 17, color = "#10342d", size = 3)+
  theme_minimal() + theme(axis.text.y = element_blank(), text = element_text(size = 18)) +
  
  labs(x = "Seed yield (kg/ha)")

#Amostra com 100 indivíduos
set.seed(8970011)
gen_100 <- rnorm(n=100, mean = 2000, sd = 276)
env_100 <- rnorm(n=100, mean = 859, sd = 32)

pheno_100 <- gen_100 + env_100

dados_100 <- data.frame(pheno_100, gen_100, env_100)

mean (gen_100)
var)gen_100

hist_100 <- ggplot(data = data.frame(pheno_100)) + 
  geom_histogram(aes(x = pheno_100), bins = 20, color = 'white', fill = "#10342d") + 
  labs(x = "Pheno 100", y = "Count") + 
  theme_minimal() + theme(text = element_text(size = 18))
hist_100

density_100 <- ggplot(data = data.frame(pheno_100)) + 
  geom_density(aes(x = pheno_100),color = 'white', fill = '#10342d', alpha = .75) + 
  labs(x = "Pheno 100", y = "Density") + 
  theme_minimal() + theme(text = element_text(size = 18))
density_100

box_100 <- ggplot(data = data.frame(pheno_100)) +
  geom_boxplot(aes(x = sam), fill = "#b2bc63") + 
  annotate(geom = 'point', x = mean(sam), y = 0, shape = 17, color = "#100342d", size = 3)+
  theme_minimal() + theme(axis.text.y = element_blank(), text = element_text(size = 18)) +
  
  labs(x = "Seed yield (kg/ha)")

#Exercício de exploração de dados
dat = read.csv("https://raw.githubusercontent.com/mauricioaraujj/Pan_African_Trials_Network/refs/heads/main/data/data.csv", 
               sep = ';')
dim(dat)

dat_sub = subset(dat, env == "E007")
dim(dat_sub)