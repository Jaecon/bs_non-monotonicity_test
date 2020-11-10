# Just monotonicity

rm(list=ls())

library(lfe)
library(dplyr)
library(parallel)
library(boot)

mydata <- read.csv("db1b.csv")
mydata <- mydata %>% mutate(rc = as.factor(rc))
mydata <- mydata %>% mutate(tc = as.factor(tc))
#mydata <- mydata %>% mutate(phase2_tr = as.factor(phase2_tr))

df <- mydata[mydata$oligo_route_r==1,] %>% na.omit()
df_ph2 <- df[df$phase2_tr==1, ]
df_ph1 <- df[df$phase2_tr==0, ]

# Number of bootstrapped sample
B = 5000

status.mon = 0
status.nonmon = 0


# Objective function for the statistic

cores = detectCores()-1

obj <- function(theta, data, type){
	# data will be bs.ph2
  lower_index = round(theta[1])
  upper_index = round(theta[2])
  sorted = sort(data[, "pred_response"])
  lb = sorted[lower_index]
  ub = sorted[upper_index]
  bdd_bs_dt <- rbind(data[  (data$pred_response>=lb & data$pred_response<=ub),  ], df_ph1) %>% as.data.frame()
  
  fit <- lfe::felm(log_avgfare_trc ~ phase2_tr + phase2_tr:pred_response + distance:avg_fuel_price + number_concarr_tr + number_dircarr_tr 
                   | rc + tc, data=bdd_bs_dt)
  
  bdd_x <- bdd_bs_dt[bdd_bs_dt$phase2_tr==1, "pred_response"]
  Q <- sqrt( sum( (bdd_x-mean(bdd_x))^2 ) )
  beta <- unname( coef(fit)[names(coef(fit))=="phase2_tr:pred_response"] )
  # monotonicity test like EE
  if(type == "maxPlus"){
    return(Q*beta)
  }
  else if(type == "maxMinus"){
    return(-Q*beta)
  }
  else
    stop("type mismatched")
}

t.fun <- function(ph2, test_type){
	# Note that ui %*% theta - ci >= 0.
    #            lb ub  
    ui = rbind(c(-1, 1),
               c( 1, 0),
               c( 0,-1) 
               )
    # Length of phase2==1 observations for the initial range
    M = length(ph2[, "pred_response"])
    # Smoothing parameter for minimum observations in the range
    m = round(1.5*log(M)^2)
  	# Constants on constraints
    ci = c(m-1, 0.99, -(M+0.01))
    
    res_maxPlus <- constrOptim(c(1,M), function(x) obj(theta=x, data=ph2, type="maxPlus"), 
                               NULL, ui=ui, ci=ci, control = list(fnscale = -1))
    res_maxMinus <- constrOptim(c(1,M), function(x) obj(theta=x, data=ph2, type="maxMinus"), 
                               NULL, ui=ui, ci=ci, control = list(fnscale = -1))
    if(test_type == "nonmon"){
      return(
        min(res_maxPlus$value, res_maxMinus$value)
        )
    }
    else if(test_type == "mon"){
      return(res_maxPlus$value)
    }
    else{
      stop("test_type mismatched")
    }
}

bs_x <- function(ph2, test_type){
    bs.r.index <- sample(unique(ph2$r), replace=TRUE)
    bs.ph2 <- do.call(
    	rbind, lapply(bs.r.index, function(x) ph2[ph2$r==x,])
    	)
    t.fun(bs.ph2, test_type)
  }

RNGkind("L'Ecuyer-CMRG")
set.seed(3656)
t.mon <- mclapply(1:B, 
					function(i){
						return(
							bs_x(df_ph2, test_type = "mon")
							)						
					},
					mc.cores = cores
					) %>% unlist()

save(list=ls(), file="01-mon.Rdata")

RNGkind("L'Ecuyer-CMRG")
set.seed(7911)
t.nonmon <- mclapply(1:B, 
					function(i){
						return(
							bs_x(df_ph2, test_type = "nonmon")
							)
					},
					mc.cores = cores
					) %>% unlist()


save(list=ls(), file="02-nonmon.Rdata")

#---------------------------#
# compute original theta, t #
#---------------------------#
t0.mon <- t.fun(df_ph2, "mon")
t0.nonmon <-t.fun(df_ph2, "nonmon")

#----------------#
# simple perc CI #
#----------------#
mon.perc.ci <- quantile(t.mon, 0.95)
nonmon.perc.ci <- quantile(t.nonmon, c(0.025,0.975))


#---------#
# BCa CI  #
#---------#
z0.fun <- function(t, t0){
	top <- sum(t<t0)
	qnorm(top/B)
}

jack <- mclapply(unique(df_ph2$r), 
                 function(i){ 
                   jack.ph2 <- df_ph2[df_ph2$r!=i, ]
                   c(t.fun(jack.ph2, "mon"), t.fun(jack.ph2, "nonmon"))
                 }
                 , mc.cores = cores 
) %>% do.call(rbind.data.frame, args = .)

colnames(jack) <- c("mon", "nonmon")


bca.fun <- function(t, t0, type, alpha=0.05){
  if(type == "nonmon"){
  	alpha <- alpha/2
  }
  z0 = z0.fun(t, t0)
  L = mean(jack[,type])-jack[,type]
  
  atop = sum(L^3)
  abot = 6*(sum(L^2))^(3/2)
  ahat = atop/abot
  
  alpha1 = pnorm(
  	z0 + (z0+qnorm(alpha))/(1-ahat*(z0+qnorm(alpha)))
  	)
  alpha2 = pnorm(
  	z0 + (z0+qnorm(1-alpha))/(1-ahat*(z0+qnorm(1-alpha)))
  	)
  if(type == "nonmon"){
  	confpoint = quantile(t, probs=c(alpha1,alpha2))
  }
  else if(type == "mon"){
  	confpoint = quantile(t, probs=alpha2)
  }
  else{}
  list(theta=t0, CI=confpoint, z0=z0, acc=ahat, alpha=alpha)
}

bca.mon <- bca.fun(t.mon, t0.mon, "mon")
bca.nonmon <- bca.fun(t.nonmon, t0.nonmon, "nonmon")

save(list=ls(), file="03-bca.Rdata")

load("03-bca.Rdata")

mycolor <- RColorBrewer::brewer.pal(n=3,"Dark2")

windows()
par(mfrow=c(1,2))
hist(t.nonmon, breaks = 100, xlim = c(-10,10), freq=F,
     main = "Histogram of Estimates from Bootstrap Samples",
     xlab = expression(tilde("T")[bold(HH)])
     )
abline(v = t0.nonmon, col=mycolor[1], lty=2, lwd=2)
abline(v = bca.nonmon$CI, col=mycolor[2], lty=3, lwd=2)
#abline(v = nonmon.perc.ci)
legend("topright", legend = c(expression(hat("T")[bold(HH)]),"BCa CI"), col = mycolor[1:2], lty=2:3, cex=0.65)
qqnorm(t.nonmon, frame=F, pch=1)
qqline(t.nonmon, col=mycolor[3], lwd=2)
