
#Simulate the asymmetric issue.
asymmetry.sim <- function (num.nodes=400, scale=3, offset=0,
                           y.noise=0.02, friend.noms=1,
                           response="linear", nominate.by="distance",
                           time.trend=0.3) {
  #num.nodes=400; scale=3; offset=0; y.noise=0.02; response="not-linear"; nominate.by="not-distance";friend.noms=1; time.trend=0.3
  
  #covariates.
  invlogit <- function(cc) exp(cc)/(1+exp(cc))
  xx <- runif(num.nodes)
  distances <- as.matrix(dist(xx, diag=TRUE, upper=TRUE))
  
  #potential friends.
  adj.probs <- array(rbinom(length(distances), 1, invlogit(offset-scale*distances)), dim(distances))
  diag(adj.probs) <- 0
  adj.probs[lower.tri(adj.probs)] <- t(adj.probs)[lower.tri(adj.probs)]
  
  fixed.x.dist <- distances
  diag(fixed.x.dist) <- Inf
  nominees <- t(rbind(sapply(1:num.nodes, FUN=function(ii) {
    possibles <- which(adj.probs[ii,] == 1)
    return(ifelse(rep(nominate.by=="distance",friend.noms),
                  possibles[which(order(fixed.x.dist[ii,possibles])<=friend.noms)],
                  sample(possibles, size=friend.noms, prob=invlogit(-abs((xx[possibles] - 0.5))))
    ))
  })))
  #num.nodes-by-friend.nom should be the size of the object.
  
  #true adjacency matrix.
  aa.mat <- array(0, dim(distances))
  for (ii in 1:num.nodes) aa.mat[ii, nominees[ii,]] <- 1
  #backwards matrix.
  rev.aa.mat <- t(aa.mat)
  #reciprocated matrix.
  #recip.aa.mat <- aa.mat*rev.aa.mat
  
  y1 = (xx-0.5)^3+rnorm(num.nodes,0,y.noise)
  y2 = y1+rnorm(num.nodes,time.trend*xx,y.noise)
  
  infl.y1 <- aa.mat%*%y1
  back.y1 <- rev.aa.mat%*%y1
  #recip.y1 <- recip.aa.mat%*%y1
  
  XX <- cbind(1, y1, infl.y1, back.y1)
  
  trial <- lm(y2 ~ y1 + infl.y1 + back.y1)
  #summary(trial)
  coef.table <- summary(trial)$coefficients[,1]
  v.plus.c <- c(diag(summary(trial)$cov.unscaled), summary(trial)$cov.unscaled[3,4])*summary(trial)$sigma^2
  #print(summary(trial)$cov.unscaled)
  z.stat <- (coef.table[3]-coef.table[4])/sqrt(v.plus.c[3]+v.plus.c[4]+2*v.plus.c[5])
  
  out <- cbind(c(coef.table, v.plus.c, z.stat))
  rownames(out) <- c("int.b", "auto.b", "infl.b", "back.b",
                     "int.s", "auto.s", "infl.s", "back.s",
                     "cov.infl.back",
                     "z.stat")
  
  return(out)
  
}

result <- replicate(5000, asymmetry.sim(friend=1, time.trend=0.3))

dim(result)

pdf("timeseriesmodel-act-5000.pdf", width=12, height=6)

par(mfrow=c(1,2))
hist(result[3,], main="Effect of Phantom `Influencer' on `Influenced' in Time Series", xlab="Regression Coefficient"); hist(result[10,], main="z-score of Directional Difference", xlab=paste("Proportion greater than 0:", mean(result[10,]>0)))
dev.off()


hist(result[2,])
