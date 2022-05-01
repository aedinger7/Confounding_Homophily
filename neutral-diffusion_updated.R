# ACT's attempt at the voter model diffusion problem.
# From https://www.stat.cmu.edu/~cshalizi/homophily-confounding/neutral-diffusion-act.R
# This has been slightly modified from the original code to fix some bugs and run smoothly

# library(ElectroGraph)
# electrograph is not up to date with R, but is only used for network visualization
# functions dependent on electrograph have been commented out for testing

create.network <- function(nodes=100, homoph.factor=1,
                           baseline.density=0.01,
                           same.boost=0.05) {
  
  prop <- rbinom(nodes, 1, 0.5)
  prob.matrix <- 1*(array(prop,c(nodes,nodes)) == t(array(prop,c(nodes,nodes))))*same.boost + baseline.density
  sociomat <- array(rbinom(length(prob.matrix), 1, prob.matrix), dim(prob.matrix))
  sociomat <- 1*(sociomat+t(sociomat)>0)
  diag(sociomat) <- 0
  
  return(list(sociomatrix=sociomat, traits=prop))
}


evolve.network <- function(sociomatrix, choice=rbinom(dim(sociomatrix)[1], 1, 0.5),
                           save.points=1000, iterations.per.save=10) {
  
  nn <- dim(sociomatrix)[1]
  choice.out <- array(NA, c(nn, save.points))
  for (kk in 1:(save.points*iterations.per.save)) {
    pick <- sample(nn, 1)
    nominees = choice[which(sociomatrix[pick,]==1)]
    if (length(nominees)>0)  {
      choice[pick] <- sample(nominees, 1)
    }
    print(pick, choice)
    if (kk %% iterations.per.save == 0) choice.out[,kk/iterations.per.save] <- choice
  }
  return(choice.out)
}



nettest <- create.network(baseline.density=0.002)
print(mean(nettest$sociomatrix))

diag(nettest$socio) <- 0
# plot.obj <- electrograph(nettest$socio)
# initial.positions <- plot(plot.obj, label.cex=0.1, max.t=0.5, edge.color.f=c("gray","blue","red"))

# plot(plot.obj, node.colors=1+nettest$traits, manual.coord=initial.positions$coordinates,
#      label.cex=0.1, max.t=0.5, edge.color.f=c("gray","blue","red"))
# save(nettest, plot.obj, initial.positions, file="dont-lose-this.RData")


initial.choice <- rbinom(dim(nettest$socio)[1], 1, 0.5)
# pdf("diffusion-initial-act.pdf", width=10, height=6)
# plot(plot.obj, node.colors=1+initial.choice, manual.coord=initial.positions$coordinates,
#      label.cex=0.1, max.t=0.5, edge.color.f=c("gray","blue","red"), component.bor=0)
# dev.off()

choice.prog <- evolve.network(nettest$socio, choice=initial.choice)

latest.choice <- choice.prog[,1000]

# pdf("diffusion-midway-act.pdf", width=10, height=6)
# plot(plot.obj, node.colors=1+choice.prog[,297], manual.coord=initial.positions$coordinates,
#      label.cex=0.1, max.t=0.5, edge.color.f=c("gray","blue","red"), component.bor=0)
# dev.off()



#Density of current plot?
dens <- mean(nettest$sociomatrix)
test.zero <- array(rbinom(length(nettest$socio), 1, dens/2), dim(nettest$socio))
test.zero <- test.zero + t(test.zero)

choice.zero <- evolve.network(test.zero, choice=initial.choice)

prog.effect <- zero.effect <- array(NA, c(2,dim(choice.prog)[2]))
for (kk in 1:(dim(choice.prog)[2]/10)) {
  res.e <- glm(choice.prog[,kk] ~ nettest$trait, family=binomial(link=logit))
  prog.effect[,kk] <- summary(res.e)$coef[2,1:2]
  
  res.z <- glm(choice.zero[,kk] ~ nettest$trait, family=binomial(link=logit))
  zero.effect[,kk] <- summary(res.z)$coef[2,1:2]
}

prog.e <- prog.effect[,1:100]
zero.e <- zero.effect[,1:100]
lims <- range(c(prog.e[1,]-2*prog.e[2,], prog.e[1,]+2*prog.e[2,],
                zero.e[1,]-2*zero.e[2,], zero.e[1,]+2*zero.e[2,]))

# pdf("correlation-homoph.pdf", width=10, height=5)
plot(c(0,1000), lims, ty="n", main="Homophilous Network", xlab="Step", ylab="Effect Size")
x.seq <- seq(10, 1000, by=10)
points(x.seq, prog.e[1,], col=2, pch=19)
segments(x.seq, prog.e[1,]-2*prog.e[2,], x.seq, prog.e[1,]+2*prog.e[2,], col=2)
abline(h=0, col=4)
# dev.off()

# pdf("correlation-neutral.pdf", width=10, height=5)
plot(c(0,1000), lims, ty="n", main="Neutral-Model Network", xlab="Step", ylab="Effect Size")
points(x.seq, zero.e[1,], col=1, pch=19)
segments(x.seq, zero.e[1,]-2*zero.e[2,], x.seq, zero.e[1,]+2*zero.e[2,], col=1)
abline(h=0, col=4)
# dev.off()

