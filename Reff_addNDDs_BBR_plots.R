#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#This code adapted from Arathi's "schisto_halstead_2pops_mda_bouncebackRate_cloudParam.R"

##################################################################################################

require(deSolve)
require(graphics)
require(tidyverse)

params['lamda'] = fin_pars95ci$lamda_twa[which(fin_pars95ci$negLL == min(fin_pars95ci$negLL))]
k<-params['k']
mu_W<-params["mu_W"]
mu_H<-params["mu_W"]

#Vectors to fill for base model (no added NDDs) ##########
  W.seq = c(seq(from=0.01, to=1, by = 0.0025), seq(from=2, to=200, by = 0.05))
  Reff.seq = as.numeric()
  r0.seq = as.numeric()
  dwdt.seq = as.numeric()
  dwdt.seq2 = as.numeric()

#Fill vectors across range of W values    
for(w in 1:length(W.seq)){
    W = W.seq[w]
    Reff.seq[w] = getReff_addNDD(params, W = W, k = k)[1]
    r0.seq[w] = getReff_addNDD(params, W = W, k = k)[2]

    dwdt.seq[w] = (Reff.seq[w]-1)*W*(mu_W+mu_H) 
    dwdt.seq2[w] = (r0.seq[w]-1)*W*(mu_W+mu_H) 
    
}
  
opar<-par()
    
#plot Reff profile ##########
#windows(width = 19, height = 11)    
plot(log(W.seq), Reff.seq, type = 'l', lwd = 3, 
     ylim = c(0,3), xlim = log(c(min(W.seq)+0.003, max(W.seq))),
     xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')  
  abline(h = 1, lty=3, lwd = 2)
  #abline(h = max(r0.seq), lty = 3, col = 2)
  lines(log(W.seq), r0.seq, lty = 2, lwd = 3)
  
    segments(x0 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]),
             x1 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]), #col = 3,
             y0 = 1, y1 = -1, lty = 2, lwd = 2)
    segments(x0 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][8]),
             x1 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][8]), #col = 3,
             y0 = 1, y1 = -1, lty = 2, lwd = 2)
    segments(x0 = log(W.seq[which(Reff.seq == max(Reff.seq))]),
             x1 = log(W.seq[which(Reff.seq == max(Reff.seq))]), #col = 3,
             y0 = max(Reff.seq), y1 = -1, lty = 2, lwd = 2)
    
    axis(1, at = c(log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]), 
                   log(W.seq[which(Reff.seq == max(Reff.seq))]),
                   log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][8])),
         labels = c(expression(paste(italic('W'[bp]), sep = '')),
                    expression(paste(italic(' W'[peak]), sep = '')),
                    expression(paste(italic('    W'[eq]), sep = ''))), #col = 3,
         lty = 2, cex.axis = 1.2)
    axis(1, at = c(-4,0,4), mgp = c(3,1.2,0), cex.axis = 1.2)
    axis(2, at = max(r0.seq), labels = expression(paste(italic('R'[0]), sep = '')), #col = 2, 
         lty=2, las = 2, mgp = c(3,0.75,0), cex.axis = 1.2)
    axis(2, at = c(0,1,2,3), labels = c('0', '1', '2', '3'), cex.axis = 1.2)
    
    mtext(side = 2, text = expression(paste('R'[eff], sep = '')), line = 2.5, cex = 1.4)
    mtext(side = 1, text = expression(paste('     log (', italic('W'), ')', sep = '')),
          line = 3, cex = 1.4)
    
#Plot BBR profile #######
#windows(height = 11, width = 19)
par(mar=c(4.6,6.6,1.1,0.6))    

layout(matrix(c(1,1,1,1,1,1,1,1,1,1,
                1,1,1,1,1,1,1,1,1,1,
                1,1,1,1,1,1,1,1,1,1,
                1,1,1,1,1,1,1,1,1,1,
                2,2,2,2,2,2,2,2,2,2,
                2,3,3,3,3,2,2,2,2,2,
                2,3,3,3,3,2,2,2,2,2,
                2,3,3,3,3,2,2,2,2,2,
                2,2,2,2,2,2,2,2,2,2,
                2,2,2,2,2,2,2,2,2,2), 
                10, 10, byrow = TRUE))
  
plot(log(W.seq), (dwdt.seq*365)/W.seq, type='l', lwd=3, ylim=c(-1, 1.2), xlim =c(-4.3, log(max(W.seq))),
     ylab="", xaxt = 'n', yaxt ='n', xlab = '')
    abline(h = 0, lty=3, lwd = 2)
    lines(log(W.seq), (dwdt.seq2*365)/W.seq, lwd = 3, lty = 2)
    
    segments(x0 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]),
             x1 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]), #col = 3,
             y0 = 0, y1 = -1, lty = 2, lwd = 2)
    segments(x0 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][8]),
             x1 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][8]), #col = 3,
             y0 = 0, y1 = -1, lty = 2, lwd = 2)
    segments(x0 = log(W.seq[which(Reff.seq == max(Reff.seq))]),
             x1 = log(W.seq[which(Reff.seq == max(Reff.seq))]), #col = 3,
             y0 = max((dwdt.seq*365)/W.seq), y1 = -1, lty = 2, lwd = 2)
    
    axis(1, at = c(log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]), 
                   log(W.seq[which(Reff.seq == max(Reff.seq))]),
                   log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][8])),
         labels = c(expression(paste(italic('W'[bp]), sep = '')),
                    expression(paste(italic(' W'[peak]), sep = '')),
                    expression(paste(italic('    W'[eq]), sep = ''))), #col = 3,
         lty = 2,  mgp = c(3,1.2,0), cex.axis = 1.25)
    axis(1, at = c(-4,0,4), mgp = c(3,1.4,0), cex.axis = 1.25)
    axis(2, at = c(-1, -0.5,  -0.25, 0, 0.25,  0.5, 1), cex.axis = 1.25, 
         labels = c('-1.0', '-0.5',  '', '0', '',  '0.5', '1.0'), mgp = c(3,0.75,0))
    
    mtext(side = 2, text = expression(italic(BBR)), line = 2.75, cex = 1.5)
    text(x=5.25, y=0.99, labels='a', pos=1, cex=3)
    
        
#plot dwdt trajectories ############
plot(log(W.seq), dwdt.seq, type = 'l', lwd = 3, ylim = c(-0.002, 0.035), xlim = c(min(log(W.seq)), 4.5),
     xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')
    
    mtext(side = 1, text = expression(paste(' log (', italic('W'), ')', sep = '')), line = 3.5, cex = 1.5)
    mtext(side = 2, text = expression(paste(italic(frac(dW, dt)))), line = 1.75, cex = 1.5)
    
    axis(1, at = c(-4,0,4), mgp = c(3,1,0), cex.axis = 1.25)
    axis(2, at = c(-0.002,0,0.002,0.004,0.006,0.008, 0.01, 0.02, 0.03),
         labels = c('-0.002','0','0.002', '', '0.006','0.008','0.01', '0.02', '0.03'), 
         mgp = c(3,0.75,0), cex.axis = 1.25)
    
    lines(log(W.seq), dwdt.seq2, lwd = 3, lty = 2)
    text(x=4.4, y=0.03, labels='b', pos=1, cex=3)
    legend(x = 1.5, y = 0.0075, legend = c(expression(italic(PDD-free)), expression(italic(PDD))),
           lwd = 2, lty = c(2,1), bty = 'n', cex = 1.25)
    
#zoom to y axis ##########
plot(log(W.seq), dwdt.seq, type = 'l', lwd = 3, xaxt = 'n', yaxt = 'n', bty = 'n', xlab = '',
     ylab = '', #xlab = expression(paste('log (', italic('W'), ')', sep = '')),
     ylim = c(-0.0002, 0.0002), xlim = c(-4,-1))

    #mtext(side = 2, text = expression(paste(italic(frac(dW, dt)))), line = 1.5)
    axis(1, at = c(-4,-3,-1), mgp = c(3,1,0), cex.axis = 1.25)
    axis(2, at = c(-0.0002, 0,0.0002),
         labels = c('-0.0002',  '0', '0.0002'),
         mgp = c(3,0.75,0), cex.axis = 1.25)
    
    
    segments(x0 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]),
             x1 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]), #col = 3,
             y0 = 0, y1 = -1, lty = 2, lwd = 2)
    
    axis(1, at = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]),
         labels = expression(paste(italic('   W'[bp]), sep = '')),
         lty = 2,  mgp = c(3,1.4,0), cex.axis = 1.25)
    
    abline(h = 0, lty=3, lwd = 2)
    lines(log(W.seq), dwdt.seq2, lwd = 2, lty = 2)
    
par(opar)    

