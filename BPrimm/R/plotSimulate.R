`plotSimulate` <-
function(chr, pos, b.true, main, b, ylab="Coefficient",
 b.cut=0.05, lty=1, ylim=c(-1,1), mar=c(3,4,3,1), chr.lab=1){
  
  if(!is.character(chr) && !is.numeric(chr)){
    stop("chr must be a character/numerical vector\n")
  }
  
  chr1 = chr
  if(is.character(chr)){
    chr1[chr=="X" | chr=="x"] = "99"
    chr1[chr=="Y" | chr=="y"] = "100"
    chr1 = as.numeric(chr1)
  }

  nchr = length(unique(chr))
  
  if(any(!(chr1 %in% 1:100))){
    stop("unrecognized chromsome\n")
  }

  if(!is.numeric(pos)){
    stop("pos must be a numeric vector\n")
  }
  
  cL = 1.01*tapply(pos, chr1, max)
  ppos = rep(NA, length(pos))
  chrL = c(0, cumsum(cL))
  
  uchr = unique(chr1)
  for(i in 1:length(uchr)){
    whichi = which(chr1==uchr[i])
    ppos[whichi] = pos[whichi] + chrL[i]
  }
  
  if(is.null(ylim)){
    ylim = range(c(b.true, b))
  }
  par(mar=mar, xaxt="n")
  plot(ppos, b.true, type="n", xlab="", ylab=ylab, col="darkblue", 
    ylim=ylim, main = main)
  abline(h=0, col="gray")
  ats=0.5*(chrL[-1]+chrL[-length(chrL)])

  if(chr.lab==1){
    mtext(1:nchr, side=1, line=0.5, at=ats)
  }else if(chr.lab==2){
    mtext(seq(1,nchr,by=2), side=1, line=0.5, at=ats[seq(1,nchr,by=2)])
  }
  abline(v=chrL, lty=2, col="gray")
  
  ww0 = which(abs(b) > 1e-10)
  if(length(ww0) > 0)
    segments(ppos[ww0], 0, ppos[ww0], b[ww0], lty=lty, col="red")
  wl = which(abs(b) > b.cut)
  if(length(wl) > 0)
    points(ppos[wl], b[wl], pch=22, col="red", bg="red", cex=0.5)
  
  wl = which(b.true !=0)
  segments(ppos[wl], 0, ppos[wl], b.true[wl], lty=1, col="darkblue")
  points(ppos[wl], b.true[wl], pch=20, col="darkblue")
}
