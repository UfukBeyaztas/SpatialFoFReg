predict_sff_qr <- function(object, xnew, Wnew)
{
  
  b0hat <- object$surfaces$b0hat
  bhat <- object$surfaces$bhat
  rhohat <- object$surfaces$rhohat
  gpy <- object$surfaces$grid_x
  gpx <- object$surfaces$grid_y
  
  temp1 <- (xnew %*% bhat) * (gpx[2] - gpx[1])
  for(i in 1:dim(temp1)[1])
    temp1[i,] <- temp1[i,] + b0hat
  
  fn <- function(f) as.matrix(Wnew %*% f %*% rhohat * (gpy[2] - gpy[1]))
  yhat0 <- temp1
  yhat1 <- temp1 + fn(yhat0)
  yhat.diff <- max(abs(yhat0 - yhat1))
  Lfy <- 1
  while(yhat.diff > 0.001 & Lfy < 1000){
    yhat1 <- temp1 + fn(yhat0)
    yhat.diff <- max(abs(yhat0 - yhat1))
    yhat0 <- yhat1
    Lfy <- Lfy + 1
  }
  yhat_pred <- yhat1
  
  return(yhat_pred)
}
