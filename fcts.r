## All functions

getModelComponents <- function(b) { 
  # A function that calculates all components needed to calculate the bias
  # correction as in Greven & Kneib (2010)
  #
  # Args: 
  #   b        = gamObject
  #   analytic = FALSE if the numeric hessian of the (restricted) marginal log-
  #              likelihood from the lmer optimization procedure should be used.
  #              Otherwise (default) TRUE, i.e. use a analytical version that 
  #              has to be computed.
  #
  # Returns:
  #   model = List of components needed to calculate the bias correction
  #   
  
  model             <- list()
  model$df          <- sum(b$edf)
  ppq               <- length(coef(b))
  s                 <- length(b$smooth)
  S0                <- S <- matrix(0, ppq, ppq)
  model$sp          <- b$sp
  model$db.drdy     <- list()
  model$s           <- s
  model$B           <- matrix(0, nrow = s, ncol = s)  
  
  ## Write everything into a return list
  model$X           <- model.matrix(b)   
  if(b$method == "REML") {
    model$p         <- b$min.edf
    #model$B         <- b$outer.info$hess[1:s, 1:s]
	#model$data		<- b$data
	model$dlogS 	<- obtainHessianREML(b)
  }
  if(b$method == "ML") {
    model$p         <- 0
    #model$B         <- b$outer.info$hess[1:s, 1:s]
	#model$data		<- b$data
	model$dlogS 	<- obtainHessianREML(b)
  }

  model$n           <- length(b$y)        
  model$S           <- S
  model$sig2        <- b$sig2
  model$residuals   <- b$residuals
  model$method      <- b$method
  model$Vb          <- b$Vp/b$sig2
  model$db.dy       <- model$Vb %*% t(model$X)
  model$db.dr       <- b$db.drho
  model$b           <- coef(b)
  model$y           <- as.vector(b$y)
  model$D           <- b$deviance
  
  if (model$method %in% c("REML", "ML")) {
  
    model$db.dr       <- b$db.drho
    
    for(k in 1:s) {
      
      model$Slist[[k]] <- S0
      
      if(length(b$smooth[[k]]$S) != 1) {
        stop("Can't handle more then one penalty matrix per smooth")
      }
      
      model$Slist[[k]][b$smooth[[k]]$first.para:b$smooth[[k]]$last.para, 
                      b$smooth[[k]]$first.para:b$smooth[[k]]$last.para] <- 
                      b$smooth[[k]]$S[[1]] #* b$smooth[[k]]$S.scale
            
      model$sp[k]  <- model$sp[k]#/b$smooth[[k]]$S.scale
      
      model$S <- model$S + model$Slist[[k]] * model$sp[k]
      
      model$db.drdy[[k]] <- - model$sp[k] * model$Vb %*% 
                              (model$Slist[[k]] %*% model$db.dy)
                              
  
      model$D <- as.numeric(b$deviance + coef(b) %*% model$S %*% coef(b))                              
      
    }
  } else {
  
    model$db.dr       <- matrix(0, nrow = length(coef(b)), ncol = s)
    
    for(k in 1:s) {
      
      model$Slist[[k]] <- S0
      
      if(length(b$smooth[[k]]$S) != 1) {
        stop("Can't handle more then one penalty matrix per smooth")
      }
      
      model$Slist[[k]][b$smooth[[k]]$first.para:b$smooth[[k]]$last.para, 
                      b$smooth[[k]]$first.para:b$smooth[[k]]$last.para] <- 
                      b$smooth[[k]]$S[[1]] #* b$smooth[[k]]$S.scale
      
      model$sp[k]  <- model$sp[k]#/b$smooth[[k]]$S.scale
      
      model$S <- model$S + model$Slist[[k]] * model$sp[k]
      
      model$db.drdy[[k]] <- - model$sp[k] * model$Vb %*% 
                              (model$Slist[[k]] %*% model$db.dy)
                              
      model$db.dr[,k] <- - model$sp[k] * model$Vb %*% model$Slist[[k]] %*%
                         model$b
      
    }
  }

  return(model)
}

calculateGaussianBc <- function(model) {
  # A function that calculates the analytic representation of the bias 
  # corrections in linear mixed models, see Greven & Kneib (2010).
  #
  # Args: 
  #   model    = From getAllModelComponents()
  #
  # Returns:
  #   df = Bias correction (i.e. degrees of freedom) for a linear mixed model.
  #           0
  G <- matrix(0, nrow = model$s, ncol = model$n) 
  
  if (model$method %in% c("REML","ML")) {
    for (k in 1:model$s) {  
      G[k,] <- mREMLdrdy(k, model)      
			#if (is.null(model$data)) {
			#stop("For REM/ML data must be provided, i.e. add control = list(keepData = TRUE) to gam call")
			#}
        for(j in k:model$s) {   
          model$B[k, j] <- model$B[j, k] <- mREMLdrdr(k , j, model)        
        }
    }
  } else {
    if (model$method == "UBRE") {    
      for (k in 1:model$s) {  
        G[k,] <- UBREdrdy(k, model)    

        for(j in k:model$s) {   
          model$B[k, j] <- model$B[j, k] <- UBREdrdr(k, j, model)        
        }
      }
    } else {
      for (k in 1:model$s) {  
        G[k,] <- GCVdrdy(k, model)    

        for(j in k:model$s) {   
          model$B[k, j] <- model$B[j, k] <- GCVdrdr(k, j, model)        
        }
      }   
    }
  }  
  
  #Rchol   <- chol(model$B)
  #L1      <- backsolve(Rchol, G, transpose = TRUE)
  #Lambday <- - backsolve(Rchol, L1)
  
  Lambday <- solve(- model$B, G)

  df <- numeric(model$s)
  
  for (j in 1:model$s) {
    df[j] <- - model$sp[j] * Lambday[j,] %*% (t(model$db.dy) %*%(model$Slist[[j]] %*% (model$db.dy %*% model$y)))
  }
  
  return(df)
}

mREMLdrdy <- function(k, model) {
 
  dD.drdy <- 2 * (as.vector(model$db.dr[,k]) %*% t(model$X) %*% 
             (model$X %*% model$db.dy - diag(model$n)) +
             t(model$X %*% model$b - model$y) %*% model$X %*% model$db.drdy[[k]] + 
             model$b %*% model$S %*% model$db.drdy[[k]] + 
             as.vector(model$db.dr[,k]) %*% model$S %*% model$db.dy + 
             model$sp[k] * model$b %*% model$Slist[[k]] %*% model$db.dy)
             
  dD.dr   <- 2 * t(model$X %*% model$b - model$y) %*% 
             model$X %*% as.vector(model$db.dr[,k]) + 
             2 * as.vector(model$db.dr[,k]) %*% model$S %*% model$b + 
             model$sp[k] * model$b %*% model$Slist[[k]] %*% model$b
             
  dD.dy   <- 2 * t(model$X %*% model$b - model$y) %*% 
             (model$X %*% model$db.dy - diag(model$n)) + 
             2 * model$b %*% model$S %*% model$db.dy
  
  return((model$n - model$p)/(2 * model$D) * 
         (dD.drdy - as.numeric(dD.dr) /model$D * dD.dy))

}

mREMLdrdr <- function(k, j,  model){

  

#    dK.drdr  <- (-model$sp[k] * model$sp[j] * sum(diag(model$Vb %*% model$Slist[[k]] %*% model$Vb %*% model$Slist[[j]])) - ncol(model$S)/model$sp[k]^2)/2

  if (k == j) {
  
    dlogXtXS.drdr <- model$sp[k] * sum(diag(model$Vb %*% model$Slist[[k]])) - model$sp[k] * model$sp[j] * sum(diag(model$Vb %*% model$Slist[[k]] %*% model$Vb %*% model$Slist[[j]]))
  	
	#dK.drdr  <- dK.drdr + (model$dlogS$det1 + model$sp[k] * sum(diag(model$Vb %*% model$Slist[[k]])))/2 #dim det1 check for multiple sp
	
	#dK.drdr  <- (dlogXtXS.drdr - model$dlogS$det2[k, k])/2 
	
	dK.drdr  <- (dlogXtXS.drdr)/2 
    
	db.drdr <- as.vector(model$db.dr[,j]) - 2 * model$sp[k] * model$Vb %*% (model$Slist[[k]] %*% as.vector(model$db.dr[,j]))
  
  } else {
  
    dlogXtXS.drdr <- - model$sp[k] * model$sp[j] * sum(diag(model$Vb %*% model$Slist[[k]] %*% model$Vb %*% model$Slist[[j]]))
  
    #dK.drdr  <- (dlogXtXS.drdr - model$dlogS$det2[k,j])/2
  
    dK.drdr  <- (dlogXtXS.drdr)/2
  
    db.drdr <- -  model$Vb %*% (model$sp[k] * model$Slist[[k]] %*% as.vector(model$db.dr[,j]) + model$sp[j] * model$Slist[[j]] %*% as.vector(model$db.dr[,k]))
  
  }

  dD.drk   <- 2 * t(model$X %*% model$b - model$y) %*% 
             model$X %*% as.vector(model$db.dr[,k]) + 
             2 * as.vector(model$db.dr[,k]) %*% model$S %*% model$b + 
             model$sp[k] * model$b %*% model$Slist[[k]] %*% model$b
			 
  dD.drj   <- 2 * t(model$X %*% model$b - model$y) %*% 
             model$X %*% as.vector(model$db.dr[,j]) + 
             2 * as.vector(model$db.dr[,j]) %*% model$S %*% model$b + 
             model$sp[j] * model$b %*% model$Slist[[j]] %*% model$b		

#  dD.drkdrj <- 2 * (crossprod(model$X %*% as.vector(model$db.dr[,k])) + model$residuals %*% model$X %*% as.vector(db.drdr) +
#				as.vector(db.drdr) %*% model$S %*% model$b + as.vector(model$db.dr[,k]) %*% model$S %*% as.vector(model$db.dr[,j]) + 
#				model$sp[k] * model$b %*% model$Slist[[k]] %*% as.vector(model$db.dr[,j]) + 
#				model$sp[j] * model$b %*% model$Slist[[j]] %*% as.vector(model$db.dr[,k]))

  dD.drkdrj <- 2 * (as.vector(model$db.dr[,j]) %*% t(model$X) %*% model$X %*% as.vector(model$db.dr[,k]) + t(model$X %*% model$b - model$y) %*% model$X %*% as.vector(db.drdr) +
				as.vector(db.drdr) %*% model$S %*% model$b + as.vector(model$db.dr[,k]) %*% model$S %*% as.vector(model$db.dr[,j]) + 
				model$sp[k] * model$b %*% model$Slist[[k]] %*% as.vector(model$db.dr[,j]) + 
				model$sp[j] * model$b %*% model$Slist[[j]] %*% as.vector(model$db.dr[,k]))

  if (k == j) {
  
	dD.drkdrj <- dD.drkdrj + model$sp[k] * model$b %*% model$Slist[[k]] %*% model$b
  
  }  

  return((model$n - model$p)/(2 * model$D) * (dD.drkdrj - dD.drk * dD.drj /model$D) + dK.drdr)#dim problem, negative?
}

UBREdrdy <- function(k, model) {
 
  return(2 * (as.vector(model$db.dr[,k]) %*% t(model$X) %*% 
             (model$X %*% model$db.dy - diag(model$n)) -
             model$residuals %*% model$X %*% model$db.drdy[[k]]))

}

UBREdrdr <- function(k, j, model) {
 
  if(k == j) {
  
    dt.drdr <- - model$sp[j]  * 
              sum(diag(t(model$db.dy) %*% model$Slist[[k]] %*% model$db.dy)) +
              2 * model$sp[j]^2 * 
              sum(diag(t(model$db.dy) %*% model$Slist[[k]] %*% model$Vb %*%
              model$Slist[[k]] %*% model$db.dy))
              
    db.drdr <- model$db.dr[,k] - 
               2 * model$sp[k] * (model$Vb %*% model$Slist[[k]] %*% model$db.dr[,k])
               
  } else {
  
    dt.drdr <- 2 * (model$sp[j] * model$sp[k] * 
              sum(diag(t(model$db.dy) %*% model$Slist[[k]] %*% model$Vb %*%
              model$Slist[[j]] %*% model$db.dy)))

    db.drdr <- - model$Vb %*% (model$sp[k] * model$Slist[[k]] %*% model$db.dr[,j] +
                               model$sp[j] * model$Slist[[j]] %*% model$db.dr[,k])
  }
  
  return(2 * (as.vector(model$db.dr[,k]) %*% t(model$X) %*% 
             model$X %*% as.vector(model$db.dr[,j]) -
             model$residuals %*% model$X %*% db.drdr) + 2 * model$sig2 * dt.drdr) 
}

GCVdrdy <- function(k, model) {
 
  return(
    2 * model$n / (model$n - model$df)^2 * (
    model$db.dr[,k] %*% t(model$X) %*% (model$X %*% model$db.dy - diag(model$n))
    - model$residuals %*% model$X %*% model$db.drdy[[k]]) - 
    4 * model$n / (model$n - model$df)^3 * model$sp[k] * 
    sum(diag(t(model$db.dy) %*% model$Slist[[k]] %*% model$db.dy))  *
    model$residuals %*% (diag(model$n) - model$X %*% model$db.dy)
  )
  
}

GCVdrdr <- function(k, j, model) {
 
  if(k == j) {
  
    dt.drdr <- - model$sp[j]  * 
              sum(diag(t(model$db.dy) %*% model$Slist[[k]] %*% model$db.dy)) +
              2 * model$sp[j]^2 * 
              sum(diag(t(model$db.dy) %*% model$Slist[[k]] %*% model$Vb %*%
              model$Slist[[k]] %*% model$db.dy))
              
    db.drdr <- model$db.dr[,k] - 
               2 * model$sp[k] * (model$Vb %*% model$Slist[[k]] %*% model$db.dr[,k])
               
  } else {
  
    dt.drdr <- 2 * (model$sp[j] * model$sp[k] * 
              sum(diag(t(model$db.dy) %*% model$Slist[[k]] %*% model$Vb %*%
              model$Slist[[j]] %*% model$db.dy)))

    db.drdr <- - model$Vb %*% (model$sp[k] * model$Slist[[k]] %*% model$db.dr[,j] +
                               model$sp[j] * model$Slist[[j]] %*% model$db.dr[,k])
  }
  
  dD.drdr <- 2 * (model$db.dr[,k] %*% t(model$X) %*% model$X %*% model$db.dr[,j]
             - model$residuals %*% model$X %*% db.drdr)
             
  dD.drj   <- - 2 * model$residuals %*% model$X %*% model$db.dr[,j]
  
  dD.drk   <- - 2 * model$residuals %*% model$X %*% model$db.dr[,k]
  
  dt.drj  <- - model$sp[j] * 
               sum(diag(t(model$db.dy) %*% model$Slist[[j]] %*% model$db.dy))
  
  dt.drk  <- - model$sp[k] * 
               sum(diag(t(model$db.dy) %*% model$Slist[[k]] %*% model$db.dy))
  
  dGCV.drdr <- 2 * model$n / (model$n - model$df)^3 *
               ((model$n - model$df) / 2 * dD.drdr + dD.drk * dt.drj
               + dD.drj * dt.drk + 3 / (model$n - model$df) * model$D * dt.drj * dt.drk
               + model$D * dt.drdr)
               
  return(dGCV.drdr)
}

cAIC <- function(object) {
  # A function that calls the bias correction functions.
  #
  # Args: 
  #   object = Object of class gam.
  #
  # Returns:
  #   list   = The list contains the conditional log-likelihood; the estimated 
  #            conditional prediction error (degrees of freedom); If a new model
  #            was fitted, the new model and an boolean indicating weather a new
  #            model was fitted; the conditional Akaike information, caic.
  #

  if (!inherits(object, c("gam"))) {
    stop("Class needs to be gam")
  }
  
  if (family(object)$family != "gaussian") {
    stop("Method only applicable for Gaussian responses")
  }

  if (object$boundary == TRUE) {
    stop("Parameters are estimated on the boundary. These need to be excluded")
  }

  if (!(object$method %in% c("REML", "ML", "UBRE", "GCV"))) {
   stop("Method not supplied")
  }  
  
  model <- getModelComponents(object)
                    
  df <- calculateGaussianBc(model)

  return(df)
}

obtainHessianREML <- function(b) {
	
#  G <-mgcv:::gam.setup(formula.gam(b), b$pterms, data = b$data)
#  
#  G$rS <- mgcv:::mini.roots(G$S, G$off, ncol(G$X), G$rank)
#  
#  Ssp <- mgcv:::totalPenaltySpace(G$S, G$H, G$off, ncol(G$X))
#  
##  Ssp <- totalPenaltySpace(G$S, G$H, G$off, ncol(G$X))
#  G$Eb <- Ssp$E
#  G$U1 <- cbind(Ssp$Y, Ssp$Z)
#  G$Mp <- ncol(Ssp$Z)
#  G$UrS <- list()
#  if (length(G$S) > 0) 
#      for (i in 1:length(G$S)) G$UrS[[i]] <- t(Ssp$Y) %*% G$rS[[i]]
#  else i <- 0
#  if (!is.null(G$H)) {
#      G$UrS[[i + 1]] <- t(Ssp$Y) %*% mroot(G$H)
#  }
  
#  asd <- gam.reparam(G$UrS, log(b$sp), 2)
  
  asd <- 0*b$fitted.values[1] 
  
  return(asd)
}
	