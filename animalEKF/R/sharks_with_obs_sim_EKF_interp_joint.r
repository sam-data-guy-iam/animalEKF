sharks_with_obs_sim_EKF_interp_joint <- function(env_obj) {

	for (s in env_obj$sharks_with_obs) {

		#print("tis values")
		#print(s)
		#print(Xpart_history[1:i,c("time_in_state","lambda"),1,s])
		#print(shark_intervals[[ s ]])
		
		#sequence of regular steps for which we will update values
		#should be different for each shark
		#sequence of step since last observation
		env_obj$steps_to_resamp[[ s ]] <- min((max(env_obj$steps_to_resamp[[ s ]]) + 1), env_obj$i-1):(env_obj$i-1)
		
		
			
		for (k in 1:env_obj$nstates) {
			
			
		  env_obj$sigma_draw[,k,s] <- pmax(1e-15, MCMCpack::rinvgamma(n=env_obj$npart, shape=env_obj$sigma_pars[,2*k-1,s], scale=env_obj$sigma_pars[,2*k,s]))
			env_obj$tau_draw[,k,s] <- pmax(1e-15, MCMCpack::rinvgamma(n=env_obj$npart, shape=env_obj$tau_pars[,2*k-1,s],  scale=env_obj$tau_pars[,2*k,s]))
			
		}
		
		
		for (p in 1:env_obj$npart) { 
	  
			#particle covariance, one for each state. only second part depends on state though
			
			env_obj$Qt[1:2, 1:2,, p,s] <- MCMCpack::riwish(v=env_obj$Particle_errvar[[ s ]][[ p ]]$dof, S=env_obj$Particle_errvar[[ s ]][[ p ]]$sig)
			
							
			#draw values for mu_alpha, mu_beta
			#here beta is the mean of the log-transformed angle, which we have to change to 
			
			for (k in 1:env_obj$nstates) {
			
				#this is the error for xt| x_{t-1}: first part is the same for each state
				diag(env_obj$Qt[3:4,3:4,k,p, s]) <- c(env_obj$sigma_draw[p,k,s], env_obj$tau_draw[p,k,s])
		   
				#draw block covariance matrices, same as D before, depends on the state of xt ->yt
				#this is R_k
				env_obj$XY_errvar_draw[,,k,p,s] <- MCMCpack::riwish(v=env_obj$XY_errvar[[ s ]][[ p ]][[ k ]]$dof, S=env_obj$XY_errvar[[ s ]][[ p ]][[ k ]]$sig)
				
				#print(XY_errvar_draw[,,k,p,s])
				#mu has diagonal matrix
				#on first step dont have any turn, assume you know it

				
				#draw logV and normalized turn for each state from current (or prior estimates)
			
				#logv_angle_mu_draw[p,,k] <- rnorm(n=2, mean=mu[[ s ]][[ p ]][[ k ]], sd=sqrt(diag(Qt[3:4, 3:4,k,p, s])*mu[[ s ]][[ p ]]$V[[ k ]]))
				env_obj$logv_angle_mu_draw[p,,k, s] <- rnorm(n=2, mean=env_obj$mu[k,,"mu",p,s], sd=sqrt(diag(env_obj$Qt[3:4, 3:4,k,p, s]) * env_obj$mu[k,,"V",p,s]))
				#logv_angle_draw[p,,k] <-    rnorm(n=2, mean=logv_angle_mu_draw[p,,k], sd=sqrt(diag(Qt[3:4, 3:4,k,p, s])))
			
			}
		}#loop over part and k	
		#print(XY_errvar_draw[,,,1:5,s])
		
		env_obj$logv_angle_mu_draw[,"turn",,s] <- normalize_angle(env_obj$logv_angle_mu_draw[,"turn",,s])		
		
		env_obj$cov_err_hist["X","particle",, env_obj$i,s,"orig"] <- env_obj$cov_err_hist["X","particle",, env_obj$i,s,"resamp"] <-  env_obj$Qt["X","X",1,,s]
		env_obj$cov_err_hist["Y","particle",, env_obj$i,s,"orig"] <- env_obj$cov_err_hist["Y","particle",, env_obj$i,s,"resamp"] <- env_obj$Qt["Y","Y",1,,s]
		env_obj$cov_err_hist["cov","particle",, env_obj$i,s,"orig"] <- env_obj$cov_err_hist["cov","particle",, env_obj$i,s,"resamp"] <- env_obj$Qt["X","Y",1,,s]

		env_obj$cov_err_hist["X",-1,, env_obj$i,s,"orig"] <- env_obj$cov_err_hist["X",-1,, env_obj$i,s,"resamp"] <- env_obj$XY_errvar_draw["X","X",,,s]
		#print(apply(cov_err_hist["X",-1,, env_obj$i,s,"orig"], 1, summary))
		env_obj$cov_err_hist["Y",-1,, env_obj$i,s,"orig"] <- env_obj$cov_err_hist["Y",-1,, env_obj$i,s,"resamp"] <- env_obj$XY_errvar_draw["Y","Y",,,s]
		#print(apply(cov_err_hist["Y",-1,, env_obj$i,s,"orig"], 1, summary))
		env_obj$cov_err_hist["cov",-1,, env_obj$i,s,"orig"] <- env_obj$cov_err_hist["cov",-1,, env_obj$i,s,"resamp"] <- env_obj$XY_errvar_draw["X","Y",,,s]	
		#print(apply(cov_err_hist["cov",-1,, env_obj$i,s,"orig"], 1, summary))
		
		#multiply by gradient since later will be variance of theta
		
		#print(	logv_angle_mu_draw[,"turn",])			
		#Qt[4,4,,,s] <- Qt[4,4,,,s]*t(grad_ginv(psi=logv_angle_mu_draw[,"turn",])^2)
		#logv_angle_draw[,"turn",] <- normalize_angle(log2rad(logv_angle_draw[,"turn",]))
			
		#re-loop
		#print(round(apply(logv_angle_draw[,"turn",], 2, summary), digits=3))
		
		#mk_prev is starting X,Y coordinate and logv, bearing from the starting point of the interval
		#from that, take fractions of dt for each observation

		for (p in 1:env_obj$npart) {
			for (k in 1:env_obj$nstates) {
			
				#mk actual is x_{t-1}, where do we go from here with each lambda_t, to each observation?
				env_obj$mk_prev[,k,p,s] <- env_obj$f(mk=env_obj$mk_actual[,p,s], new_logv=env_obj$logv_angle_mu_draw[p,"logv",k, s], 
													 theta=env_obj$logv_angle_mu_draw[p,"turn",k, s], dtprev=env_obj$reg_dt) #a_{t+1}
				
				
				#mk_prev[4,k,p,s] <- normalize_angle(mk_prev[4,k,p,s])
				Fx_tmp <- env_obj$Fx(mk=env_obj$mk_actual[,p,s], dtprev=env_obj$reg_dt)
				env_obj$Pk_prev[,,k,p,s] <- as.matrix(Matrix::nearPD(Fx_tmp %*% env_obj$Pk_actual[,,p,s] %*% t(Fx_tmp) + env_obj$Qt[,,k,p, s], ensureSymmetry=TRUE)$mat) #R_{t+1}
								
			}
			
			#generate X-Y coordinates for starting at time t given actual movement (mk_actual) at time t-1

			#Xpart_history[ i, c("X","Y"),p, s] <- reject_sampling(mu=mk_prev[c("X","Y"),1,p,s], cmat=Qt[1:2,1:2,1,p], prev_val=mk_actual[,p,s])$val	
			
			#do this instead of just using the mean	
			
			#Xpart[,c("X","Y"),k,"next_t",s]  <- t(apply(Xpart[,,k,"curr",s], 1, function(x) h(mk=x, dtprev=reg_dt)))
			#difference in x-y coordinates
			#Xpart[,c("X","Y"),k,"d",s] <- Xpart[,c("X","Y"),k,"next_t",s] - Xpart[,c("X","Y"),k,"curr",s]
			
		}
		
		
		
		
		#add a little fudge so always have interpolation, never occur at the boundary
		env_obj$j_list[[ s ]][[ env_obj$i ]] <- pmin(pmax((env_obj$ynext[rownames(env_obj$ynext)==s,"date_as_sec"] - env_obj$t_reg[env_obj$i])/env_obj$reg_dt, 1e-5), 1-(1e-5))
				
		print(paste("j:", paste(round(env_obj$j_list[[ s ]][[ env_obj$i ]], digits=4), collapse=", ")))
		
		env_obj$MuY[[ s ]] <- array(NA, dim=c(2, env_obj$yobs_sharks[ s ], env_obj$nstates, env_obj$npart), dimnames=list(c("X","Y"), 1:env_obj$yobs_sharks[ s ], env_obj$state_names, env_obj$pnames)) 
		#rep(list(rep(list(rep(list(matrix(NA, ncol=1, nrow=2)), env_obj$yobs_sharks[ s ])), env_obj$nstates)), env_obj$npart)
		env_obj$SigY[[ s ]] <- array(NA, dim=c(2,2, env_obj$yobs_sharks[ s ], env_obj$nstates, env_obj$npart), dimnames=list(c("X","Y"), c("X","Y"), 1:env_obj$yobs_sharks[ s ], env_obj$state_names, env_obj$pnames)) 
		#rep(list(rep(list(rep(list(diag(2)), env_obj$yobs_sharks[ s ])), env_obj$nstates)), env_obj$npart)
		#Pk_prev_interp[[ s ]] <- array(NA, dim=c(4,4, env_obj$yobs_sharks[ s ], env_obj$nstates, env_obj$npart), dimnames=list(1:4,1:4, 1:env_obj$yobs_sharks[ s ], state_names, pnames)) 
		#rep(list(rep(list(rep(list(diag(4)), env_obj$yobs_sharks[ s ])), env_obj$nstates)), env_obj$npart)	
		env_obj$Kgain[[ s ]] <- array(NA, dim=c(4,2, env_obj$yobs_sharks[ s ], env_obj$nstates, env_obj$npart), dimnames=list(1:4,c("X","Y"), 1:env_obj$yobs_sharks[ s ], env_obj$state_names, env_obj$pnames)) 
		#rep(list(rep(list(rep(list(matrix(0, ncol=2, nrow=4)), env_obj$yobs_sharks[ s ])), env_obj$nstates)), env_obj$npart)	

		
		
		
		#print(env_obj$yobs_sharks)
		#prediction of y is the direct interpolation
		#j_tmp i the differences, not the original fractions
		j_tmp <- diff(c(0, env_obj$j_list[[ s ]][[ env_obj$i ]]))
		
		for (p in 1:env_obj$npart) {
			for (k in 1:env_obj$nstates) {
				for (y in 1:env_obj$yobs_sharks[ s ]) {
					
					if (y==1) {
					
						#Fx_tmp <- Fx(mk=mk_prev[,k,p,s], dtprev=j_tmp[ y ]*env_obj$reg_dt)
						env_obj$MuY[[ s ]][,y,k,p] <- keep_finite(env_obj$h(mk=env_obj$mk_prev[,k,p,s], dtprev=j_tmp[ y ] * env_obj$reg_dt))
						Hx_tmp <- env_obj$Hx(mk=env_obj$mk_prev[,k,p,s], dtprev=j_tmp[ y ] * env_obj$reg_dt)
						
					}
					else {
					
						#take previous x-y vaulues and starting logv and bearing
					
						mk_tmp <- c(env_obj$MuY[[ s ]][,y-1,k,p], env_obj$mk_prev[3:4,k,p,s])
						#Fx_tmp <- Fx(mk=mk_tmp, dtprev=j_tmp[ y ]*env_obj$reg_dt)
						env_obj$MuY[[ s ]][,y,k,p] <- keep_finite(env_obj$h(mk=mk_tmp, dtprev=j_tmp[ y ] * env_obj$reg_dt))
						Hx_tmp <- env_obj$Hx(mk=mk_tmp, dtprev=j_tmp[ y ]*env_obj$reg_dt)
					
					}
					#print(paste("k=",k))
					#print(Hx_tmp)
					
					#Pk_prev_interp[[ s ]][,,y,k,p] <- as.matrix(Matrix::nearPD(Fx_tmp%*%Pk_actual[,,p,s]%*%t(Fx_tmp) + (j_tmp[ y ])*Qt[,,k,p], ensureSymmetry=TRUE)$mat) #R_{t+1}
								
					#multiply the result by j_tmp multiplies the matrix by j_tmp^2, which is what we want			
					env_obj$SigY[[ s ]][,,y,k,p] <- keep_finite(as.matrix(Matrix::nearPD(Hx_tmp %*% env_obj$Pk_prev[,,k,p,s] %*% t(Hx_tmp)  + (j_tmp[ y ])* env_obj$XY_errvar_draw[,,k,p,s], ensureSymmetry=TRUE)$mat))
				
					
					# Fx_tmp <- Fx(mk=mk_prev[,k,p,s], dtprev=j_list[[ s ]][[ i ]][ y ]*env_obj$reg_dt)
					# Pk_prev_interp[[ s ]][,,y,k,p] <- as.matrix(Matrix::nearPD(Fx_tmp%*%Pk_actual[,,p,s]%*%t(Fx_tmp) + (j_list[[ s ]][[ i ]][ y ]^2)*Qt[,,k,p], ensureSymmetry=TRUE)$mat) #R_{t+1}
								
					# MuY[[ s ]][,y,k,p] <- keep_finite(h(mk=Xpart[p,,k,"curr",s], dtprev=j_list[[ s ]][[ i ]][ y ]*env_obj$reg_dt))
					# Hx_tmp <- Hx(mk=Xpart[p,,k,"curr",s], dtprev=j_list[[ s ]][[ i ]][ y ]*env_obj$reg_dt)
					# SigY[[ s ]][,,y,k,p] <- keep_finite(as.matrix(Matrix::nearPD(Hx_tmp%*%Pk_prev_interp[[ s ]][,,y,k,p]%*%t(Hx_tmp)  + (j_list[[ s ]][[ i ]][ y ]^2)*XY_errvar_draw[,,k,p,s], ensureSymmetry=TRUE)$mat))
				
				}
			}
		}

	}
	
	invisible(NULL)

}
	