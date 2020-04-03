sharks_first_obs_EKF_1d_interp_joint <- function(env_obj) {
	
	
	for (s in env_obj$sharks_first_obs) {
	
		ids <- (env_obj$d[,"t_intervals"] == env_obj$i) & (env_obj$tags == s) 
					
		env_obj$y_first <- env_obj$d[ ids, c("X","logvelocity","date_as_sec","t_intervals","state.guess2","shark_obs_index"), drop=FALSE]
		env_obj$j_list[[ s ]][[ env_obj$i ]] <- pmin(pmax((env_obj$y_first[,"date_as_sec"] - env_obj$t_reg[env_obj$i])/env_obj$reg_dt, 1e-5), 1-(1e-5))
		jtmp <- env_obj$j_list[[ s ]][[ env_obj$i ]]	
		
		print(paste("j:", paste(round(env_obj$j_list[[ s ]][[ env_obj$i ]], digits=4), collapse=", ")))

		#choose what the behavior at beginning of time i is, based on what was observed before.
		#ignore this, really, dont store anything for it
		print(env_obj$y_first)
		
	
		
		prev_z <- low_var_sample(wts=table(factor(env_obj$y_first[,"state.guess2"], levels=1:env_obj$nstates)), M=env_obj$npart)
						
		env_obj$Xpart_history[env_obj$i,"lambda",,s]  <- env_obj$lambda_matrix[,env_obj$i,s]  <- prev_z
		#Xpart_history[i,"region",,s] <- y_first["region"]
		
		#initial bearing approximation		
		nobs <- nrow(env_obj$y_first)

		first_pos <- env_obj$y_first[1, "X"]
		anchor_loc <- env_obj$y_first[nobs, "X"]

		if (nobs >= 3) {
			anchor_dist <- first_pos - anchor_loc #purposely reverse, so we are going in the opposite direction, to the beginning
			anchor_time_diff <- diff(env_obj$y_first[c(1, nobs), "date_as_sec"])
			anchor_vel <- anchor_dist/anchor_time_diff
		}

		anchor_time_back <- jtmp[nobs] * env_obj$reg_dt
		

		#speed and turn (unused) angles given start XY and given state
		sigma_draw <- MCMCpack::rinvgamma(n=env_obj$npart, env_obj$sigma_pars[,,s][ cbind(1:env_obj$npart, 2*prev_z -1)], env_obj$sigma_pars[,,s][ cbind(1:env_obj$npart, 2*prev_z)])
		
		env_obj$Qt[2,2,,,s][ cbind(prev_z, 1:env_obj$npart) ] <- sigma_draw
		
		for (p in 1:env_obj$npart) { 
			
			#particle covariance, one for each state. only second part depends on state though
			
			env_obj$Pk_actual[,,p,s]  <- env_obj$Pk_prev[,,prev_z[ p ],p,s] 	  
			#print(Qt[1,1,p])
			env_obj$Qt[1,1,prev_z[ p ],p,s]  <- MCMCpack::riwish(v=env_obj$Particle_errvar[[ s ]][[ p ]]$dof, S=env_obj$Particle_errvar[[ s ]][[ p ]]$sig)
			#print(Qt[1,1,p])	
			#print(Particle_errvar[[ s ]][[ p ]])
			#draw values for mu_alpha, mu_beta
			#here beta is the mean of the log-transformed angle, which we have to change to 
								
			#draw block covariance matrices, same as D before, depends on the state of xt ->yt
				
			env_obj$logv_angle_mu_draw[p,"logv",prev_z[ p ], s] <- as.numeric(mvtnorm::rmvnorm(n=1, mean=env_obj$mu[prev_z[p], "mu", p, s], 
																							sigma=as.matrix(env_obj$Qt[2, 2, prev_z[ p ],p,s] * env_obj$mu[prev_z[p], "V", p, s])))
			
		
			#take the logvelocity mu back fraction of sectonds to beginning of interval, then use that log-velocity to go forwards next.
			
			mk_tmp <- c(anchor_loc, -1 * env_obj$logv_angle_mu_draw[p,"logv",prev_z[ p ], s])

			env_obj$mk_actual[,p,s]  <- env_obj$f(mk=mk_tmp, new_logv= env_obj$logv_angle_mu_draw[p,"logv",prev_z[ p ], s], dtprev=anchor_time_back) #a_{t+1}
			#print(mk_actual[,p,s])
			Fx_tmp <- env_obj$Fx(mk=mk_tmp, dtprev=anchor_time_back)
			#print(Fx_tmp)
			#print(Fx_tmp%*%Pk_actual[,,p,s]%*%t(Fx_tmp) + Qt[,,p])
			
			env_obj$Pk_actual[,,p,s]  <- as.matrix(Matrix::nearPD(Fx_tmp %*% env_obj$Pk_actual[,,p,s]  %*% t(Fx_tmp) + jtmp[nobs] * env_obj$Qt[,,prev_z[ p ],p,s] , ensureSymmetry=TRUE)$mat) #R_{t+1}
			

			#Xpart_history[ i, "logv",p,s ] <- logv_angle_draw[p,"logv",prev_z[ p ]]
			#simulate location if traveled backwards (negative velocity)
			#Xpart_history[ i, "X", p, s ] <- rnorm(n=1, mean=h(mk=c(y_first[ "X" ], -1*Xpart_history[i, "logv", p, s]), dtprev=dt_tmp), sd=sqrt(Qt[1,1,p]))
			
			#go back distance from y_first
			#mk_actual[2,p,s] <- -1*mk_actual[2,p,s]

			env_obj$Xpart_history[ env_obj$i, c("X","logv"), p, s]  <- mvtnorm::rmvnorm(n=1, mean=env_obj$mk_actual[,p,s] , sigma=env_obj$Pk_actual[,,p,s] )
			
			
		}#loop over part
		

		env_obj$Xpart_history[ env_obj$i, "time_in_state",,s ] <- 1
		
		# correct this in the history, since are actually going the reverse direction
		# simulated X is actually the simulated starting point going back a fraction from the first observed location
		
		
		env_obj$error_beforesamp_quantiles[env_obj$i,,s]  <-  quantile(apply(env_obj$Xpart_history[env_obj$i, c("X","logv"),,s], 2, function(x) sum(abs(env_obj$h(mk=x, dtprev=jtmp) - env_obj$y_first[,"X"]))), p=c(.1, .5, .9))
	}
	invisible(NULL)

}		