sharks_first_obs_EKF_interp_joint <- function(env_obj) {

	for (s in env_obj$sharks_first_obs) {
	

		ids <- (env_obj$d[,"t_intervals"] == env_obj$i) & (env_obj$tags ==s) 
				
		env_obj$y_first <- env_obj$d[ ids, c("X","Y","logvelocity","bearing.to.east.tonext.rad","date_as_sec","t_intervals","state.guess2","shark_obs_index", "region"), drop=FALSE]
		
		env_obj$j_list[[ s ]][[ env_obj$i ]] <- pmin(pmax((env_obj$y_first[,"date_as_sec"] - env_obj$t_reg[env_obj$i])/env_obj$reg_dt, 1e-5), 1-(1e-5))
		jtmp <- env_obj$j_list[[ s ]][[ env_obj$i ]]				

		
		print(paste("j:", paste(round(env_obj$j_list[[ s ]][[ env_obj$i ]], digits=4), collapse=", ")))
		
		#choose what the behavior at beginning of time i is, based on what was observed before.
		#ignore this, really, dont store anything for it
		print(env_obj$y_first)
	

		num_reject <- c()
		#approximate a draw of the behavior of the first one	
		prev_z <- low_var_sample(wts=table(factor(env_obj$y_first[,"state.guess2"], levels=1:env_obj$nstates)), M=env_obj$npart)
		
		env_obj$Xpart_history[env_obj$i,"lambda",,s] <- env_obj$lambda_matrix[,env_obj$i,s] <- prev_z
		
		#initial bearing approximation		
		nobs <- nrow(env_obj$y_first)
		first_pos <- env_obj$y_first[1, c("X","Y")]
		initial_bearing_forward <- env_obj$y_first[1, "bearing.to.east.tonext.rad"]
		initial_time_back <- jtmp[1] * env_obj$reg_dt
		
		# in case start time and first observed time don't match up
		time0_estimate <- env_obj$h(mk=c(env_obj$y_first[, c("X", "Y", "logvelocity")], normalize_angle(initial_bearing_forward - pi)), dtprev=jtmp[1] * env_obj$reg_dt)
		anchor_pos <- env_obj$y_first[nobs, c("X", "Y")]
	
		#use interpolation estimate location at end of first interval. 
		#go forward a given number of steps
		
		end_of_interval_estimate <- env_obj$h(mk=env_obj$y_first[nobs, c("X", "Y", "logvelocity", "bearing.to.east.tonext.rad")], dtprev=(1-jtmp[nobs]) * env_obj$reg_dt)
		dxy <- apply(rbind(end_of_interval_estimate, time0_estimate), 2, diff)
		
		bearing_est_to_first <- normalize_angle(atan2(x=dxy[1], y=dxy[2]))
		
				
				
		#speed and turn (unused) angles given start XY and given state
		env_obj$sigma_draw[,,s][ cbind(1:env_obj$npart, prev_z) ] <- MCMCpack::rinvgamma(n=env_obj$npart, env_obj$sigma_pars[,,s][ cbind(1:env_obj$npart, 2*prev_z -1)], env_obj$sigma_pars[,,s][ cbind(1:env_obj$npart, 2*prev_z)])
		
		#tau_draw   <- MCMCpack::rinvgamma(n=env_obj$npart, env_obj$tau_pars[,,s][ cbind(1:env_obj$npart, 2*prev_z -1)], env_obj$tau_pars[,,s][ cbind(1:env_obj$npart, 2*prev_z)])
		env_obj$tau_draw[,prev_z,s] <- rep(1e-1, env_obj$npart)	
		
		#this will become env_obj$Pk actual
		env_obj$Pk_prev[4,4,,,s][ cbind(prev_z, 1:env_obj$npart) ] <- 1e-3
		
		#take overall bearing or just the one observed
		
				
		for (p in 1:env_obj$npart) { 
			
			#particle covariance, one for each state. only second part depends on state though
			env_obj$Pk_actual[,,p,s] <- env_obj$Pk_prev[,,prev_z[ p ],p,s]	  

			env_obj$Qt[1:2, 1:2,,p,s] <- MCMCpack::riwish(v=env_obj$Particle_errvar[[ s ]][[ p ]]$dof, S=env_obj$Particle_errvar[[ s ]][[ p ]]$sig)
			env_obj$cov_err_hist[,"particle",p, env_obj$i,s,"orig"] <- env_obj$Qt[, ,1,p,s][ c(1,4,2) ]
			env_obj$cov_err_hist[,"particle",p, env_obj$i,s,"resamp"] <- env_obj$Qt[, ,1,p,s][ c(1,4,2) ]							
			
			#draw values for mu_alpha, mu_beta
			#here beta is the mean of the log-transformed angle, which we have to change to 
			
			#this is the error for xt| x_{t-1}: first part is the same for each state
			
			env_obj$Qt[3:4, 3:4,,p,s] <- diag(c(env_obj$sigma_draw[p, prev_z[ p ],s], env_obj$tau_draw[p, prev_z[ p ],s]))
			
			#draw block covariance matrices, same as D before, depends on the state of xt ->yt
				
			env_obj$logv_angle_mu_draw[p,,prev_z[ p ], s] <- rnorm(n=2, mean=env_obj$mu[prev_z[ p ],,"mu",p,s], sd=sqrt(diag(env_obj$Qt[3:4,3:4,prev_z[ p ],  p,s]) * env_obj$mu[prev_z[ p ],,"V",p,s]))
			#logv_angle_draw[p,,prev_z[ p ]] <-  as.numeric(mvtnorm::rmvnorm(n=1, mean=logv_angle_mu_draw[p,,prev_z[ p ]], sigma=env_obj$Qt[3:4, 3:4, p]))		   
								
			#dont turn
			env_obj$logv_angle_mu_draw[p,"turn",,s] <- 0 
			
			#multiply by gradient since later will be variance of theta
			env_obj$Qt[4,4,,p,s] <- 1e-3 
			#else { env_obj$Qt[4,4]*(grad_ginv(psi=logv_angle_mu_draw[p,"turn",prev_z[ p ]])^2) }
			#logv_angle_draw[p,"turn", prev_z[ p ]] <- normalize_angle(log2rad(logv_angle_draw[p,"turn",prev_z[ p ]]))
			
			#go back to beginning of interval using speed and reverse the bearing.  then take that going forward.
			
			mk_tmp <- c(first_pos, env_obj$logv_angle_mu_draw[p,"logv",prev_z[ p ],s], normalize_angle(initial_bearing_forward - pi))
			
			env_obj$mk_actual[,p,s] <- env_obj$f(mk=mk_tmp, new_logv=env_obj$logv_angle_mu_draw[p,"logv",prev_z[ p ], s], 
												 theta=0, dtprev=initial_time_back) #a_{t+1}
			
			Fx_tmp <- env_obj$Fx(mk=mk_tmp, dtprev=initial_time_back)

			env_obj$Pk_actual[,,p,s] <- as.matrix(Matrix::nearPD(Fx_tmp %*% env_obj$Pk_actual[,,p,s] %*% t(Fx_tmp) + jtmp[1] * env_obj$Qt[,,prev_z[ p ],p,s], ensureSymmetry=TRUE)$mat) #R_{t+1}
			
					
			if (env_obj$truncate_to_map) {
				
				tmp <- reject_sampling(mu=env_obj$mk_actual[,p,s], cmat=env_obj$Pk_actual[,,p,s], prev_val=c(time0_estimate, env_obj$logv_angle_mu_draw[p,"logv",prev_z[ p ], s], normalize_angle(initial_bearing_forward - pi)), obj=env_obj)
				env_obj$Xpart_history[ env_obj$i, c("X","Y","logv","bearing_rad"), p, s] <- tmp$val
				num_reject[ p ] <- tmp$iter
			
			}
			else {
				env_obj$Xpart_history[ env_obj$i, c("X","Y","logv","bearing_rad"), p, s] <-  mvtnorm::rmvnorm(n=1, mean=env_obj$mk_actual[,p,s], sigma=env_obj$Pk_actual[,,p,s])  
			}

			# now reset Fx and Pk to refer to movement forward
			
			env_obj$mk_actual["bearing_rad",p,s] <- normalize_angle(env_obj$mk_actual["bearing_rad",p,s] + pi)
			
			Fx_tmp <- env_obj$Fx(mk=env_obj$mk_actual[,p,s], dtprev=env_obj$reg_dt)

			env_obj$Pk_actual[,,p,s] <- as.matrix(Matrix::nearPD(Fx_tmp %*% env_obj$Pk_actual[,,p,s] %*% t(Fx_tmp) + env_obj$Qt[,,prev_z[ p ],p,s], ensureSymmetry=TRUE)$mat) #R_{t+1}
	

		}#loop over part
		
		#turn 180 degrees around to go the correct original direction next time

		
		#normalize angle and correct since this is the direction moving forward
		env_obj$Xpart_history[ env_obj$i, "bearing_rad", , s] <-  normalize_angle(env_obj$Xpart_history[ env_obj$i, "bearing_rad", , s] + pi)

		if (env_obj$truncate_to_map) {
			env_obj$reject_samp_hist[ env_obj$i,,s] <- c(mean(num_reject), median(num_reject))
		}
						

		env_obj$Xpart_history[ env_obj$i, "region",,s ] <- apply(env_obj$Xpart_history[ env_obj$i , c("X","Y"),,s], 2, function(x) which_region(newcoord=x, centroid=env_obj$centroids))

		env_obj$Xpart_history[ env_obj$i, "time_in_state",,s ] <- 1
		
		err_tmp <- c()
		
		for (y in 1:nrow(env_obj$y_first)) {
			err_tmp <- cbind(err_tmp,  dist_func(center=env_obj$y_first[ y, c("X","Y"), drop=FALSE], otherXY=t(apply(env_obj$Xpart_history[env_obj$i, c("X","Y","logv","bearing_rad"),,s], 2, function(x) env_obj$h(mk=x, dtprev=jtmp[y] * env_obj$reg_dt)))))
		}
		env_obj$error_beforesamp_quantiles[ env_obj$i,,s] <-  quantile(rowSums(err_tmp), p=c(.1, .5, .9))
		
	}
	
	invisible(NULL)
	
}
	