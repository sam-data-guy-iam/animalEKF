fix_data_EKF_interp_joint <- function(env_obj) {

	env_obj$d <- env_obj$d[ order(env_obj$d$date_as_sec),]
	env_obj$first_time <- min(env_obj$d$date_as_sec)

	env_obj$t_reg <- seq(from=env_obj$first_time, by=env_obj$reg_dt, length.out=min(env_obj$maxStep + 1, ceiling(1 + (env_obj$d$date_as_sec[ nrow(env_obj$d) ] - env_obj$first_time)/env_obj$reg_dt)))
	#do this so the first interval captures the first observation at t=0
	env_obj$t_reg[ 1 ] <- env_obj$t_reg[ 1 ] - env_obj$reg_dt * 0.0001  #.Machine$double.eps

	env_obj$N <- length(env_obj$t_reg) 

	env_obj$d <- env_obj$d[ env_obj$d$date_as_sec <= env_obj$t_reg[ env_obj$N ] ,]

	env_obj$tags <- env_obj$d$tag

	env_obj$included_intervals <- 1:(env_obj$N - 1)

	#calculate which regular step each observation falls into
	#we do right=TRUE so that if this is regular intervals, each observation is the next step (j= 1 rather than 0).
	env_obj$d$t_intervals <- as.numeric(as.character(cut(x=env_obj$d$date_as_sec, breaks=env_obj$t_reg, labels=env_obj$included_intervals, right=TRUE)))

	print("t_intervals")
	print(env_obj$d$t_intervals)

	env_obj$shark_names <- as.character(sort(unique(env_obj$tags)))
	print(paste("shark names are",paste(env_obj$shark_names, collapse=" ")))
	env_obj$shark_intervals <- list()
	env_obj$shark_valid_steps <- list()



	for (s in env_obj$shark_names) {
		env_obj$shark_intervals[[ s ]] <- unique(env_obj$d$t_intervals[ env_obj$d$tag==s ])

		#keep steps where there are less than a certain gap between observations 
		#excluded_steps <- env_obj$shark_intervals[[ s ]][ -1 ][ diff(env_obj$shark_intervals[[ s ]]) > max_int_wo_obs ]
		
		#env_obj$shark_intervals[[ s ]] <- c(env_obj$shark_intervals[[ s ]][ 1 ], env_obj$shark_intervals[[ s ]][ -1 ][ diff(env_obj$shark_intervals[[ s ]]) <= max_int_wo_obs ])
		
		#env_obj$shark_valid_steps[[ s ]] <- min(env_obj$shark_intervals[[ s ]]):max(env_obj$shark_intervals[[ s ]])
		#env_obj$shark_valid_steps[[ s ]] <- env_obj$shark_valid_steps[[ s ]][ ! env_obj$shark_valid_steps[[ s ]] %in% excluded_steps ]
		#env_obj$shark_valid_steps[[ s ]] <- env_obj$shark_valid_steps[[ s ]][ ! env_obj$shark_valid_steps[[ s ]] %in% excluded_steps ]
	  
	  
		if (length(env_obj$shark_intervals[[ s ]]) > 1) {
	  
			
			tmp <- c()
			tmp1 <- c()
		
	  
			print(env_obj$shark_intervals[[ s ]])
			for (jj in 1:(length(env_obj$shark_intervals[[ s ]])-1)){
				if(diff(env_obj$shark_intervals[[ s ]][ jj:(jj+1) ]) <= env_obj$max_int_wo_obs) {
					tmp <- c(tmp, (env_obj$shark_intervals[[ s ]][ jj ]):(env_obj$shark_intervals[[ s ]][ jj+1 ]))
					tmp1 <- c(tmp1, env_obj$shark_intervals[[ s ]][ jj:(jj+1) ])  
				}
			}
			
			#valid intervals to simulate for
			env_obj$shark_valid_steps[[ s ]] <- sort(unique(tmp))
			#only the intervals with observations
			env_obj$shark_intervals[[ s ]] <- sort(unique(tmp1))
		}
		else {
			env_obj$shark_names[ env_obj$shark_names==s ] <- NA
		}
		
	}

	env_obj$shark_names <- env_obj$shark_names[ ! is.na(env_obj$shark_names) ]

	#env_obj$shark_intervals <- env_obj$shark_intervals[ sapply(env_obj$shark_intervals, function(x) length(x) >1) ]
	#env_obj$shark_names <- names(env_obj$shark_intervals)
	env_obj$nsharks <- length(env_obj$shark_names)
	#env_obj$shark_valid_steps <- env_obj$shark_valid_steps[ env_obj$shark_names ]

	env_obj$included_intervals <- sort(unique(unlist(env_obj$shark_valid_steps))) 

	print("valid observations per shark:")
	print(env_obj$shark_valid_steps)

	print(paste("sharks:", paste(env_obj$shark_names, collapse=" ")))


	env_obj$first_intervals <- lapply(env_obj$shark_valid_steps, function(x) x[ !((x-1) %in% x) ])


	names(env_obj$first_intervals) <- env_obj$shark_names
	print("starting observations per shark:")
	print(env_obj$first_intervals)

	#last interval with a valid observation
	env_obj$shark_final_obs <- sapply(env_obj$shark_intervals, max)
	names(env_obj$shark_final_obs) <- env_obj$shark_names


	if (env_obj$nsharks==1) {
		env_obj$interact <- FALSE
	}

	if (env_obj$nstates==1) {
		env_obj$next_states <- env_obj$states <- rep(1, nrow(env_obj$d))
		env_obj$interact <- FALSE
	}

	if (is.null(env_obj$area_map)) {
		env_obj$truncate_to_map <- FALSE
	}

	if (env_obj$truncate_to_map==FALSE) {
		env_obj$do_trunc_adjust <- FALSE
	}
	
	
	if (env_obj$update_params_for_obs_only) {
		env_obj$update_eachstep <- FALSE		
	}

	print(paste("env_obj$nstates:", env_obj$nstates))
	print(paste("env_obj$interactions", env_obj$interact))


	#keep the initial interval, plus any that come after more than max_int_wo_obs observations
	#env_obj$first_intervals <- lapply(env_obj$shark_intervals, function(x) c(x[1], x[-1][ diff(x) > max_int_wo_obs ]) )
	#second interval
	#second_intervals <- sapply(env_obj$shark_intervals, function(x) x[2])

	#names(env_obj$first_intervals) <- env_obj$shark_namesenv_obj$
	#print("starting observations per shark:")
	#print(env_obj$first_intervals)

	env_obj$shark_symbols <- 1:env_obj$nsharks
	names(env_obj$shark_symbols) <- env_obj$shark_names

	#regions
	env_obj$XY <- as.matrix(env_obj$d[,c("X","Y")])

	env_obj$centroids <- as.matrix(env_obj$centroids, ncol=2)
	rownames(env_obj$centroids) <- NULL
	colnames(env_obj$centroids) <- NULL
	env_obj$regions <- apply(env_obj$XY[,,drop=FALSE], 1, function(x) which_region(x, centroid=env_obj$centroids))
	env_obj$nregions <- nrow(env_obj$centroids)

	env_obj$d <- env_obj$d[, colnames(env_obj$d) != "region" ]

	if(env_obj$nregions > 1) { env_obj$d <- cbind(env_obj$d, region=env_obj$regions) }
	else { env_obj$d <- cbind(env_obj$d, region=1) }


	env_obj$pnames <- paste("p", 1:env_obj$npart, sep="")
	env_obj$state_names <- paste("state", 1:env_obj$nstates, sep="")
	env_obj$rnames <- paste("r", 1:env_obj$nregions, sep="")
	env_obj$Nnames <- paste("N", 1:env_obj$N, sep="")

	#lookup vector for 11 12 21 22
	#trans_names <- 1:(env_obj$nstates^2)
	#names(trans_names) <- do.call(paste, c(expand.grid(1:env_obj$nstates,1:env_obj$nstates)[,2:1], sep=""))
	env_obj$trans_names <- as.character(do.call(paste, c(expand.grid(1:env_obj$nstates, 1:env_obj$nstates)[,2:1], sep="")))

	env_obj$d <- cbind(env_obj$d, shark_obs_index= NA)
	for (s in env_obj$shark_names) {
		ss <- which(env_obj$tags==s)
		env_obj$d[ss,"shark_obs_index"] <- 1:length(ss)
	}

	#print(d$shark_obs_index)
	#if (env_obj$nsharks>1) { print(env_obj$tags) }  

	env_obj$d <- cbind(env_obj$d, rowid=1:nrow(env_obj$d))

	print("intervals with observations per shark:")
	print(env_obj$shark_intervals)


	print("intervals to be simulated per shark:")
	print(env_obj$shark_valid_steps)
	#drop env_obj$tags and convert to matrix for easier computation


	#if we want to model it as one state, so be it
	
	if (env_obj$nstates == 1) {
		env_obj$d$state.guess2[ ! is.na(env_obj$d$state.guess2)] <- 1		
		env_obj$d$next.guess2[ ! is.na(env_obj$d$next.guess2)] <- 1
		env_obj$d$lambda[ ! is.na(env_obj$d$lambda)] <- 1

	}


	env_obj$d$state.guess2 <- as.numeric(env_obj$d$state.guess2)
	env_obj$d$next.guess2 <- as.numeric(env_obj$d$next.guess2)
	env_obj$d$lambda <- as.numeric(env_obj$d$lambda)


	env_obj$states <- as.numeric(env_obj$d$state.guess2)
	env_obj$next_states <- as.numeric(env_obj$d$next.guess2)


	#env_obj$nstates <- max(length(unique(states)), env_obj$nstates)

	print(sort(colnames(env_obj$d)))

	env_obj$d <- env_obj$d[,c("rowid","shark_obs_index","X","Y","logvelocity","date_as_sec",
							   "angle_velocity","turn.angle.rad","bearing.to.east.tonext.rad","region","time_to_next","dx_to_next",
							   "dy_to_next","lambda","state.guess2","next.guess2","t_intervals")]


	rownames(env_obj$d) <- 1:nrow(env_obj$d)
	env_obj$d <- as.matrix(env_obj$d)
	
	env_obj$obs_XY_bbox <- apply(env_obj$d[,c("X","Y"), drop=FALSE], 2, function(x) range(x, na.rm=TRUE))
	env_obj$obs_XY_bbox[,"X"] <- env_obj$obs_XY_bbox[,"X"] + c(-0.5, 0.5)*diff(env_obj$obs_XY_bbox[,"X"])  
	env_obj$obs_XY_bbox[,"Y"] <- env_obj$obs_XY_bbox[,"Y"] + c(-0.5, 0.5)*diff(env_obj$obs_XY_bbox[,"Y"])  
	
	env_obj$obs_XY_bbox[1,"X"] <- max(env_obj$obs_XY_bbox[1,"X"], env_obj$bbox[1,1])
	env_obj$obs_XY_bbox[2,"X"] <- min(env_obj$obs_XY_bbox[2,"X"], env_obj$bbox[1,2])

	env_obj$obs_XY_bbox[1,"Y"] <- max(env_obj$obs_XY_bbox[1,"Y"], env_obj$bbox[2,1])
	env_obj$obs_XY_bbox[2,"Y"] <- min(env_obj$obs_XY_bbox[2,"Y"], env_obj$bbox[2,2])
	
	invisible(NULL)
	
	nus <- length(unique(env_obj$states, na.rm=TRUE))
	nust <- env_obj$nstates
	if (env_obj$compare_with_known) {
		nust <- length(unique(env_obj$known_regular_step_ds$state.guess2, na.rm=TRUE))
	}
	
	
	if (nus != env_obj$nstates || nust != env_obj$nstates) {
		print(paste("Observed/true data has", nus, "and", nust, "behaviors, but choose to model with", env_obj$nstates, "behaviors"))
	}
	
}