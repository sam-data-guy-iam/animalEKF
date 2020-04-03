#http://stackoverflow.com/questions/34384907/how-can-put-multiple-plots-side-by-side-in-shiny-r



cdlm_robot_twostate <- function() {

shinyApp(

ui=shinyUI(fluidPage(
  titlePanel("Particle filter example of estimating mean velocities and location of 1-D robot\nfor two behavior types"),
  sidebarLayout(position = "left",
                             sidebarPanel("Simulation parameters", width=3,
							   sliderInput('iter',"Progress of simulation",value=1, min=1, max=50, round=TRUE, step=1,
				                 animate=animationOptions(interval=7000, loop=FALSE)),
		      				    actionButton('run',"Accept settings, press play above"),
								sliderInput('max_iter','Maximum number of iterations',value=20, min=1, max=30, round=TRUE, step=1),
								checkboxInput('sep_col',"Color particles differently (not recommended if more than 30 particles)",value=TRUE),								
								sliderInput('npart','Number of particles',value=6, min=2, max=40, round=TRUE, step=1),
								sliderInput('vel_mu1',"Unknown true mean of velocity (units/sec), type 1",value=4, min=-10, max=10, round=FALSE, step=0.1),
								sliderInput('vel_mu2',"Unknown true mean of velocity (units/sec), type 2",value=8, min=-10, max=10, round=FALSE, step=0.1),
								sliderInput('vel_var','Known true variance of velocity',value=1, min=0.1, max=5, round=FALSE, step=1),
								sliderInput('delta_t','Time step (sec)',value=1, min=0.1, max=10, round=FALSE, step=1),
								sliderInput('mu_mean1','Prior mean (units/sec) of velocity mean, type 1',value=4, min=-10, max=10, round=FALSE, step=1),
								sliderInput('mu_mean2','Prior mean (units/sec) of velocity mean, type 2',value=8, min=-10, max=10, round=FALSE, step=1),
								sliderInput('mu_var','Prior variance on velocity means',value=2, min=0.1, max=5, round=FALSE, step=1),
								sliderInput('p1to2','Transition probability between type 1 and 2',value=0.2, min=0, max=1, round=FALSE, step=0.05),
								sliderInput('p2to1','Transition probability between type 2 and 1',value=0.2, min=0, max=1, round=FALSE, step=0.05),
								checkboxInput('tp_known','Are transition probabilities known?',value=FALSE),
								textInput('dir_prior',"Dirichlet prior values on transition probability.\nEnter 4 positive numbers, comma separated.", value="8,2,2,8"),
							    textInput('xt_var',"Covariance matrix of position and velocity.\nEnter 4 numbers (comma delimited) for a covariance matrix, where the 4th should be the same as the velocity variance.  Coerced to PD matrix.",
							        value="0.1,0,0,0.1"),
								sliderInput('yt_var','Measurement error variance of true location',value=0.05, min=0.001, max=1, round=FALSE, step=0.1),
								textInput('Pk_init',"Covariance matrix of particles.\nEnter 4 numbers (comma delimited) for a covariance matrix.   Coerced to PD matrix.",
							        value="0.25,0,0,0.25")
							   
                              ),
                 mainPanel( 
                            splitLayout(cellWidths = c("50%", "50%"), plotOutput("densplot"), plotOutput("pred_loc")),
							splitLayout(cellWidths = c("50%", "50%"), plotOutput("wt_dist_sep"), plotOutput("pred_loc_resamp")),
							splitLayout(cellWidths = c("50%", "50%"), plotOutput("rugplot"), plotOutput("loc_hist")),
  						    splitLayout(cellWidths = c("50%", "50%"), plotOutput("trans_prob"),plotOutput("agree_hist"))

							 )
						)
               )#end of layout
  
)#end of UI definition
,

server=shinyServer(function(input, output, session) 
{
  
   observeEvent( input$run, {
        
		
	    prev_pars <- list(mai=par()$mai)
	    #plot.new()	
        updateSliderInput(session, "iter",  label="Progress of simulation", value=1, min=1, max=input$max_iter, step=1)
		#try(rm(params))
	    tmat <- matrix(c(1-input$p1to2, input$p1to2, input$p2to1, 1-input$p2to1), byrow=TRUE, ncol=2)
		
	   
		
		dir_prior <- abs(as.numeric(numextractall(input$dir_prior)))
	    dir_prior[ dir_prior==0 ] <- 1e-10
	
        yt <- matrix(0, ncol=4, nrow=input$max_iter +1)
        colnames(yt) <- c("true_v","true_loc","state","obs_loc")
		
		yt[1,"state"] <- sample(x=1:2, size=1, prob=c(.5,.5))#sample(x=1:2, size=1, prob=c(sum(dir_prior[c(1,3)]), sum(dir_prior[c(2,4)])))
		for (ii in 1:input$max_iter) {
		  yt[ii+1,"state"] <- sample(x=1:2, size=1, prob=tmat[yt[ii,"state"],])
		}
        state_table <- table(factor(yt[,"state"], levels=1:2))
		
		#at time t, state=current state, and true_v is the velocity
		
        if (state_table[1]>0) { yt[yt[,"state"]==1,"true_v"] <- rnorm(n=state_table[1], mean=input$vel_mu1, sd=sqrt(input$vel_var)) }
		if (state_table[2]>0) { yt[yt[,"state"]==2,"true_v"] <- rnorm(n=state_table[2], mean=input$vel_mu2, sd=sqrt(input$vel_var)) }
	    yt[,"true_loc"] <- cumsum(c(0, input$delta_t*yt[1:input$max_iter ,"true_v"]))
		
		#add measurement noise
		yt[-1,"obs_loc"] <- yt[-1,"true_loc"] + rnorm(n=input$max_iter, mean=0, sd=sqrt(input$yt_var))
		

	
		xt_var <- as.numeric(numextractall(input$xt_var))
		xt_var[1] <- abs(xt_var[1])
        xt_var[4] <- input$vel_var			 
		xt_var <- as.matrix(Matrix::nearPD(matrix(xt_var,ncol=2, byrow=TRUE), ensureSymmetry=TRUE)$mat)
		
		Pk_init <- as.numeric(numextractall(input$Pk_init))
        Pk_init[c(1,4)] <- abs(Pk_init[c(1,4)])
		Pk_init <- as.matrix(Matrix::nearPD(matrix(Pk_init,ncol=2, byrow=TRUE), ensureSymmetry=TRUE)$mat)
	
			
		params<- reactiveValues(dnames= list(c("orig","resamp"), paste("i",1:(input$max_iter+1),sep=""), paste("state",1:2,sep=""), paste("p",1:input$npart,sep=""))) 
		     
		params$Pk = array(Pk_init, dim=c(2,2,2,input$max_iter+1,2,input$npart), dimnames=c(list(1:2,1:2),params$dnames))
		params$mk = array(0, dim=c(2,1,2,input$max_iter+1,2,input$npart), dimnames=c(list(1:2,1),params$dnames))
		params$xpart = array(NA, dim=c(2,3,input$max_iter+1,input$npart), dimnames=c(list(1:2), list(c("orig","resamp","resamp_hist")), params$dnames[ -c(1,3) ]))
        params$Kgain =array(NA, dim=c(2,1,2,input$max_iter+1,2,input$npart), dimnames=c(list(1:2,1),params$dnames))
	    params$mk_prev = array(NA, dim=c(2,1,2,input$max_iter+1,2,input$npart), dimnames=c(list(1:2,1),params$dnames))
        params$Pk_prev = array(Pk_init, dim=c(2,2,2,input$max_iter+1,2,input$npart), dimnames=c(list(1:2,1:2),params$dnames))
		params$Ydist = array(NA, dim=c(2,3,input$max_iter+1,2,input$npart), dimnames=c(list(c("MuY","VarY"),c("orig","resamp","resamp_hist")), params$dnames[ -1 ]))
		params$Ydist_actual = array(NA, dim=c(3,2,input$max_iter+1,input$npart), dimnames=c(list(c("MuY","VarY","state"),c("resamp","resamp_hist")), params$dnames[c(2,4)]))
		params$curr_state=array(NA, dim=c(2,input$max_iter+1,input$npart), dimnames=params$dnames[ -3 ])
		params$mu_pars= array(NA, dim=c(2,2, input$max_iter+1,2,input$npart), dimnames=c(list(c("mu_mean","mu_var")),params$dnames))
		params$indices=array(NA,dim=c(input$max_iter+1, input$npart))
		params$wts=array(NA,dim=c(input$max_iter+1,2,input$npart), dimnames=params$dnames[-1])
		params$mu_pars["mu_mean","orig",,"state1",] <- input$mu_mean1
   		params$mu_pars["mu_mean","orig",,"state2",] <- input$mu_mean2
		params$mu_pars["mu_var","orig",,,] <- input$mu_var
		params$dir_params =array(NA, dim=c(4,2,input$max_iter+1,input$npart), dimnames=c(list(c("a11","a12","a21","a22")), params$dnames[ -3 ]))
		params$trans_draws =array(NA, dim=c(2,2,2,input$max_iter+1,input$npart), dimnames=c(list(c("1to","2to"),1:2),params$dnames[ -3 ]))
		
				
		params$resamp_colors = array(NA, dim=c(2,input$max_iter+1, input$npart), dimnames=params$dnames[ -3 ])
		params$is_new_step = rep(TRUE, input$max_iter)
	
  	    #initial guesses
			 
		params$mu_guess = array(NA, dim=c(2,input$max_iter+1,2,input$npart), dimnames=params$dnames)
		params$mu_guess["orig","i1",,] <- rnorm(n=2*input$npart, mean=as.vector(params$mu_pars["mu_mean","orig","i1",,]), sd=sqrt(as.vector(params$mu_pars["mu_var","orig","i1",,])))
		params$dir_params[,"orig","i1",] <- dir_prior
		
		#if known just use values
		if (input$tp_known) {
			params$trans_draws["1to",2,,,] <- input$p1to2
			params$trans_draws["1to",1,,,] <- 1-input$p1to2
						
			params$trans_draws["2to",1,,,] <- input$p2to1
			params$trans_draws["2to",2,,,] <- 1-input$p2to1
		}
		else {
			params$trans_draws["1to",1:2,"orig","i1",] <- apply(params$dir_params[c("a11","a12"),"orig","i1",], 2, function(x) MCMCpack::rdirichlet(n=1, alpha=x))
			params$trans_draws["2to",1:2,"orig","i1",] <- apply(params$dir_params[c("a21","a22"),"orig","i1",], 2, function(x) MCMCpack::rdirichlet(n=1, alpha=x))
		}
		
		trans_dens_mode <- c(params$dir_params[c("a12","a21"),"orig",1,1]-1)/c(params$dir_params[c("a12","a21"),"orig",1,1] + params$dir_params[c("a11","a22"),"orig",1,1]-2)
		trans_dens_mode_dens <- c(dbeta(x=trans_dens_mode[ 1 ], shape1=params$dir_params["a12","orig",1,1], shape2=params$dir_params["a11","orig",1,1]),
								  dbeta(x=trans_dens_mode[ 2 ], shape1=params$dir_params["a21","orig",1,1], shape2=params$dir_params["a22","orig",1,1]))
		trans_dens_range <- c(0, 2*max(trans_dens_mode_dens))						  
		 	
		
		
		
		
		#know which one you start in
		params$curr_state["orig","i1",] <- yt[1,"state"]
		
   	    if (input$sep_col==TRUE) { rainbow_cols <- colorspace::rainbow_hcl(n=input$npart +1, c=200,l=seq(40, 100, length.out=input$npart+1), alpha=0.5)[1:input$npart] }
        else { rainbow_cols <- rep("lightgray",input$npart) }			
        
		
	 	params$vel_guess = array(NA, dim=c(2,input$max_iter+1,2,input$npart), dimnames=params$dnames)
		params$vel_guess["orig","i1",,] <- rnorm(n=2*input$npart, mean=as.vector(params$mu_guess["orig","i1",,]), sd=sqrt(input$vel_var))
	    #for some reason it gives you a dumb error because it keeps input$iter as the last one for some reason  
     	   for (ii in 2:(input$max_iter+1)) { 
			    params$mu_guess["orig",ii,,] <- params$mu_guess["orig","i1",,] 
				params$vel_guess["orig",ii,,] <- params$vel_guess["orig","i1",,] 
				
			}
			 
		 #initialize the first time.  after 
	 
	    params$mk[2,1,"orig","i1",,] <- params$vel_guess ["orig","i1",,]
		params$xpart[,"orig","i1",] <- params$mk[,1,"orig","i1",yt[1,"state"],]
		
		
				
		loc_axis <- seq(0, ceiling(yt[input$max_iter+1,"obs_loc"]), by=round(max(input$max_iter/8, 1)*mean(c(input$vel_mu1, input$vel_mu2))))    
        xrange <- nice_range(x=yt[,"obs_loc"], ep=0.1)
		
		
		xdens_range1 <- input$vel_mu1 + c(-2,2)*sqrt(input$mu_var)
		xdens_range2 <- input$vel_mu2 + c(-2,2)*sqrt(input$mu_var)

			
		xdens_range <- nice_range(range(c(xdens_range1, xdens_range2)), ep=0.4)
	
		ydens_range <- c(0, 1.33*dnorm(0, 0, sqrt(input$vel_var)))
		part_axis <- seq(1, input$npart, by=min(5,round(input$npart/2)))
		
		f <- function(xt, dt=input$delta_t, newV) {
         c(xt[1]+dt*xt[2], newV)
        }
		Fx <- function(xt, dt=input$delta_t, newV) {
         matrix(c(1, dt, 0, 0), ncol=2, byrow=TRUE)
        } 
		#just keep one since just constant
		Fxmat <- Fx(xt=matrix(c(1,1), ncol=1), newV=1)
		
		h <- function(xt, dt=input$delta_t) {
         xt[1]+dt*xt[2]
        }
		Hx <- function(xt, dt=input$delta_t) {
        matrix(c(1, dt), nrow=1)
		}
    	#just keep one since just constant
		Hxmat <- Hx(xt=c(1,1))
		dens_pts <- seq(xdens_range[1], xdens_range[2], length.out=150) 

        dens_pts1 <- seq(xdens_range1[1], xdens_range1[2], length.out=75) 
		dens_pts2 <- seq(xdens_range2[1], xdens_range2[2], length.out=75) 
		plot.new()
        par(xpd=TRUE, las=1) 		
				
		 observeEvent( input$iter, { 
				    
					req(input$run)
					req(input$iter>1)	
					
					#print(paste("iter",input$iter))              
					#print(length(params$is_new_step))
					#print(dim(params$wts))
					if (params$is_new_step[ input$iter ] & input$iter >1) {
					   #a quick fix to set future values to TRUE
				       params$is_new_step[ (input$iter+1):(input$max_iter)] <- TRUE
					
					   #params$resamp_colors["orig",input$iter,order(apply(params$mu_guess["orig",input$iter,,],2,min))] <- rainbow_cols 
					   params$resamp_colors["orig",input$iter,] <- rainbow_cols
					
				    for (nn in 1:input$npart) {
                        for (kk in 1:2) {					
				      
							 params$mk_prev[,,"orig",input$iter,kk,nn] <- f(xt=params$mk[,,"orig",input$iter-1,kk,nn], newV=params$vel_guess["orig",input$iter,kk,nn]) 
							 params$Pk_prev[,,"orig",input$iter,kk,nn] <- as.matrix(Matrix::nearPD(Fxmat%*%params$Pk[,,"orig",input$iter,kk,nn]%*%t(Fxmat) + xt_var, ensureSymmetry=TRUE)$mat)
							 params$Ydist["MuY","orig",input$iter,kk,nn] 	<- h(xt=params$mk_prev[,,"orig",input$iter,kk,nn])
							 params$Ydist["VarY","orig",input$iter,kk,nn] 	<- as.numeric(Matrix::nearPD(Hxmat%*%params$Pk_prev[,,"orig",input$iter,kk,nn]%*%t(Hxmat) + input$yt_var, ensureSymmetry=TRUE)$mat)
        					}							
						}
					
                    }
					
					#density of particle velocities	
					output$densplot <- renderPlot({ 
					    #par(mfrow=c(1,2))
					    if (input$iter >0) {
							densplot_twostate(dpts=dens_pts, mu_guess=params$mu_guess["orig",input$iter,,], norm_sds=sqrt(input$vel_var), known_mean=c(input$vel_mu1, input$vel_mu2),
									colors=params$resamp_colors["orig",input$iter,], ylims=ydens_range, xlims=xdens_range, before_after="before",
									sep_col=input$sep_col, npart=input$npart)
									
							#densplot(dpts=dens_pts, mu_guess=params$mu_guess["orig",input$iter,"state1",], norm_sds=sqrt(input$vel_var), known_mean=input$vel_mu1,
							#		colors=params$resamp_colors["orig",input$iter,], ylims=ydens_range, xlims=xdens_range, before_after="before",
							#		sep_col=input$sep_col, npart=input$npart, iter=input$iter,
							#		main_title="Velocity guess type 1, t=")
											
									
							#densplot(dpts=dens_pts, mu_guess=params$mu_guess["orig",input$iter,"state2",], norm_sds=sqrt(input$vel_var), known_mean=input$vel_mu2,
							#		colors=params$resamp_colors["orig",input$iter,], ylims=ydens_range, xlims=xdens_range, before_after="before",
							#		sep_col=input$sep_col, npart=input$npart, iter=input$iter,
							#		main_title="Velocity guess type 2, t=")
								
						}
					})
					
					
					
					
					if (params$is_new_step[ input$iter ] & input$iter >1) {			
					
						#resample
						
						ynext <- yt[input$iter+1,"obs_loc"]
						params$wts[input$iter,"state1",] <- dnorm(x=yt[input$iter+1,"obs_loc"], mean=params$Ydist["MuY","orig",input$iter,"state1",], sd=sqrt(params$Ydist["VarY","orig",input$iter,"state1",]))
						params$wts[input$iter,"state2",] <- dnorm(x=yt[input$iter+1,"obs_loc"], mean=params$Ydist["MuY","orig",input$iter,"state2",], sd=sqrt(params$Ydist["VarY","orig",input$iter,"state2",]))
						#multiply by transition probabilities
						 #print(params$trans_draws[params$curr_state["orig",input$iter,1],1:2,"orig",input$iter,1])
							
							#note that for curr_state 
						for (nn in 1:input$npart) {
							 params$wts[input$iter,,nn] <- params$wts[input$iter,,nn]*params$trans_draws[ params$curr_state["orig",input$iter-1,nn],, "orig", input$iter-1,nn]
						}
						#params$wts[ input$iter, is.na(params$wts[input$iter,])] <- 1#1e-15
					
												

						qwt <- quantile(x=params$wts[input$iter,,], probs=0.85)
									
						params$indices[input$iter,] <- low_var_sample(wts=colSums(params$wts[input$iter,,]), M=input$npart)
						##sample(x=1:input$npart, size=input$npart, prob=colSums(params$wts[input$iter,,]), replace=TRUE)
						
						if (! input$sep_col) { params$resamp_colors["orig",input$iter, colSums(params$wts[input$iter,,]) >qwt ] <- "red" }
						
					}
				
                    
					# output$trans_prob <- renderPlot({
						# if (input$tp_known==FALSE) {
							# probability_trans(true_probs=c(input$p1to2, input$p2to1), dir_params=params$dir_params[,"orig",input$iter,], 
					                   # npart=input$npart, colors=params$resamp_colors["resamp",input$iter,])
						# }
					# })
					
					output$pred_loc <- renderPlot({
						if (input$iter >1) {	
						 pred_loc(Ydist1=params$Ydist[,"orig",input$iter,"state1",], wts=params$wts[input$iter,,], xlims=xrange, xticks=loc_axis, before_after="before", 
								  colors=params$resamp_colors["orig",input$iter,], npart=input$npart, sep_col=input$sep_col,  indices=1:input$npart,
								  yt=yt[ input$iter:(input$iter+1),"obs_loc"], Ydist2=params$Ydist[,"orig",input$iter,"state2",],
								  Yindex1=1:input$npart, Yindex2=1:input$npart)
						}	
					})
					
					
					
				
					
					# #next location
					
					#main problem: this is reacting too late because the locations are the same in both before/after plots
				   		
					if (params$is_new_step[ input$iter ] & input$iter >1) {
						ord <- params$indices[input$iter, ]
						
						
						#draw state according to weights
						
						
						
						#previous behavior			
						params$curr_state["resamp",input$iter-1, ] <- params$curr_state["orig",input$iter-1, ord]
						
						#new behavior
						params$curr_state["orig",input$iter, ] <- apply(params$wts[input$iter,,params$indices[input$iter, ]], 2, function(y) low_var_sample(wts=y, M=1))
						
						
						
						params$mk[,,"resamp",input$iter-1,,] <- params$mk[,,"orig",input$iter-1,, ord]
						params$mk_prev[,,"resamp",input$iter,,] <- params$mk_prev[,,"orig",input$iter,, ord]
						params$Pk[,,"resamp",input$iter-1,,] <- params$Pk[,,"orig",input$iter-1,, ord]
						params$Pk_prev[,,"resamp",input$iter,,] <- params$Pk_prev[,,"orig",input$iter,, ord]
						params$mu_guess["resamp",input$iter,,] <- params$mu_guess["orig",input$iter,, ord]
						params$vel_guess["resamp",input$iter,,] <- params$vel_guess["orig",input$iter,, ord]
						params$Kgain[,,"resamp",input$iter-1,,] <- params$Kgain[,,"orig",input$iter-1,, ord]
						
						params$Ydist[,"resamp",input$iter,,] <- params$Ydist[,"resamp_hist",input$iter,,] <- params$Ydist[,"orig",input$iter,, ord]      
						params$Ydist[,"resamp",1:(input$iter-1),,] <- params$Ydist[,"resamp",1:(input$iter-1),,ord]
						
						params$xpart[,"resamp",input$iter-1,] <- params$xpart[,"resamp_hist", input$iter-1,] <- params$xpart[,"orig",input$iter-1, ord]
											
						
						if (input$iter >2 ) {
							
							params$xpart[,"resamp", 1:(input$iter-2),] <- params$xpart[,"resamp",1:(input$iter-2), ord]

						}	
						
						params$mu_pars[,"resamp",input$iter,,] <- params$mu_pars[,"orig",input$iter,, ord]
						
						params$resamp_colors["resamp",input$iter,] <- params$resamp_colors["orig",input$iter, ord]
						
						params$dir_params[,"resamp",input$iter-1,] <- params$dir_params[,"orig",input$iter-1,ord]
						if (input$tp_known==FALSE) { params$trans_draws[,,"resamp",input$iter-1,] <- params$trans_draws[,,"orig",input$iter-1,ord] }
						
						
						curr_state <- params$curr_state["orig",input$iter,]
						
						
						params$Ydist_actual["state","resamp_hist",input$iter,] <- params$Ydist_actual["state","resamp",input$iter,] <- curr_state 

						for (nn in 1:input$npart) {
							  params$Ydist_actual[c("MuY","VarY"),"resamp_hist",input$iter,nn] <- params$Ydist_actual[c("MuY","VarY"),"resamp",input$iter,nn] <- params$Ydist[c("MuY","VarY"),"orig",input$iter,curr_state[nn],nn]
						}
						params$Ydist_actual[,"resamp",1:(input$iter-1),] <- params$Ydist_actual[,"resamp",1:(input$iter-1),ord]
						
						
						
					}
					
					wts_sep_counter <<- locs_counter <<- 0 
					
					
					#distribution of weights
					output$wt_dist_sep <- renderPlot({
					   #par(mfrow=c(1,2))
					   #should have a panel for each state and combined
					   if (input$iter >1) {	
							
					   
					   
							if (params$is_new_step[ input$iter ]==FALSE & params$is_new_step[ input$iter+1 ]==TRUE) {
								#print(params$curr_state["resamp", input$iter, ])
								
								wts_sep_counter <<- wts_sep_counter +1
																
								
								if ( wts_sep_counter >=1 & wts_sep_counter< input$npart) {  invalidateLater(millis=3000/input$npart) }
									# a_index <- rep(NA, 2)
								
								
									# a_index[ params$curr_state["resamp", input$iter, wts_sep_counter] ] <- wts_sep_counter
									
									# #print(ifelse(params$curr_state["resamp", input$iter, wts_counter]==1, wts_counter, NULL))
									# m1 <- wt_dist_loop(wts=params$wts[input$iter,"state1",], xlims=xdens_range, ylims=ydens_range, known_mean=input$vel_mu1, mu_guess=params$mu_guess["orig", input$iter,"state1",], 
													# npart=input$npart, colors=params$resamp_colors["orig",input$iter,], sep_col=input$sep_col,
													# main_title=ifelse(input$sep_col, "Resamp. weights, type 1", "Weights, type 1\n"), 
													# index=a_index[ 1 ])
													
									# m1
									# #print(ifelse(params$curr_state["resamp", input$iter, wts_counter]==2, wts_counter, NULL))				
									# m2 <- wt_dist_loop(wts=params$wts[input$iter,"state2",], xlims=xdens_range, ylims=ydens_range, known_mean=input$vel_mu2, mu_guess=params$mu_guess["orig", input$iter,"state2",], 
													# npart=input$npart, colors=params$resamp_colors["orig",input$iter,], sep_col=input$sep_col,
													# main_title=ifelse(input$sep_col, "Resamp. weights, type 2", "Weights, type 2\n"), 
													# index=a_index[ 2 ]) 
							
										
									# m2
									
									m <- wt_dist_twostate_loop(wts=params$wts[input$iter,,], xlims=xdens_range, ylims=ydens_range, known_mean=c(input$vel_mu1, input$vel_mu2), 
																mu_guess=params$mu_guess["orig", input$iter,,], iter=input$iter,
																npart=input$npart, colors=params$resamp_colors["orig",input$iter,], sep_col=input$sep_col, 
																index=params$indices[ input$iter, wts_sep_counter], var_name="velocity",
																behavior=params$curr_state["orig", input$iter, wts_sep_counter])
													
									m
									
							}	
							else {
								
								 wt_dist_twostate(wts=params$wts[input$iter,,], xlims=xdens_range, ylims=ydens_range, known_mean=c(input$vel_mu1, input$vel_mu2), 
																mu_guess=params$mu_guess["orig", input$iter,,], iter=input$iter,  var_name="velocity",
																npart=input$npart, colors=params$resamp_colors["orig",input$iter,], sep_col=input$sep_col,
																behavior=yt[input$iter,"state"])
													
										
														
							
							}
						
						}	
					})
					
					
					# #weight sums
					# output$wt_dist <- renderPlot({
					   # #par(mfrow=c(1,2), mai = c(1.02, 0.4, 0.82, 0.4))
					   # #should have a panel for each state and combined
					   # if (input$iter >1) {	
					   
							
							# if (params$is_new_step[ input$iter ]==FALSE & params$is_new_step[ input$iter+1 ]==TRUE) {
								# #print(params$curr_state["resamp", input$iter, ])
								
								# wts_counter <<- wts_counter +1
																
								
								# if ( wts_counter >=1 & wts_counter< input$npart) {  pause_between(3000/input$npart, session) }
									
																		
									# #print(ifelse(params$curr_state["resamp", input$iter, wts_counter]==1, wts_counter, NULL))
									# m <- wt_dist_loop(wts=colSums(params$wts[input$iter,,]), xlims=xdens_range, ylims=ydens_range, known_mean=c(input$vel_mu1, input$vel_mu2)[ yt[input$iter,"state"] ], 
													  # mu_guess=colSums(params$mu_guess["orig", input$iter,,]*apply(params$wts[input$iter,,], 2, function(x) x/sum(x))), 
													  # npart=input$npart, colors=params$resamp_colors["orig",input$iter,], sep_col=input$sep_col,
													  # main_title=paste("Distribution of TOTAL PARTICLE resampling weights\nbased on predictions of next location, t=",input$iter), 
													  # index=params$indices[ input$iter, wts_counter])
													
									# m
									
							# }	
							# else {
								
								# wt_dist(wts=colSums(params$wts[input$iter,,]), xlims=xdens_range, ylims=ydens_range, known_mean=c(input$vel_mu1, input$vel_mu2)[ yt[input$iter,"state"] ], 
										 # mu_guess=colSums(params$mu_guess["orig", input$iter,,]*apply(params$wts[input$iter,,], 2, function(x) x/sum(x))), 
										 # npart=input$npart, colors=params$resamp_colors["orig",input$iter,], sep_col=input$sep_col,
										 # main_title=paste("Distribution of TOTAL PARTICLE resampling weights\nbased on predictions of next location, t=",input$iter))									
							
							# }
						
						# }	
					# })
					
					
						# output$wt_dist <- renderPlot({
					   # # #should have a panel for each state and combined
						# if (input$iter >0 ) {
					   
							# par(mfrow=c(1,1))
								# wt_dist(wts=colSums(params$wts[input$iter,,]), xlims=xdens_range, known_mean=c(input$vel_mu1, input$vel_mu2)[ yt[input$iter,"state"] ], 
										# mu_guess=colSums(params$mu_guess["orig", input$iter,,]*apply(params$wts[input$iter,,], 2, function(x) x/sum(x))), 
										# npart=input$npart, colors=params$resamp_colors["orig",input$iter,], sep_col=input$sep_col)				
						# }
					 # })
					
					
					
					
					output$pred_loc_resamp <- renderPlot({
						if (input$iter >1) {	
						
							
							curr_state <- params$curr_state["orig",input$iter, ]
							
							if (params$is_new_step[ input$iter ]==FALSE & params$is_new_step[ input$iter+1 ]==TRUE) {
							
								locs_counter <<- locs_counter +1
							
								if ( locs_counter >=1 & locs_counter< input$npart) { invalidateLater(millis=3000/input$npart) }
							
									m <- pred_loc_loop(Ydist1=params$Ydist[,"resamp",input$iter,"state1",], wts=params$wts[input$iter,,], 
													xlims=xrange, xticks=loc_axis, before_after="after", colors=params$resamp_colors["resamp",input$iter,], 
													npart=input$npart, sep_col=input$sep_col, yt=yt[ input$iter:(input$iter+1),"obs_loc"], 
													Ydist2=params$Ydist[,"resamp",input$iter,"state2",], Yindex1=which(curr_state==1),
													indices=params$indices[input$iter,], index=locs_counter, Yindex2=which(curr_state==2))

									m
								
							
							}
							else {
								
								pred_loc(Ydist1=params$Ydist[,"resamp",input$iter,"state1",], wts=params$wts[input$iter,,], xlims=xrange, xticks=loc_axis, before_after="after", 
								  colors=params$resamp_colors["resamp",input$iter,], npart=input$npart, sep_col=input$sep_col, 
								  yt=yt[ input$iter:(input$iter+1),"obs_loc"], Ydist2=params$Ydist[,"resamp",input$iter,"state2",],
								  indices=params$indices[input$iter,], Yindex1=which(curr_state==1), Yindex2=which(curr_state==2))
													
							
							}
										
						   
						
						}		
					 })
					
					
	              			
					if (params$is_new_step[ input$iter ] & input$iter >1) {
					
					
						
						#update distributions
					
						for (nn in 1:input$npart) {
							#now draw for next time
							
							z <- params$curr_state["orig",input$iter,nn]
							#only update for the predicted state that actually happened
						
							params$Kgain[,,"orig",input$iter,z, nn ] <- params$Pk_prev[,,"resamp",input$iter,z, nn ]%*%t(Hxmat)%*%(1/params$Ydist["VarY","resamp",input$iter,z, nn ])
							params$mk[,,"orig",input$iter,z, nn ] <- params$mk_prev[,,"resamp",input$iter,z, nn ] + params$Kgain[,,"orig",input$iter,z, nn ]*(ynext-h(xt=params$mk_prev[,,"resamp",input$iter,z, nn ]))
							params$Pk[,,"orig",input$iter,z, nn ] <- as.matrix(Matrix::nearPD(params$Pk_prev[,,"resamp",input$iter,z, nn ] - params$Kgain[,,"orig",input$iter,z, nn ]%*%(params$Ydist["VarY","resamp",input$iter, z, nn ]%*%t(params$Kgain[,,"orig",input$iter,z, nn ])), ensureSymmetry=TRUE)$mat)
							
							z_oth <- which(1:2 !=z)
							params$Kgain[,,"orig",input$iter,z_oth, nn ] <- params$Kgain[,,"resamp",input$iter,z_oth, nn ]
							
							#for mk, take the actual predicted state and use it to update the other, so that start at the same x location. 
							#the velocity will be done to make mk_prev different
							#for Pk and Kgain only update the one that happened
							params$mk[,,"orig",input$iter,z_oth, nn ] <- params$mk[,,"orig",input$iter,z, nn ]
							params$Pk[,,"orig",input$iter,z_oth, nn ] <- params$Pk[,,"orig",input$iter,z, nn ]
							
						
						}
						
						#now propagate next values
						
						curr_state <- params$curr_state["orig",input$iter,]
						prev_state <- params$curr_state["resamp",input$iter-1,]
						alpha_index <- paste("a",prev_state,curr_state,sep="")
					
						#set mu pars for next to current, then overwrite with the one observed
						
						params$mu_pars["mu_mean","orig",input$iter+1,,] <- params$mu_pars["mu_mean","resamp",input$iter,,]
						params$mu_pars["mu_var","orig",input$iter+1,,] <-  params$mu_pars["mu_var","resamp",input$iter,,]
					
					
						for (nn in 1:input$npart) {
						  
						   z <- curr_state[ nn ]
						   params$xpart[,"orig",input$iter, nn ] <- mvtnorm::rmvnorm(n=1, mean=params$mk[,,"orig",input$iter,z,nn], sigma=params$Pk[,,"orig",input$iter,z,nn])
                           # #now update values of parameters for state observed
					       
				           params$mu_pars["mu_var","orig",input$iter+1,z,nn] <- 1/(1/params$mu_pars["mu_var","resamp",input$iter,z,nn] + 1/input$vel_var)
					       params$mu_pars["mu_mean","orig",input$iter+1,z,nn] <- (params$mu_pars["mu_mean","resamp",input$iter,z,nn]/params$mu_pars["mu_var","resamp",input$iter,z,nn] + params$xpart[2,"orig",input$iter,nn]/input$vel_var)*params$mu_pars["mu_var","orig",input$iter+1,z,nn]
					       #copy and update only value
						   
						    if (input$tp_known==FALSE) {
							   params$dir_params[,"orig",input$iter,nn] <- params$dir_params[,"resamp",input$iter-1, nn ]
							   #params$dir_params[ alpha_index[ nn ],"orig",input$iter+1,nn] <- params$dir_params[alpha_index[ nn ],"orig",input$iter+1,nn] +1
							   params$dir_params[ alpha_index[ nn ],"orig",input$iter,nn] <- params$dir_params[alpha_index[ nn ],"orig",input$iter-1,nn] +1
						    }
						}
                
                    }
					
												 
							
				    # #new parameter draws
					for (kk in 1:2) {
						params$mu_guess["orig",input$iter+1,kk,] <- rnorm(n=input$npart, mean=params$mu_pars["mu_mean","orig",input$iter+1,kk,], sd=sqrt(params$mu_pars["mu_var","orig",input$iter+1,kk,]))
						params$vel_guess["orig",input$iter+1,kk,] <- rnorm(n=input$npart, mean=params$mu_guess["orig",input$iter+1,kk,], sd=sqrt(input$vel_var))      
					}
					
				    if (input$tp_known==FALSE) {
						params$trans_draws["1to",1:2,"orig",input$iter,] <- apply(params$dir_params[c("a11","a12"),"orig",input$iter,], 2, function(x) MCMCpack::rdirichlet(n=1, alpha=x))
						params$trans_draws["2to",1:2,"orig",input$iter,] <- apply(params$dir_params[c("a21","a22"),"orig",input$iter,], 2, function(x) MCMCpack::rdirichlet(n=1, alpha=x))
					}
					
					output$rugplot <- renderPlot({
					 
						if (input$iter >0) {	
							
						 #par(mfrow=c(1,2))
						  convergence_rugplot_twostate(xlims=xdens_range, known_mean=c(input$vel_mu1, input$vel_mu2), max_iter=input$max_iter, iter=input$iter, 
									 mu_guess=params$mu_guess, npart=input$npart, colors=params$resamp_colors) 
						 
							# convergence_rugplot(xlims=xdens_range, known_mean=input$vel_mu1, max_iter=input$max_iter, iter=input$iter, 
									# mu_guess=params$mu_guess[,,"state1",], npart=input$npart, 
									# main_title="Convergence of velocity\nmean, type 1") 
									
							# convergence_rugplot(xlims=xdens_range, known_mean=input$vel_mu2, max_iter=input$max_iter, iter=input$iter, 
									# mu_guess=params$mu_guess[,,"state2",], npart=input$npart,
									# main_title="Convergence of velocity\nmean, type 2") 
						}					
					})
					output$loc_hist <- renderPlot({
					    if (input$iter >1) {	
						
							omega_tmp <- params$xpart[1,"resamp",1:(input$iter-1),]
							omega_tmp <- rbind(omega_tmp, params$xpart[1,"orig",input$iter,])
							
							location_history_v2(omega=omega_tmp, xticks=loc_axis, xlims=xrange, iter=input$iter, 
					                     npart=input$npart, yt=yt[,"true_loc"], colors=params$resamp_colors["resamp",,])
						
						#location_history(Ydist=params$Ydist_actual[,"resamp_hist",,], xticks=loc_axis, xlims=xrange, iter=input$iter, 
					    #                 npart=input$npart, colors=params$resamp_colors["resamp",,], yt=yt)
					  }
					})	
                   
				    output$agree_hist <- renderPlot({
						if (input$iter >1) {		
										       
							#state_agreement(particle_states=params$Ydist_actual["state","resamp_hist",,], actual_states=yt[,"state"], 
							#				iter=input$iter, npart=input$npart) 
							
							states_tmp <- params$curr_state["resamp",1:(input$iter-1),]
							states_tmp <- rbind(states_tmp, params$curr_state["orig",input$iter:(input$max_iter+1),])
							
							
							state_agreement(particle_states=states_tmp, actual_states=yt[,"state"], 
											iter=input$iter, npart=input$npart) 
						}		
					})

					output$trans_prob <- renderPlot({
						if (input$iter >1 & input$tp_known==FALSE) {		
									
							probability_trans(true_probs=c(input$p1to2, input$p2to1), dir_params=params$dir_params[,"orig", input$iter,],
  										      npart=input$npart, colors=params$resamp_colors["resamp",input$iter,], ylims=trans_dens_range)
						}		
					})

									
				
  				params$is_new_step[ input$iter ] <- FALSE
					
				 })#end of iteration 
	
	 
		  }) #if run
		  
 }
)
)

}

