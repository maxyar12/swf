!----------------------------------------------------------------------------------------------------------------
!
!  Square well fluid  
!	(length L=Rsw*Nx); Usw= square well depth/T
!       mu= chem potential /T + step_sc*rank.    
!       Linear bin length = lambda (by constuction)
!----------------------------------------------------------------------------------------------------------------     
	PROGRAM MAIN 
        USE MUMBERS      
                                         
	IMPLICIT NONE
	include 'mpif.h' 
     
!----------------------------Global_Vars-----------------------------------------------------------------------
	integer                              ::  Nx,Ny,Nparts,Npairs,Nmax,p,n_j,i,j,ii,iii,jj,jjj
	integer                              ::  Nbins,bin_max,lin_Nbins,cbin,rbin,bin,tag,bin_sw,reset_code
	logical                              ::  bounce
	integer, parameter                   ::  d=2
	double precision                     ::  mu_c
	double precision, parameter          ::  usw_c = 1.805054d0
	double precision                     ::  Usw,mu,Rhc,Rsw,step_sc,L,Z,x_j,y_j,acc_r,lin_bin_size,l_cx,phi,radphi
	double precision, allocatable        ::  Rx(:),Ry(:),Hist_N(:),Hist_P(:),Hist_N_P(:,:)
	integer, allocatable           ::  bin_count(:),bin_parts(:,:),neighbors(:,:)
	double precision       ::  dx, dy, dist,r_attempt,r_accept,c_accept,c_attempt,Q_Np2_Npts,Q_Np3_Npts,Q_Np4_Npts
	double precision          ::  Q_15,Q_16,Q_17,Q_18,Q_19,Q_20,a_1,a_2,Q_21,Q_22,Q_23,Q_24,Q_25,l_cx_check
	double precision        ::  amax, tmax, amin, tmin,step,s_step,p_step,t_step,i_t,i_p,i_s,Q_1,Q_2,Q_3
	double precision          ::  Q_4,Q_5,Q_6,Q_7,Q_Nparts,Q_Nparts_sq,Q_Npairs,Q_Npairs_sq,Q_Nparts_Npairs
	double precision       ::  Q_Npts2_Np,Q_Npts3_Np,Q_Npts4_Np,Q_11,Q_12,Q_14,Q_i_Nparts,Q_i_Npairs
	double precision          ::  Q_Nparts_c,Q_Nparts_q,Q_Npairs_c,Q_Npairs_q,Q_Nparts_f,Q_Npairs_f,Q_8,Q_9,Q_10
    integer                 ::  config_mode,prnt_number
    integer, parameter     ::  prt_m=10000                                          
   double precision     :: Q_1_prntouts(prt_m),Q_2_prntouts(prt_m),Q_3_prntouts(prt_m) 
   double precision     :: Q_6_prntouts(prt_m),Q_7_prntouts(prt_m),Q_10_prntouts(prt_m),Q_14_prntouts(prt_m)
   double precision       :: Q_4_prntouts(prt_m),Q_5_prntouts(prt_m),Q_17_prntouts(prt_m),Q_20_prntouts(prt_m) 
    character*256    :: results_name,cnf_name,stats_name,fname,error_name,mystr,histname  
          integer :: rank, comm_size, ierr, outproc, status(MPI_STATUS_SIZE)

        call MPI_INIT( ierr )
        call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, comm_size, ierr )   

!------------INPUTS from PAR file-------------------------------------------------------------
         CALL READ_INPUT_FILE         
!---------------------------------------------------------------------------    
          write(6,*) Nx
          write(6,*) phi
          write(6,*) l_cx, ' ', l_cx_check
		 write(6,*) Usw
		 write(6,*) mu
         write(6,*) cos(radphi), ' ', a_1
         write(6,*) sin(radphi), ' ',  a_2
!---------Initialize-----------------------------------------------------------------------
	
        if (phi < 100.d0) then
       	 if (Nx > 99) then 
             write(results_name,'(A,f8.6,A,I3,A,F5.2,A)') 'output/l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat' 
			 write (histname,"(A,f8.6,A,I3,A,F5.2,A)") 'output/Hist_l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat'
             write(stats_name,'(A,f8.6,A,I3,A,F5.2,A)') 'stats/l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat'     
             write(cnf_name,'(A,f8.6,A,I3,A,F5.2,A)') 'cnfs/l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat' 
             write(fname,'(A,f8.6,A,I3,A,F5.2,A)') 'cnfs/RND_l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat'             
         else if (Nx > 9 .AND. Nx < 100 ) then         
             write(results_name,'(A,f8.6,A,I2,A,F5.2,A)') 'output/l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat'  
			 write (histname,"(A,f8.6,A,I2,A,F5.2,A)") 'output/Hist_l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat'  
             write(stats_name,'(A,f8.6,A,I2,A,F5.2,A)') 'stats/l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat' 
             write(cnf_name,'(A,f8.6,A,I2,A,F5.2,A)') 'cnfs/l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat' 
             write(fname,'(A,f8.6,A,I2,A,F5.2,A)') 'cnfs/RND_l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat' 
         else 
       	       write(results_name,'(A,f8.6,A,I1,A,F5.2,A)') 'output/l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat' 
              write (histname,"(A,f8.6,A,I1,A,F5.2,A)") 'output/Hist_l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat'
               write(stats_name,'(A,f8.6,A,I1,A,F5.2,A)') 'stats/l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat'          
               write(cnf_name,'(A,f8.6,A,I1,A,F5.2,A)') 'cnfs/l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat' 
                 write(fname,'(A,f8.6,A,I1,A,F5.2,A)') 'cnfs/RND_l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat'                
       	end if     
           else
            if (Nx > 99) then 
             write(results_name,'(A,f8.6,A,I3,A,F6.2,A)') 'output/l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat' 
			 write (histname,"(A,f8.6,A,I3,A,F6.2,A)") 'output/Hist_l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat' 
             write(stats_name,'(A,f8.6,A,I3,A,F6.2,A)') 'stats/l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat' 
             write(cnf_name,'(A,f8.6,A,I3,A,F6.2,A)') 'cnfs/l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat' 
             write(fname,'(A,f8.6,A,I3,A,F6.2,A)') 'cnfs/RND_l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat' 
         else if (Nx > 9 .AND. Nx < 100 ) then         
             write(results_name,'(A,f8.6,A,I2,A,F6.2,A)') 'output/l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat'             
			 write (histname,"(A,f8.6,A,I2,A,F6.2,A)") 'output/Hist_l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat'  
             write(stats_name,'(A,f8.6,A,I2,A,F6.2,A)') 'stats/l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat' 
             write(cnf_name,'(A,f8.6,A,I2,A,F6.2,A)') 'cnfs/l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat' 
             write(fname,'(A,f8.6,A,I2,A,F6.2,A)') 'cnfs/RND_l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat' 
         else 
       	       write(results_name,'(A,f8.6,A,I1,A,F6.2,A)') 'output/l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat' 
              write (histname,"(A,f8.6,A,I1,A,F6.2,A)") 'output/Hist_l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat' 	
              write(stats_name,'(A,f8.6,A,I1,A,F6.2,A)') 'stats/l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat' 
              write(cnf_name,'(A,f8.6,A,I1,A,F6.2,A)') 'cnfs/l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat'    
              write(fname,'(A,f8.6,A,I1,A,F6.2,A)') 'cnfs/RND_l_',l_cx,'_Nx_',Nx,'_phi_',phi,'.dat'              
       	end if  
           end if            

!---------Initialize-----------------------------------------------------------------------
       
     if (config_mode == 1) then      
         CALL READ_CONFIG 
         CALL READ_STATS
         CALL ASSOCIATIONS         
      else      
         CALL INIT_CNF
         CALL INIT_STAT  
         CALL ASSOCIATIONS         
         !----------------------anneal
            	DO	
		
                    if (Nparts == 0) then			
                        CALL CREATE	
                    else if ( rndm() < 0.5d0 ) then                     			
                        CALL CREATE			
                    else    	        	
                        CALL REMOVE	        	 	        	
                    end if
	        
                    i_t = i_t + 1.d0
	        
                    if (i_t > t_step) then
                        exit
                    end if
	        
                END DO
         !------------------------------
      end if    
!---------------Main loop--------------------------------------------------------------------------------------------
	
	DO	
		
		if (Nparts == 0) then			
			CALL CREATE	
		else if ( rndm() < 0.5d0 ) then                     			
		        CALL CREATE			
	        else    	        	
	               	CALL REMOVE	        	 	        	
	        end if

	        CALL MEASURE
	            
	            step = step + 1.d0	          
	            i_p = i_p + 1.d0
	            i_s = i_s + 1.d0	               
                    
              if (i_p >= p_step) then                            !=> write the results, including errors                   
                  i_p = 0.d0                                  
                  CALL PRNT_RESULTS                 
              end if
              
               if (i_s  >= s_step)  then                          !=> save the configuration and statistics            
                   i_s = 0.d0                
                   CALL SAVE_CONFIG 
                   CALL SAVE_STATS                     
               end if                          		
        END DO
		
		 call MPI_FINALIZE(ierr)
	   
    CONTAINS
    
    !---------------------------------------------------------------------------------------------	
	    SUBROUTINE READ_INPUT_FILE 
	    
        OPEN(1, FILE='par')
        READ(1,*) mu_c
        READ(1,*) Nx        
	    READ(1,*) Ny
	    READ(1,*) phi
        READ(1,*) l_cx 
        READ(1,*) Rhc
        READ(1,*) Rsw
        READ(1,*) step_sc
        READ(1,*) t_step                       ! thermolization
        READ(1,*) p_step                       ! print
        READ(1,*) s_step                       ! save
        READ(1,*) config_mode
        read(1,*) reset_code
        
        close(1)
        
        radphi = 3.1415926535d0*phi/180.d0 
        
        L = Rsw*Nx
        lin_bin_size = Rsw                                   ! careful about setting params here
        lin_Nbins = Nx
        Nbins = lin_Nbins*Lin_Nbins
        bin_max = 8
        Nmax = INT(2.6*Nx*Nx)
        
        mu = mu_c - l_cx*sin(radphi)
        Usw = usw_c + l_cx*cos(radphi)
		
		l_cx_check = sqrt((usw - usw_c)*(usw - usw_c) + (mu - mu_c)*(mu - mu_c))
		a_1 = (usw - usw_c)/l_cx
		a_2 = (mu_c - mu)/l_cx
		
        
        ALLOCATE(Rx(1:Nmax))
        ALLOCATE(Ry(1:Nmax))
        ALLOCATE(bin_count(1:Nbins))
        ALLOCATE(bin_parts(1:Nbins,1:bin_max))
        ALLOCATE(neighbors(1:Nbins,1:8))
        
        ALLOCATE(Hist_N(0:Nmax))
	ALLOCATE(Hist_P(0:4*Nmax))
	ALLOCATE(Hist_N_P(0:Nmax,0:4*Nmax))
        
        t_step = t_step*d*L*L
        p_step = p_step*d*L*L
        s_step = s_step*d*L*L
            
        END SUBROUTINE READ_INPUT_FILE
!---------------------------------------------------------------------------------------------------------------------
        SUBROUTINE INIT_STAT       
        ! initialize statistical quantities
        
        Z   = 0.d0 
        
        r_attempt = 0.d0
        c_attempt = 0.d0
        r_accept = 0.d0
        c_accept = 0.d0
        
        Hist_N = 0.d0
	    Hist_P = 0.d0
	    Hist_N_P = 0.d0

        
        Q_Nparts = 0.d0
        Q_Npairs = 0.d0
        Q_Nparts_sq = 0.d0
        Q_Npairs_sq = 0.d0
        Q_Nparts_Npairs = 0.d0
        Q_Nparts_c  = 0.d0
        Q_Nparts_q = 0.d0
        Q_Npairs_c  = 0.d0
        Q_Npairs_q = 0.d0
        Q_Npairs_f  = 0.d0
        Q_Npairs_f = 0.d0 
        
        Q_Npts2_Np = 0.d0
        Q_Npts3_Np = 0.d0
        Q_Npts4_Np = 0.d0
        Q_Np2_Npts = 0.d0
        Q_Np3_Npts = 0.d0
        Q_Np4_Npts = 0.d0
        
        Q_1 = 0.d0
        Q_2 = 0.d0
        Q_3 = 0.d0
        Q_4 = 0.d0
        Q_5 = 0.d0
        Q_6 = 0.d0
        Q_7 = 0.d0
        Q_8 = 0.d0
        Q_9 = 0.d0
        Q_10 = 0.d0
        Q_11 = 0.d0
        Q_12 = 0.d0
        Q_14 = 0.d0
        Q_15 = 0.d0
        Q_16 = 0.d0
        Q_17 = 0.d0
        Q_18 = 0.d0
        Q_19 = 0.d0
        Q_20 = 0.d0
		Q_21 = 0.d0
        Q_22 = 0.d0
		Q_23 = 0.d0
        Q_24 = 0.d0
        Q_25 = 0.d0
		
        
        step = 0.d0
        i_p  = 0.d0
        i_s  = 0.d0
        i_t  = 0.d0
        
        prnt_number = 0

        Q_1_prntouts = 0.d0
        Q_2_prntouts = 0.d0
        Q_3_prntouts = 0.d0 
        Q_4_prntouts = 0.d0
        Q_5_prntouts = 0.d0 
        Q_6_prntouts = 0.d0
        Q_7_prntouts = 0.d0
        Q_10_prntouts = 0.d0
        Q_14_prntouts = 0.d0
        
        END SUBROUTINE INIT_STAT
        !-------------------------------------------------------------------------
        	SUBROUTINE CREATE
        	
        	        x_j = L*rndm()        !  propose to create particle at (x_j,y_j)
			y_j = L*rndm()
			bounce = .FALSE.
			n_j = 0               ! initialize var (new pairs from creation update)
			jj = 0
			cbin = 1 + int(x_j/lin_bin_size) + lin_Nbins*int(y_j/lin_bin_size)
			
			c_attempt = c_attempt + 1.d0
					
		   DO WHILE (jj <= 8)   				! loop through the 9 bins
			
			 if (bounce .EQV. .TRUE.)then
			 	 exit
			 end if                	 
	                         if (jj == 0) then
	                         	 bin = cbin
	                         else
	                         	 bin = neighbors(cbin,jj)
	                         end if	      			
	                         	 p = 1                 !  for looping through existing particles (IN BINS)                         	 
			DO WHILE ( p  <= bin_count(bin))                  ! loop through particles in bin
			                                                  ! calculate interaction number n_j, set acc_r to zero if hard core interaction occurs			
			    dx = abs(x_j - Rx(bin_parts(bin,p)))			    
			    if ( dx > L/2.d0 ) then
			    	     dx = L - dx
			    end if					    
			    dy = abs(y_j - Ry(bin_parts(bin,p)))			    
			    if ( dy > L/2.d0) then
			    	    dy = L - dy
			    end if			    
			    dist = sqrt(dx*dx + dy*dy)			    
			       if ( dist < Rhc ) then
			       	        bounce = .TRUE.      ! hard-core interaction. exit loop
			       	        exit
			       end if			       
			       if ( dist <= Rsw) then
			       	       n_j = n_j + 1
			       end if 	
			       
			    p = p + 1				    
			    
			END DO
			  jj = jj + 1
			END DO
			
			if ( bounce .EQV. .TRUE.) then
				acc_r = 0.d0				
		        else 
		        	acc_r = (L*L/(Nparts +1))*exp(n_j*Usw + mu)						
		        end if
		        
		        if (rndm() < acc_r) then
		        	c_accept = c_accept + 1.d0
		        	Nparts = Nparts + 1
		        	Rx(Nparts) = x_j
		        	Ry(Nparts) = y_j
		        	Npairs = Npairs + n_j	
		        	bin_count(cbin) = bin_count(cbin) + 1
		        	bin_parts(cbin,bin_count(cbin)) = Nparts
		        	
			end if
        	
        	END SUBROUTINE CREATE
       !--------------------------------------------------------------------------
       	       SUBROUTINE REMOVE 
	           
	        	p = RN(Nparts)       ! for this update p is fixed
	        	n_j = 0
	        	jj = 0
	        	rbin = 1 + int(Rx(p)/lin_bin_size) + lin_Nbins*int(Ry(p)/lin_bin_size)	
	        	
	        	r_attempt = r_attempt + 1.d0
	        	
	        	DO WHILE ( jj <= 8 )
	        	
	        	if ( jj == 0) then
	        		
	        		bin = rbin
	                else 
	                	bin = neighbors(rbin,jj)
	                end if
	        		i = 1
	        		
	        	DO WHILE (i <= bin_count(bin))
	      
	        	 if (bin_parts(bin,i) == p) then
	        	 	 tag = i 
	        	 	 i = i + 1
	        	 	 CYCLE      
	        	 end if
	        	 
	        	    dx = abs(Rx(p) - Rx(bin_parts(bin,i)))			    
			    if ( dx > L/2.d0 ) then
			    	     dx = L - dx
			    end if					    
			    dy = abs(Ry(p) - Ry(bin_parts(bin,i)))			    
			    if ( dy > L/2.d0) then
			    	    dy = L - dy
			    end if			    
			    dist = sqrt(dx*dx + dy*dy)
			    	if ( dist <= Rsw) then
			       	       n_j = n_j + 1
			       end if 	
	        	 i = i + 1
	        	 
	        	 END DO
	        	 jj = jj + 1
	        	 END DO
	        	 
	        	 acc_r = (Nparts/(L*L))*exp(-n_j*Usw - mu)
	        	 	        	 
	        	 if (rndm() < acc_r) then
	        	 	r_accept = r_accept + 1.d0
	        	 	Rx(p) = Rx(Nparts)
		        	Ry(p) = Ry(Nparts)
		  ! actually here you need to find the bin where particle labeled by Nparts is and relabel that as p in the bin also!!
		        	bin_sw = 1 + int(Rx(Nparts)/lin_bin_size) + lin_Nbins*int(Ry(Nparts)/lin_bin_size)
		        	DO jjj = 1, bin_count(bin_sw)
		        	    if(bin_parts(bin_sw,jjj) == Nparts) then
		        	    	    bin_parts(bin_sw,jjj) = p
		        	    	    exit
		        	    end if
		        	END DO
		        	Nparts = Nparts - 1
		        	Npairs = Npairs - n_j
		        	bin_parts(rbin,tag) = bin_parts(rbin,bin_count(rbin))
		        	bin_count(rbin) = bin_count(rbin) - 1
			end if
			
       	       END SUBROUTINE REMOVE
        	
       !---------------------------------------------------------------------------
	    SUBROUTINE INIT_CNF
	    
        integer :: ij,kl,jj
        
!       Initializes configuration, sets all the configuration variables to 0  
                  
        Rx = 0.d0
        Ry = 0.d0
        bin_count = 0
        bin_parts = 0
        neighbors = 0
        Nparts = 0
        Npairs = 0
	
!       Initializing random number generator  
                                                                           
         ij = 1802 + 18
         kl = 9373 + 17
         CALL SRANMAR(ij,kl)     
        
        END SUBROUTINE INIT_CNF
!-------------------------------------------------------------------------------------------------------------
	SUBROUTINE ASSOCIATIONS
	
	DO ii = 1, Nbins                      !build left/right associations first!
	
	if (mod(ii,lin_Nbins) == 0) then                           !right column
	        neighbors(ii,1) = ii - 1                           !left neighbor
	        neighbors(ii,2) = ii - lin_Nbins + 1               ! right neighbor
	else if (mod(ii,lin_Nbins) == 1) then                      ! left column
		neighbors(ii,1) = ii + lin_Nbins -1                !left neighbor
		neighbors(ii,2) = ii + 1                           !right neighbor
        else
        	neighbors(ii,1) = ii-1
        	neighbors(ii,2) = ii+1
        end if
        
        END DO
        	
        DO iii = 1, Nbins                             ! now complete the rest of the associations!
        
        if (iii <= lin_Nbins) then                              ! top row
        	 neighbors(iii,3) = iii - lin_Nbins + Nbins               !up neighbor
        	 neighbors(iii,4) = iii + lin_Nbins                       ! down neighbor
        	 neighbors(iii,5) = neighbors(iii - lin_Nbins + Nbins,1)  ! up-left
        	 neighbors(iii,6) = neighbors(iii - lin_Nbins + Nbins,2)  ! up-right
        	 neighbors(iii,7) = neighbors(iii + lin_Nbins,1)          ! down-left
        	 neighbors(iii,8) = neighbors(iii + lin_Nbins,2)          ! down-right
        else if ( iii > lin_Nbins*(lin_Nbins - 1)) then                   ! bottom row
        	neighbors(iii,3) = iii - lin_Nbins
        	neighbors(iii,4) = iii+ lin_Nbins -Nbins
        	neighbors(iii,5) = neighbors(iii - lin_Nbins,1)
        	neighbors(iii,6) = neighbors(iii - lin_Nbins,2)
        	neighbors(iii,7) = neighbors(iii+lin_Nbins - Nbins,1)
        	neighbors(iii,8) = neighbors(iii+lin_Nbins - Nbins,2)
        else
        	neighbors(iii,3) = iii - lin_Nbins
        	neighbors(iii,4) = iii + lin_Nbins
        	neighbors(iii,5) = neighbors(iii - lin_Nbins,1)
        	neighbors(iii,6) = neighbors(iii - lin_Nbins,2)
        	neighbors(iii,7) = neighbors(iii + lin_Nbins,1)
        	neighbors(iii,8) = neighbors(iii + lin_Nbins,2)
        end if
       
        END DO
	
	END SUBROUTINE ASSOCIATIONS
!-------------------------------------------------------------------------------------------------------------
        SUBROUTINE MEASURE
        
        Z = Z + 1.d0
        
        Hist_N(Nparts) = Hist_N(Nparts) + 1.d0       
        Hist_P(Npairs) = Hist_P(Npairs) + 1.d0
        Hist_N_P(Nparts,Npairs) = Hist_N_P(Nparts,Npairs) + 1.d0
                	                
	END SUBROUTINE MEASURE
	!-------------------OUTPUT------------------------------------------
        SUBROUTINE PRNT_RESULTS
        
             integer ::  skip_int,write_code,iii
		double precision :: result_set(0:11),result_set_2(0:5)                                            ! averages
		double precision :: error_set_1(0:7),error_set_2(0:7),error_set_3(0:7),error_set_4(0:7) ,error_set_5(0:7)       ! error analysis
       double precision :: error_set_6(0:7),error_set_7(0:7),error_set_8(0:7),error_set_9(0:7), error_set_10(0:7),error_set_11(0:7)
		double precision :: skip_double
	     integer :: ierr1,icerr1,ierr4,icerr4
		
		prnt_number = prnt_number + 1
      
        Q_Nparts = 0.d0
        Q_Npairs = 0.d0
        Q_Nparts_sq = 0.d0
        Q_Npairs_sq = 0.d0
        Q_Nparts_Npairs = 0.d0
        Q_Npts2_Np = 0.d0
        Q_Npts3_Np = 0.d0
        Q_Npts4_Np = 0.d0
        Q_Np2_Npts = 0.d0
        Q_Np3_Npts = 0.d0
        Q_Np4_Npts = 0.d0        
        Q_Nparts_c  = 0.d0
        Q_Nparts_q = 0.d0
        Q_Npairs_c  = 0.d0
        Q_Npairs_q = 0.d0
        Q_Npairs_f = 0.d0
        Q_Nparts_f = 0.d0
		
	DO i = 0, Nmax
             Q_Nparts = Q_Nparts + i*Hist_N(i)
             Q_Nparts_sq = Q_Nparts_sq + dble(i)*dble(i)*Hist_N(i)
             Q_Nparts_c = Q_Nparts_c + dble(i)*dble(i)*dble(i)*Hist_N(i)
             Q_Nparts_q = Q_Nparts_q + dble(i)*dble(i)*dble(i)*dble(i)*Hist_N(i)
             Q_Nparts_f = Q_Nparts_f + dble(i)*dble(i)*dble(i)*dble(i)*dble(i)*Hist_N(i)
        END DO
               
        DO i = 0, 4*Nmax
        	Q_Npairs = Q_Npairs + i*Hist_P(i)
        	Q_Npairs_sq = Q_Npairs_sq + dble(i)*dble(i)*Hist_P(i)
        	Q_Npairs_c = Q_Npairs_c + dble(i)*dble(i)*dble(i)*Hist_P(i)
        	Q_Npairs_q = Q_Npairs_q + dble(i)*dble(i)*dble(i)*dble(i)*Hist_P(i)
        	Q_Npairs_f = Q_Npairs_f + dble(i)*dble(i)*dble(i)*dble(i)*dble(i)*Hist_P(i)
        END DO		
        
        DO i = 0,Nmax
           DO j = 0,4*Nmax
               Q_Nparts_Npairs = Q_Nparts_Npairs + dble(i)*dble(j)*Hist_N_P(i,j)
               Q_Npts2_Np = Q_Npts2_Np + dble(i)*dble(i)*dble(j)*Hist_N_P(i,j)
               Q_Npts3_Np = Q_Npts3_Np + dble(i)*dble(i)*dble(i)*dble(j)*Hist_N_P(i,j)
               Q_Npts4_Np = Q_Npts4_Np + dble(i)*dble(i)*dble(i)*dble(i)*dble(j)*Hist_N_P(i,j)
               Q_Np2_Npts = Q_Np2_Npts + dble(i)*dble(j)*dble(j)*Hist_N_P(i,j)
               Q_Np3_Npts = Q_Np3_Npts + dble(i)*dble(j)*dble(j)*dble(j)*Hist_N_P(i,j)
               Q_Np4_Npts = Q_Np4_Npts + dble(i)*dble(j)*dble(j)*dble(j)*dble(j)*Hist_N_P(i,j)
               
           END DO
        END DO
		
    Q_1 = Q_Nparts/Z
	Q_2 = Q_Npairs/Z
	Q_3 = Q_Nparts_sq/Z - (Q_Nparts/Z)*Q_Nparts/Z
	Q_4 = Q_Npairs_sq/Z - (Q_Npairs/Z)*Q_Npairs/Z
	Q_5 = Q_Nparts_Npairs/Z - (Q_Nparts/Z)*(Q_Npairs/Z)
  Q_6 = ((Q_3*Q_3)/(Q_Nparts_q/Z - (4.d0*Q_Nparts_c/Z)*(Q_Nparts/Z) + 6.d0*(Q_Nparts_sq/Z)*(Q_Nparts/Z)**2 - 3.d0*(Q_Nparts/Z)**4)) !U4N
  Q_7 = ((Q_4*Q_4)/(Q_Npairs_q/Z - (4.d0*Q_Npairs_c/Z)*(Q_Npairs/Z) + 6.d0*(Q_Npairs_sq/Z)*(Q_Npairs/Z)**2 - 3.d0*(Q_Npairs/Z)**4)) !U4E
	
	Q_8 = Q_Nparts_f/Z -4.d0*(Q_Nparts_c/Z)*(Q_Nparts_sq/Z) - 5.d0*(Q_Nparts_q/Z)*(Q_Nparts/Z)+ &
	12.d0*(Q_Nparts_sq/Z)*(Q_Nparts_sq/Z)*(Q_Nparts/Z) &
	 + 14.d0*(Q_Nparts_c/Z)*(Q_Nparts/Z)*(Q_Nparts/Z)-30.d0*(Q_Nparts_sq/Z)*((Q_Nparts/Z)**3)+12.d0*(Q_Nparts/Z)**5
	Q_9 = Q_Nparts_c/Z - 3.d0*(Q_Nparts_sq/Z)*(Q_Nparts/Z)+2.d0*(Q_Nparts/Z)**3 
	Q_10 = 2.d0*Q_9/Q_3 - Q_8*(Q_6/(Q_3*Q_3)) ! dU_4N/dmu -> idiotic way of writing this lol
	
	Q_11 = Q_Npts4_Np/Z - (Q_Nparts_q/Z)*(Q_Npairs/Z) - 4.d0*(Q_Nparts_c/Z)*(Q_Nparts_Npairs/Z)- &
	         4.d0*(Q_Npts3_Np/Z)*(Q_Nparts/Z)  &
	+ 8.d0*(Q_Nparts/Z)*(Q_Nparts_c/Z)*(Q_Npairs/Z)+12.d0*(Q_Nparts_Npairs/Z)*(Q_Nparts_sq/Z)*(Q_Nparts/Z) & 
	-18.d0*(Q_Nparts_sq/Z)*(Q_Nparts/Z)*(Q_Nparts/Z)*(Q_Npairs/Z) +6.d0*(Q_Npts2_Np/Z)*(Q_Nparts/Z)*(Q_Nparts/Z) &
	-12.d0*(Q_Nparts_Npairs/Z)*(Q_Nparts/Z)**3+12.d0*(Q_Npairs/Z)*(Q_Nparts/Z)**4
 
	Q_12 = Q_Npts2_Np/Z - (Q_Nparts_sq/Z)*(Q_Npairs/Z)-2.d0*(Q_Nparts_Npairs/Z)*(Q_Nparts/Z)+ 2.d0*((Q_Nparts/Z)**2)*(Q_Npairs/Z)
	Q_14 = 2.d0*Q_12/Q_3 - Q_11*(Q_6/(Q_3*Q_3))  ! dU_4N/deps
	
	Q_15 = Q_Npairs_f/Z -4.d0*(Q_Npairs_c/Z)*(Q_Npairs_sq/Z) - 5.d0*(Q_Npairs_q/Z)*(Q_Npairs/Z)+ &
	12.d0*(Q_Npairs_sq/Z)*(Q_Npairs_sq/Z)*(Q_Npairs/Z) &
	 + 14.d0*(Q_Npairs_c/Z)*(Q_Npairs/Z)*(Q_Npairs/Z)-30.d0*(Q_Npairs_sq/Z)*((Q_Npairs/Z)**3)+12.d0*(Q_Npairs/Z)**5
	Q_16 = Q_Npairs_c/Z - 3.d0*(Q_Npairs_sq/Z)*(Q_Npairs/Z)+2.d0*(Q_Npairs/Z)**3
	Q_17 = 2.d0*Q_16/Q_4 - Q_15*(Q_7/(Q_4*Q_4))  !this is dU4E/depsilon
	
	Q_18 = Q_Np4_Npts/Z - (Q_Npairs_q/Z)*(Q_Nparts/Z) - 4.d0*(Q_Npairs_c/Z)*(Q_Nparts_Npairs/Z)-4.d0*(Q_Np3_Npts/Z)*(Q_Npairs/Z) &
	+ 8.d0*(Q_Nparts/Z)*(Q_Npairs_c/Z)*(Q_Npairs/Z)+12.d0*(Q_Nparts_Npairs/Z)*(Q_Npairs_sq/Z)*(Q_Npairs/Z) & 
	-18.d0*(Q_Npairs_sq/Z)*(Q_Npairs/Z)*(Q_Nparts/Z)*(Q_Npairs/Z) +6.d0*(Q_Np2_Npts/Z)*(Q_Npairs/Z)*(Q_Npairs/Z) &
	-12.d0*(Q_Nparts_Npairs/Z)*(Q_Npairs/Z)**3+12.d0*(Q_Nparts/Z)*(Q_Npairs/Z)**4
	Q_19 = Q_Np2_Npts/Z - (Q_Npairs_sq/Z)*(Q_Nparts/Z)-2.d0*(Q_Nparts_Npairs/Z)*(Q_Npairs/Z)+ 2.d0*((Q_Npairs/Z)**2)*(Q_Nparts/Z)
	Q_20 = 2.d0*Q_19/Q_4 - Q_18*(Q_7/(Q_4*Q_4))   ! this is dU4E/dmu
	Q_21 = a_1*a_1*(Q_Npairs_sq/Z - (Q_Npairs/Z)*Q_Npairs/Z) - 2.d0*a_1*a_2*(Q_Nparts_Npairs/Z - (Q_Nparts/Z)*(Q_Npairs/Z)) +   &
	a_2*a_2*(Q_Nparts_sq/Z - (Q_Nparts/Z)*Q_Nparts/Z)
	Q_21 = Q_21/(1.d0*L*L)
	Q_22 = a_1*Q_14 - a_2*Q_10
         
	Q_1_prntouts(prnt_number) = Q_1/(L*L)                                    
        call STAT(Q_1_prntouts,prnt_number,amax,tmax,amin,tmin)            
                
        error_set_1(0) = prnt_number
        error_set_1(1) = Q_1_prntouts(prnt_number)
        error_set_1(2) = amax-amin                                        
        error_set_1(3) = amax                                             
        error_set_1(4) = tmax
        error_set_1(5) = amin
        error_set_1(6) = tmin
       	error_set_1(7) = mu	
       	
        Q_2_prntouts(prnt_number) = Q_2/(L*L)                                    ! These assignments need to be ordered sequentially
        call STAT(Q_2_prntouts,prnt_number,amax,tmax,amin,tmin)                  ! because amax/tmax etc. are global
        
        error_set_2(0) = prnt_number
        error_set_2(1) = Q_2_prntouts(prnt_number)                        
        error_set_2(2) = amax-amin
        error_set_2(3) = amax
        error_set_2(4) = tmax
        error_set_2(5) = amin
        error_set_2(6) = tmin
       	error_set_2(7) = mu	
       	
        Q_3_prntouts(prnt_number) = Q_3/(L*L)                                       
        call STAT(Q_3_prntouts,prnt_number,amax,tmax,amin,tmin)              
        
        error_set_3(0) = prnt_number
        error_set_3(1) = Q_3_prntouts(prnt_number)                          
        error_set_3(2) = amax-amin
        error_set_3(3) = amax
        error_set_3(4) = tmax
        error_set_3(5) = amin
        error_set_3(6) = tmin
       	error_set_3(7) = mu
       	
       	Q_4_prntouts(prnt_number) = Q_4/(L*L)                                       
        call STAT(Q_4_prntouts,prnt_number,amax,tmax,amin,tmin)              
        
        error_set_4(0) = prnt_number
        error_set_4(1) = Q_4_prntouts(prnt_number)                          
        error_set_4(2) = amax-amin
        error_set_4(3) = amax
        error_set_4(4) = tmax
        error_set_4(5) = amin
        error_set_4(6) = tmin
       	error_set_4(7) = mu
       	
       	Q_5_prntouts(prnt_number) = Q_5/(L*L)                                      
        call STAT(Q_5_prntouts,prnt_number,amax,tmax,amin,tmin)              
        
        error_set_5(0) = prnt_number
        error_set_5(1) = Q_5_prntouts(prnt_number)                          
        error_set_5(2) = amax-amin
        error_set_5(3) = amax
        error_set_5(4) = tmax
        error_set_5(5) = amin
        error_set_5(6) = tmin
       	error_set_5(7) = mu
       	
       	Q_6_prntouts(prnt_number) = Q_6                                       
        call STAT(Q_6_prntouts,prnt_number,amax,tmax,amin,tmin)              
        
        error_set_6(0) = prnt_number
        error_set_6(1) = Q_6_prntouts(prnt_number)                          
        error_set_6(2) = amax-amin
        error_set_6(3) = amax
        error_set_6(4) = tmax
        error_set_6(5) = amin
        error_set_6(6) = tmin
       	error_set_6(7) = mu
       	
       	Q_7_prntouts(prnt_number) = Q_7                                       
        call STAT(Q_7_prntouts,prnt_number,amax,tmax,amin,tmin)              
        
        error_set_7(0) = prnt_number
        error_set_7(1) = Q_7_prntouts(prnt_number)                          
        error_set_7(2) = amax-amin
        error_set_7(3) = amax
        error_set_7(4) = tmax
        error_set_7(5) = amin
        error_set_7(6) = tmin
       	error_set_7(7) = mu
       	
       	Q_10_prntouts(prnt_number) = Q_10                                       
        call STAT(Q_10_prntouts,prnt_number,amax,tmax,amin,tmin)              
        
        error_set_8(0) = prnt_number
        error_set_8(1) = Q_10_prntouts(prnt_number)                          
        error_set_8(2) = amax-amin
        error_set_8(3) = amax
        error_set_8(4) = tmax
        error_set_8(5) = amin
        error_set_8(6) = tmin
       	error_set_8(7) = mu
       	
       	Q_14_prntouts(prnt_number) = Q_14                                       
        call STAT(Q_14_prntouts,prnt_number,amax,tmax,amin,tmin)              
        
        error_set_9(0) = prnt_number
        error_set_9(1) = Q_14_prntouts(prnt_number)                          
        error_set_9(2) = amax-amin
        error_set_9(3) = amax
        error_set_9(4) = tmax
        error_set_9(5) = amin
        error_set_9(6) = tmin
       	error_set_9(7) = mu     
  	
       	Q_17_prntouts(prnt_number) = Q_17                                       
        call STAT(Q_17_prntouts,prnt_number,amax,tmax,amin,tmin)              
        
        error_set_10(0) = prnt_number
        error_set_10(1) = Q_17_prntouts(prnt_number)                          
        error_set_10(2) = amax-amin
        error_set_10(3) = amax
        error_set_10(4) = tmax
        error_set_10(5) = amin
        error_set_10(6) = tmin
       	error_set_10(7) = mu  
       	
       	Q_20_prntouts(prnt_number) = Q_20                                       
        call STAT(Q_20_prntouts,prnt_number,amax,tmax,amin,tmin)              
        
        error_set_11(0) = prnt_number
        error_set_11(1) = Q_20_prntouts(prnt_number)                          
        error_set_11(2) = amax-amin
        error_set_11(3) = amax
        error_set_11(4) = tmax
        error_set_11(5) = amin
        error_set_11(6) = tmin
       	error_set_11(7) = mu  
        
            open(8,file = results_name,IOSTAT=ierr1)       
        		
		result_set(0) = mu
		result_set(1) = Q_1/(L*L)   !density
		result_set(2) = Q_2/(L*L)   !pair density
		result_set(3) = Q_3/(L*L)   ! compressibility
		result_set(4) = Q_4/(L*L)   !regular specific heat (factors beta^2)
		result_set(5) = Q_5/(L*L)   ! <N*N_pairs>
		result_set(6) = Q_6         ! U4N
		result_set(7) = Q_7         ! U4E
	        
	    
		     
     
		result_set_2(0) = Q_21        ! specific heat along coexistence line
		result_set_2(1) = Q_22
		result_set_2(2) =  Q_10     ! dU_4N/dmu
        result_set_2(3) =  Q_14     ! dU_4N/deps
		result_set_2(4) =  Q_17     ! dU4E/deps
		result_set_2(5) =  Q_20       ! dU4E/dmu

             write(8,*) ' '
             write(8,'(A,i5,A,F7.5,A,F7.2,A,F18.1)') '  Nx =',Nx,'   Usw = ',Usw,'  Rsw = ',Rsw,'   mc run steps ', step
             write(8,*) ' '
            write(8,'(A)')'     mu     dens      pair_d     CMPR      S.H.     cross     U_4N       U_4E    '
             write(8,*) ' '             
             write(8,'(F10.5,F8.3,F8.3,F9.1,F9.1,F9.1,F10.5,F10.5)') result_set
			    write(8,*) ' '
			  write(8,'(A)')'    dlogZ/dl      dU4/dl     dU_4N/dmu   dU_4N/deps    dU4E/deps   dU4E/dmu '
             write(8,*) ' '             
             write(8,'(F15.5,F15.5,F12.3,F12.3,F9.3,F9.2)') result_set_2

             write(8,*) ' '
    write(8,'(A,F4.2,A,F4.2)') 'score: ', c_accept/c_attempt, '   REMOVE accept/attempt ratio: ', r_accept/r_attempt
             write(8,*) ' '
             write(8,'(A,i12)') 'number of printouts: ',prnt_number
             write(8,*) ' '
             write(8,'(A,F12.5)') 'mu = ', mu
             write(8,*) ' '
             write(8,'(A,F12.6,A,F12.6)') ' <dens> = ', error_set_1(1), '  max-min = ', error_set_1(2)
             write(8,'(F12.6,A,F12.6)') error_set_1(3),' at -> ',error_set_1(4)
             write(8,'(F12.6,A,F12.6)') error_set_1(5),' at -> ',error_set_1(6)
             write(8,*) ' '
             write(8,'(A,F12.6,A,F12.6)') ' <pair_d> = ', error_set_2(1), '  max-min = ', error_set_2(2)
             write(8,'(F12.6,A,F12.6)') error_set_2(3),' at -> ',error_set_2(4)
             write(8,'(F12.6,A,F12.6)') error_set_2(5),' at -> ',error_set_2(6)
             write(8,*) ' '
             write(8,'(A,F12.6,A,F12.6)') ' <Compressibility> = ', error_set_3(1), '  max-min = ', error_set_3(2)
             write(8,'(F12.6,A,F12.6)') error_set_3(3),' at -> ',error_set_3(4)
             write(8,'(F12.6,A,F12.6)') error_set_3(5),' at -> ',error_set_3(6)
             write(8,*) ' '
             write(8,'(A,F12.6,A,F12.6)') " <heat cap.> = ", error_set_4(1), '  max-min = ', error_set_4(2)
             write(8,'(F12.6,A,F12.6)') error_set_4(3),' at -> ',error_set_4(4)
             write(8,'(F12.6,A,F12.6)') error_set_4(5),' at -> ',error_set_4(6)
             write(8,*) ' '
             write(8,'(A,F12.6,A,F12.6)') ' <cross> = ', error_set_5(1), '  max-min = ', error_set_5(2)
             write(8,'(F12.6,A,F12.6)') error_set_5(3),' at -> ',error_set_5(4)
             write(8,'(F12.6,A,F12.6)') error_set_5(5),' at -> ',error_set_5(6)
             write(8,*) ' '
             write(8,'(A,F12.6,A,F12.6)') ' <U_4N > = ', error_set_6(1), '  max-min = ', error_set_6(2)
             write(8,'(F12.6,A,F12.6)') error_set_6(3),' at -> ',error_set_6(4)
             write(8,'(F12.6,A,F12.6)') error_set_6(5),' at -> ',error_set_6(6)
             write(8,*) ' '
             write(8,'(A,F12.6,A,F12.6)') ' <U_4E> = ', error_set_7(1), '  max-min = ', error_set_7(2)
             write(8,'(F12.6,A,F12.6)') error_set_7(3),' at -> ',error_set_7(4)
             write(8,'(F12.6,A,F12.6)') error_set_7(5),' at -> ',error_set_7(6)
             write(8,*) ' '
             write(8,'(A,F12.6,A,F12.6)') ' <dU_4N/dmu> = ', error_set_8(1), '  max-min = ', error_set_8(2)
             write(8,'(F12.6,A,F12.6)') error_set_8(3),' at -> ',error_set_8(4)
             write(8,'(F12.6,A,F12.6)') error_set_8(5),' at -> ',error_set_8(6)
             write(8,*) ' '
             write(8,'(A,F15.5,A,F15.5)') ' <dU_4N/deps> = ', error_set_9(1), '  max-min = ', error_set_9(2)
             write(8,'(F12.6,A,F12.6)') error_set_9(3),' at -> ',error_set_9(4)
             write(8,'(F12.6,A,F12.6)') error_set_9(5),' at -> ',error_set_9(6)
             write(8,*) ' '
             write(8,'(A,F15.5,A,F15.5)') ' <dU_4E/deps> = ', error_set_10(1), '  max-min = ', error_set_10(2)
             write(8,'(F12.6,A,F12.6)') error_set_10(3),' at -> ',error_set_10(4)
             write(8,'(F12.6,A,F12.6)') error_set_10(5),' at -> ',error_set_10(6)
             write(8,*) ' '
             write(8,'(A,F15.5,A,F15.5)') ' <dU_4E/dmu> = ', error_set_11(1), '  max-min = ', error_set_11(2)
             write(8,'(F12.6,A,F12.6)') error_set_11(3),' at -> ',error_set_11(4)
             write(8,'(F12.6,A,F12.6)') error_set_11(5),' at -> ',error_set_11(6)
             write(8,*) ' '
                                 
         close(8,IOSTAT=icerr1) 
		 
         
        OPEN(32, FILE='par',IOSTAT=ierr4)
        
        READ(32,*) skip_int        
	    READ(32,*) skip_int
	    READ(32,*) skip_double
        READ(32,*) skip_double      
        READ(32,*) skip_double
        READ(32,*) skip_double
        READ(32,*) skip_double
        READ(32,*) skip_double                      
        READ(32,*) skip_double                       
        READ(32,*) skip_double                       
        READ(32,*) skip_int
        READ(32,*) reset_code   
        READ(32,*) write_code ! for histograms
        
        CLOSE(32,IOSTAT=icerr4)
         
         if (reset_code == Nx) then
         	 CALL INIT_STAT
         end if
         
         if (write_code == Nx) then
            open(56,file = histname)
               DO iii = 0, Nmax
                   write(56,*) iii, Hist_N(iii)/Z
              END DO
                close(56)
         end if
         	          
	    END SUBROUTINE PRNT_RESULTS
!-----------------stat-------------------------------------------------------------------------------------------
      SUBROUTINE STAT(prnt_storage,n,amax,tmax,amin,tmin) 
!     Analyzing 3/4 print-out
	      
      double precision :: prnt_storage, amax, tmax, amin, tmin, aa
      integer n, i
      dimension prnt_storage(n)

      amax=-1.d200; amin=1.d200

      DO i=n/4+1, n;  aa=prnt_storage(i)
         if (aa > amax) then
            amax=aa; tmax=i
         end if
         if (aa < amin) then
            amin=aa; tmin=i
         end if
      END DO

      tmax=tmax/n; tmin=tmin/n
      END SUBROUTINE STAT	    
!-------------------------------------------------------------------------------------------------------------
	SUBROUTINE READ_CONFIG
	
	 	OPEN(48, file=cnf_name)
	 	
	    READ(48,*)Rx
	 	READ(48,*)Ry
	 	READ(48,*)Nparts
	 	READ(48,*)Npairs
	 	READ(48,*)bin_count
	 	READ(48,*)bin_parts
	 	READ(48,*)neighbors
	 	
	 	CLOSE(48)
	 	
	 	CALL read_rndm(fname)
	
	END SUBROUTINE READ_CONFIG
!-------------------------------------------------------------------------------------------------------------
	SUBROUTINE READ_STATS
	
		OPEN(40,file=stats_name)
		
		READ(40,*)Z
		READ(40,*)Hist_N
		READ(40,*)Hist_P
		READ(40,*)Hist_N_P
		READ(40,*)step
		READ(40,*)i_p
		READ(40,*)i_s
		READ(40,*)i_t
		READ(40,*)prnt_number
		READ(40,*)Q_1_prntouts
		READ(40,*)Q_2_prntouts
		READ(40,*)Q_3_prntouts
		READ(40,*)Q_4_prntouts
		READ(40,*)Q_5_prntouts
		READ(40,*)Q_6_prntouts
		READ(40,*)Q_7_prntouts
		READ(40,*)Q_10_prntouts
		READ(40,*)Q_14_prntouts 		
		
		CLOSE(40)
				
	END SUBROUTINE READ_STATS
!-------------------------------------------------------------------------------------------------------------
	SUBROUTINE SAVE_CONFIG
	
	        OPEN(48, file=cnf_name)
	 	
	 	WRITE(48,*)Rx
	 	WRITE(48,*)Ry
	 	WRITE(48,*)Nparts
	 	WRITE(48,*)Npairs
	 	WRITE(48,*)bin_count
	 	WRITE(48,*)bin_parts
	 	WRITE(48,*)neighbors
	 	
	 	CLOSE(48)
	 	
	 	CALL save_rndm(fname) 
	
	END SUBROUTINE SAVE_CONFIG
!-------------------------------------------------------------------------------------------------------------
	SUBROUTINE SAVE_STATS
	
		OPEN(40,file=stats_name)
		
		WRITE(40,*)Z
		WRITE(40,*)Hist_N
		WRITE(40,*)Hist_P
		WRITE(40,*)Hist_N_P
		WRITE(40,*)step
		WRITE(40,*)i_p
		WRITE(40,*)i_s
		WRITE(40,*)i_t
		WRITE(40,*)prnt_number
		WRITE(40,*)Q_1_prntouts
		WRITE(40,*)Q_2_prntouts
		WRITE(40,*)Q_3_prntouts
		WRITE(40,*)Q_4_prntouts
		WRITE(40,*)Q_5_prntouts
		WRITE(40,*)Q_6_prntouts
		WRITE(40,*)Q_7_prntouts
		WRITE(40,*)Q_10_prntouts
		WRITE(40,*)Q_14_prntouts 
		
		CLOSE(40)
	
	END SUBROUTINE SAVE_STATS
!-------------------------------------------------------------------------------------------------------------
	
	
	END PROGRAM MAIN