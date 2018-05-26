* Organization
  - There are two directories, each correspond to a subsection in the experiments section of the paper: mixef (linear mixed effects modeling) and ml (MovieLens data).
  - Each directory has four sub directories: code, data, qsub, and result.
  - Directory 'code' has files of all the R source code that was used in the analysis. 
  - Directory 'data' has (if any) simulated/real data that was used in the analysis. This directory may be empty or absent.
  - Directory 'qsub' has SGE files (.q) that were used to submit jobs on a SGE cluster. 
  - Directory 'result' has a sub directory 'img' and stores the result (if any) produced in the analysis. This directory may be empty or absent.

* Files
  - 'simulate_data.R' contains the code to simulate the data in linear mixed-effects modeling. 
  - 'partition_data.R' contains the code to partition the data into smaller subsets that are stored across workers.   
  - 'dem_estep_sync.R' contains code for performing the E step of DEM on the K workers in parallel.
  - 'dem_mstep_sync.R' contains code for performing the M step of DEM on the managers once they have received results from \gamma-fraction of the workers.  
  - '(mixef|ml)_dem_mpi.R' contains code for fitting linear mixed-effects model using DEM and MPI in simulations and real data analysis.    
  - '(mixef|ml)_iem_mpi.R' contains code for fitting linear mixed-effects model using IEM and MPI in simulations and real data analysis.      
  - 'analyze_result.R' contains the code for analyzing the results of DEM and competing methods and making plots/tables.
  - 'lmer.R' contains the code for fitting linear mixed-effects using lme4 R package.
  - 'vandyk00.R' contains the code for fitting linear mixed-effects using ECME0 algorithm of van Dyk (2000).  
  - '(mixef|ml)_submit_dem.R' contains the  code for the R code for submitting a job on the cluster. The files in 'qsub' directory use this file for running simulations.  

* Citation
  If you use the code, then please cite the following two papers:
  - Van Dyk, D. A. (2000). Fitting mixed-effects models using efficient EM-type algorithms. Journal of Computational and Graphical Statistics, 9(1), 78-98.
  - Srivastava, S., DePalma, G.R. and Liu, C. (2018). An Asynchronous Distributed Expectation-Maximization Algorithm For Massive Data: The DEM Algorithm. Revision submitted to JCGS.
   
* Contact
  Please email Sanvesh Srivastava (<sanvesh-srivastava@uiowa.edu>) if you have any questions related to the code.

* Acknowledgment
  - Some code for linear mixed effects modeling has been borrowed from Patrick O. Perry (<http://ptrckprry.com/code/>).
  - Please email us if you think that we have missed citations to your paper/work. 
