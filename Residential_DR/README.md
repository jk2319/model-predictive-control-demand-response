# Optimization of basic residential building operation with Demand Response

The solution of the problem requires IPOPT, its Matlab interface and the Imperial College London Optimal Control Software (ICLOCS)

http://www.ee.ic.ac.uk/ICLOCS/

## Main scripts 

Input_data.m  - Define problem data and read data from the excel file data_problem.xlsx. All the data are saved in the file `case_study_sim.mat`. 
The stored data in `case_study­_sim.mat` are uploaded and used by the optimisation model.

main_res_mpc.m - It defines and solves the 					    	      optimisation problem. 

PlotClosed_loop_solution.m - It plots the closed loop solution

PlotOpen_loop_solution.m - It plots the open loop solution of 					the last problem solved				
