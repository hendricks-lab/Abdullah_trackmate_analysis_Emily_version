# Abdullah-codes, modified by Emily
Before you download any of these codes please make sure to download the General codes folders, there are 3 that must be downloaded to a single folder. Without that file many of these codes wouldn't work. Accordingly, when you download the general codes folder please make sure to add it to path in your code via addpath('XYZ\X\Y\Z'). Also, the directories need to change for all your codes. All of the codes work as k_choose == 1, 2, 3 etc where with in each k_choose you can select multiple folders and cd to them thus each k_choose is a different condition.
 
- Trackmate motility analysis (Live cell imaging)
    - STEP 1: Run Analyze_Trackmate_motility_v2. This code takes all the .csv data from trackmate (x,y positions) and computes the positional data. Before you run the code make sure to add cell center for directionality in columns 22 and 23 and change DT (Exposure time). The x,y,1D position, and origin points (cell center in x,y) are saved in a new .mat file which is then read in subsequent codes. The output file is called "1Spots....mat" and I saved it in "Processed Tracks" folder. 
    - STEP 2: Run the "Analyze_dispersion_rg_MSD_per_cell.m". To run this code cd to the output .mat files from STEP 1 ("1Spot....mat"). There's a few points to keep in mind. You have to create two folders - Rg, MSD, and dirbias. The save_dir_Rg, save_dir_MSD, and save_dirdirbias directories will save the output of this code there.
    - STEP 3: To analyze the results from STEP 2, run the "Plot_Rg_and_MSD_analysis_per_cell.m" file. Make sure to add all the paths that are on the code to your specific directory. pth1, pth2, and pth3  cd to the Rg and MSD folders so add those paths.


- Additional code to make position vs time plots:

	- Run Plot_xk_yk_per_cell.m. This code runs the output .mat files from STEP 1. To change the colours for the trajectories, check lines 129-132. To change the random arrangement along the time axis, check line 94. To increase spacing between position lines, check line 93.