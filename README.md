Replication Package: "Learning about the Long Run" (Farmer, Nakamura, & Steinsson, 2024)

HOW TO RUN THE REPLICATION (THE MASTER CONTROL PANEL)

To run this replication, open `run_replication.m`. At the very top of the script, you will find the "Master Control Panel" with three logical toggles. 

Because GitHub file-size limits prevented me from uploading the multi-gigabyte generated output files, you will need to run the mathematical estimations from scratch. Make sure that for the very first run the following settings are chosen:

1. Set `run_full_tbill_estimation = true;`
2. Set `run_full_gdp_estimation = true;`
3. Set `run_full_mc_estimation = false;` (The authors included all necessary Monte Carlo .mat files in the original package, so this does not need to be re-run).

Hit Run. 
*NOTE ON RUNTIMES:* On an Apple M2 Pro MacBook Pro, the T-bill MCMC estimation took approximately 5 hours, and the GDP estimation took approximately 3 hours. I recommend letting this run overnight. Once the run is complete and the large .mat files are generated, you can change the toggles back to `false` to instantly regenerate all figures and tables without re-running the math.

DEVIATIONS & HARDWARE LIMITATIONS

1. Hardware Limit (Figure 9): The code required to generate Figure 9 (`tblPriorEstimation_Break.m`) attempts to pre-allocate an 18.4GB matrix (`NaN(275,120,75000)`). Because this exceeds physical laptop RAM capacities, I have explicitly bypassed this specific script in the wrapper to prevent memory overflow crashes. 
2. Typo Fixes: The wrapper script includes automated string-replacements to fix a hardcoded local desktop path in the authors' `figure_5.m` file, to correct variable naming typos in `tblFigs.m` (`gam1Params` to `gamParams`), and to fix a hardcoded file load error in the Monte Carlo section.

OUTPUT MAPPING DIRECTORY

The authors' code prints raw matrices directly to the command window and generates several undocumented diagnostic figures. Below is the map connecting the code outputs to the published paper:

TABLES (Found in Tables_Output.txt):
- Tables 1 & 2: Generated during the "Data" folder execution.
- Table 5 (CBO Results): Generated at the end of the "Data" folder execution.
- Tables 3 & 4: Generated during the "T-bill" folder execution.
- Table 5 (UC Model Results): Generated during the "GDP" folder execution.
- Tables 6 & 7: Generated during the "Monte Carlo" execution (Outputs sequentially for Downward, Unbiased, and Upward priors).

MAIN FIGURES (Saved as PNGs in the root folder):
- Data_Fig_1, 2, 3, 5 correspond to Figures 1, 2, 3, 5.
- Tbill_Fig_4, 6, 7, 8 correspond to Figures 4, 6, 7, 8.
- GDP_Fig_10, 11, 12 correspond to Figures 10, 11, 12.
- MC_Fig_13, 14, 15 correspond to Figures 13, 14, 15.
*(Any extra generated figures are the authors' diagnostic plots not included in the main text).*

## Update: Epoch Sensitivity Analysis (Extension)
In addition to the base replication, I have conducted an extension analysis testing the parameter stability of the unobserved components model across different historical regimes. The following files have been added to the repository:

**1. Estimation Scripts:**
* `gdp_GreatMod.m`: Runs the MCMC Gibbs sampler on a restricted dataset covering the Great Moderation (1984-2007).
* `gdp_ModernEra.m`: Runs the estimation for the Modern Era (1984-2019).

**2. Analysis Scripts:**
* `Epoch_Analysis_Plot.m`: Generates a comparative Kernel Density plot of the posterior $\gamma$ distributions.
* `Epoch_Analysis_Table.m`: Generates the statistical summary table for the posterior estimates.

**3. Output:**
* `Epoch_Analysis_Gamma.png`: The plot produced by `Epoch_Analysis_Plot.m`.

*Note: The generated `.mat` files from the extension runs are excluded due to GitHub size constraints. Please run the two `gdp_*.m` scripts to generate the posterior draws locally before running the plot and table scripts. The main replication code needs to be run first before running this extension.* 
