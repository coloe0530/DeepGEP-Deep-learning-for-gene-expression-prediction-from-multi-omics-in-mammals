Software Usage Instructions
The software includes one shell script (preprocess.sh) and one Python script (kf_model_with0.py). Users can directly modify parameters and file paths within the scripts.

Usage instructions:

Before running the scripts, make sure Python and the required modules are properly installed, and that the scripts are executed in a compatible system environment. You will also need to prepare chromatin accessibility data (ATAC-seq) and four types of histone modification data (ChIP-seq for H3K4me3, H3K4me1, H3K27ac, and H3K27me3).

First, run the preprocess.sh script. Then, use the output file to run the Python script kf_model_with0.py. When running the scripts, only essential parameters need to be specified (as shown in Figure 2 and Figure 3); just input the corresponding information such as file paths and names. When the Python script completes successfully, the message "Run successfully!" will be displayed.
