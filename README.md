#============================================================
#      RMSD ANALYSIS & MOLECULAR ALIGNMENT TOOL
#============================================================

This tool calculates the Root Mean Square Deviation (RMSD) between 
molecular structures and generates aligned coordinates and overlay 
images for visual comparison.

------------------------------------------------------------
1. PREREQUISITES & INSTALLATION
------------------------------------------------------------
Ensure you have Python installed. Then, install the required 
libraries using pip:

    pip install numpy vtk Pillow

* numpy: For numerical and matrix calculations.
* vtk: For high-quality 3D molecule visualization.
* Pillow (PIL): For final image processing and text overlays.

------------------------------------------------------------
2. PROJECT STRUCTURE
------------------------------------------------------------
Your workspace should be organized as follows:

Codigos/
├── ALL_CODES.py           # Main script
├── READ_ME.txt            # This documentation
├── INPUT/                 # Data folder
│   ├── referencias/       # (Option 1) Place your .xyz references here
│   └── estruturas/        # (Option 1) Place your .xyz structures here
└── OUTPUT/                # Results will be generated here
    ├── RMSD.txt           # Main results report
    ├── Alinhamentos/      # Optional aligned .xyz files
    └── Imagens/           # Optional overlay .png images

------------------------------------------------------------
3. FILENAME CONVENTION
------------------------------------------------------------
The script uses a strict matching logic. A structure file will only 
be compared to a reference if its name starts with the reference 
name followed by a separator (_ or .) or the end of the string.

Example:
Reference: ORIGII.xyz
Matched Structures: ORIGII_PM6.xyz, ORIGII.01.xyz
Ignored Structures: ORIGII01.xyz (mismatch due to missing separator)

------------------------------------------------------------
4. HOW TO USE (Step-by-Step)
------------------------------------------------------------

STEP 1: Run the script
    python ALL_CODES.py

STEP 2: Choose Input Method
    Select '1' to use the INPUT/ subfolders or '2' for files in 
    the root. The script will pause to let you move your files.

STEP 3: Review Statistics
    The script will list all identified pairs. Verify that your 
    molecules were correctly matched.

STEP 4: Select Outputs
    * RMSD.txt is always generated.
    * You can choose to save Aligned XYZ files.
    * You can choose to generate Overlay Images.
    * For both options, you can select specific files from the 
      list to save space and time.

STEP 5: Collect Results
    Find your data in the OUTPUT/ folder. Aligned molecules are 
    centered on the reference, and images show the Reference 
    (Red) vs. the Structure (Green).

------------------------------------------------------------
5. CONTACT & SUPPORT
------------------------------------------------------------
For issues or feedback, please check your input file formats 
(standard .xyz) and ensure all molecules in a comparison pair 
have the SAME number of atoms.

┏┓    ┓   
┃┃┏┓┏┓┃┏┓
┣┛┗┛┣┛┗┗  
