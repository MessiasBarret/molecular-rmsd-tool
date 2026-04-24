============================================================
#      RMSD ANALYSIS & MOLECULAR ALIGNMENT TOOL
============================================================

This tool calculates the Root Mean Square Deviation (RMSD) between 
molecular structures, offering three different mathematical 
approaches, and generates aligned coordinates and overlay images.

------------------------------------------------------------
1. PREREQUISITES & INSTALLATION
------------------------------------------------------------
Ensure you have Python installed. Then, install the required 
libraries using pip:

    pip install numpy vtk Pillow

* numpy: For numerical and matrix calculations.
* vtk: For high-quality 3D molecule visualization.
* Pillow (PIL): For image post-processing and text overlays.

------------------------------------------------------------
2. PROJECT STRUCTURE
------------------------------------------------------------
The script adapts to two main workflows:

Workflow 1 (Subfolders):

    Codigos/
    ├── ALL_CODES.py
    ├── INPUT/
    │   ├── referencias/   # Place .xyz references here
    │   └── estruturas/    # Place .xyz structures here
    └── OUTPUT/            # Results (RMSD.txt, Alinhamentos/, Imagens/)

Workflow 2 (Root/Cross-match):

    Codigos/
    ├── ALL_CODES.py
    ├── *.ref              # Reference files in root or INPUT/
    ├── *.xyz              # Structure files in root or INPUT/
    └── OUTPUT/            # Results

------------------------------------------------------------
3. MATCHING LOGIC
------------------------------------------------------------
- OPTION 1 (Prefix): Matches files by name prefix.
  Ref: ORIGII.xyz matches ORIGII_01.xyz, ORIGII.test.xyz.
- OPTION 2 (Cross-match): Compares EVERY .ref file against 
  EVERY .xyz file found in the directories.

------------------------------------------------------------
4. TECHNICAL DETAILS
------------------------------------------------------------
- HYDROGENS: The script AUTOMATICALLY IGNORES Hydrogen (H) 
  atoms to focus on the molecular skeleton.
- RMSD TYPES:
    1. Normal: Simple RMSD without alignment.
    2. Kabsch: RMSD after ideal rotation (centered).
    3. Fitted: RMSD after translation + rotation optimization.
- ALIGNMENT: Images show Reference (Red) vs. Structure (Green).

------------------------------------------------------------
5. HOW TO USE (Two Modes)
------------------------------------------------------------

A. INTERACTIVE MODE (Prompt)
   Simply run: python ALL_CODES.py
   Follow the step-by-step questions in the terminal.

B. KEYWORD MODE (Fast Configuration)
   Run the script and enter all parameters in one line:
    
    Format: [FOLDER] [SAVE_XYZ] [SEL_XYZ] [SAVE_PNG] [SEL_PNG]
   
       Example: 1 y a y a
       (Subfolders, save All XYZ, generate All images)

       Example: 1 n y 1,2,3
       (Subfolders, no XYZ, generate images only for pairs 1, 2, and 3)

------------------------------------------------------------
6. CONTACT & SUPPORT
------------------------------------------------------------
Ensure all molecules in a comparison pair have the SAME number 
of non-hydrogen atoms.

____________________________________________________________

【以下内容由AI生成】

============================================================
# RMSD 分析 & 分子比对工具
============================================================

本工具用于计算分子结构之间的均方根偏差（RMSD），提供三种不同的数学方法，并生成比对后的坐标文件及叠加图像。

环境要求与安装

请确保已安装 Python，然后使用 pip 安装所需库：

pip install numpy vtk Pillow

numpy：用于数值计算与矩阵运算。

vtk：用于高质量 3D 分子可视化。

Pillow（PIL）：用于图像后处理及文字叠加。

项目结构

脚本支持两种主要工作流程：

工作流程 1（子文件夹模式）：

    Codigos/
    ├── ALL_CODES.py
    ├── INPUT/
    │ ├── referencias/ # 放置 .xyz 参考文件
    │ └── estruturas/ # 放置 .xyz 结构文件
    └── OUTPUT/ # 输出结果（RMSD.txt、Alinhamentos/、Imagens/）

    工作流程 2（根目录/交叉匹配模式）：
    Codigos/
    ├── ALL_CODES.py
    ├── *.ref # 参考文件，可位于根目录或 INPUT/
    ├── *.xyz # 结构文件，可位于根目录或 INPUT/
    └── OUTPUT/ # 输出结果

匹配逻辑

方式 1（前缀匹配）：按文件名前缀进行匹配。
例如：参考文件 ORIGII.xyz 将匹配 ORIGII_01.xyz、ORIGII.test.xyz。

方式 2（交叉匹配）：将目录中找到的每一个 .ref 文件与每一个 .xyz 文件进行比对。

技术细节

氢原子：脚本自动忽略氢原子（H），聚焦于分子骨架。

RMSD 类型：

普通 RMSD：不对齐，直接计算。
Kabsch RMSD：经理想旋转（中心对齐）后计算。
拟合 RMSD：经平移 + 旋转优化后计算。
图像比对：参考结构显示为红色，待比对结构显示为绿色。

使用方法（两种模式）

A. 交互模式（命令行提示）
直接运行：python ALL_CODES.py
按照终端中的逐步提示进行操作。

B. 关键词模式（快速配置）
运行脚本并在同一行中输入所有参数：

    格式：[FOLDER] [SAVE_XYZ] [SEL_XYZ] [SAVE_PNG] [SEL_PNG]

    示例：1 y a y a
    （使用子文件夹，保存全部 XYZ 文件，生成全部图像）

    示例：1 n y 1,2,3
    （使用子文件夹，不保存 XYZ，只为第 1、2、3 对生成图像）

联系与支持

请确保参与比对的每一对分子具有相同数量的非氢原子。

    ┏┓    ┓
    ┃┃┏┓┏┓┃┏┓
    ┣┛┗┛┣┛┗┗
