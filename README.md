#      RMSD ANALYSIS & MOLECULAR ALIGNMENT TOOL


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

# Codigos/

    ├── ALL_CODES.py       <--  Main script
    ├── READ_ME.txt        <-- This documentation
    ├── INPUT/             <-- Data folder
    │   ├── referencias/   <-- (Option 1) Place your .xyz references here
    │   └── estruturas/    <-- (Option 1) Place your .xyz structures here
    └── OUTPUT/            <-- Results will be generated here
        ├── RMSD.txt       <-- Main results report
        ├── Alinhamentos/  <-- Optional aligned .xyz files
        └── Imagens/       <-- Optional overlay .png images

------------------------------------------------------------
3. FILENAME CONVENTION
------------------------------------------------------------
The script uses a strict matching logic. A structure file will only 
be compared to a reference if its name starts with the reference 
name followed by a separator (_ or .) or the end of the string.

Example:
Reference: [name].xyz
Matched Structures: [name]_PM6.xyz, [name].01.xyz
Ignored Structures: [name]01.xyz (mismatch due to missing separator)

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


#      RMSD 分析与分子对齐工具


本工具计算分子结构之间的均方根偏差（RMSD），并生成对齐后的坐标及叠加图像，用于视觉比较。

------------------------------------------------------------
1. 环境要求与安装
------------------------------------------------------------
请确保已安装 Python，然后使用 pip 安装所需库：

    pip install numpy vtk Pillow

* numpy：用于数值与矩阵计算。
* vtk：用于高质量 3D 分子可视化。
* Pillow（PIL）：用于最终图像处理及文字叠加。

------------------------------------------------------------
2. 项目结构
------------------------------------------------------------
您的工作目录应组织如下：

# Codigos/

    ├── ALL_CODES.py       <-- 主脚本
    ├── READ_ME.txt        <-- 本文档
    ├── INPUT/             <-- 数据文件夹
    │   ├── referencias/   <-- （方式1）将参考 .xyz 文件放入此处
    │   └── estruturas/    <-- （方式1）将待比对结构 .xyz 文件放入此处
    └── OUTPUT/            <-- 结果将生成在此处
        ├── RMSD.txt       <-- 主要结果报告
        ├── Alinhamentos/  <-- 可选：对齐后的 .xyz 文件
        └── Imagens/       <-- 可选：叠加图像 .png 文件

------------------------------------------------------------
3. 文件名匹配规则
------------------------------------------------------------
脚本使用严格的匹配逻辑：只有当结构文件的名称以参考文件名开头，并且紧跟着分隔符（_ 或 .）或直接结束时，才会将该结构与参考文件进行比对。

示例：
参考文件：[name].xyz
匹配的结构：[name]_PM6.xyz、[name].01.xyz
忽略的结构：[name]01.xyz（因缺少分隔符而不匹配）

------------------------------------------------------------
4. 使用方法（逐步说明）
------------------------------------------------------------

步骤 1：运行脚本

    python ALL_CODES.py

步骤 2：选择输入方式

    选择 '1' 使用 INPUT/ 子文件夹，或选择 '2' 使用根目录下的文件。
    脚本会暂停，让您移动文件。

步骤 3：查看统计信息

    脚本会列出所有识别出的配对。请确认分子匹配正确。

步骤 4：选择输出内容

    * RMSD.txt 始终生成。
    * 可选择保存对齐后的 XYZ 文件。
    * 可选择生成叠加图像。
    * 对于上述两种可选输出，您可以从列表中指定特定文件，以节省空间和时间。

步骤 5：获取结果

    结果位于 OUTPUT/ 文件夹中。对齐后的分子以参考分子为中心进行叠合，图像中参考分子显示为红色，待比对结构显示为绿色。

------------------------------------------------------------
5. 联系与支持
------------------------------------------------------------
如有问题或反馈，请检查输入文件格式（标准 .xyz），并确保每一对比较的分子具有相同的原子数。

    ┏┓    ┓   
    ┃┃┏┓┏┓┃┏┓
    ┣┛┗┛┣┛┗┗  
