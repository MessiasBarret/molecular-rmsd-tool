"""  
┏┓    ┓  
┃┃┏┓┏┓┃┏┓
┣┛┗┛┣┛┗┗ 
    ┛ 
"""


import vtk
import glob
import os
import re
import numpy as np
from PIL import Image, ImageChops, ImageFont, ImageDraw

####################################################
# BLOCO 1 - CÁLCULOS MATEMÁTICOS E GEOMÉTRICOS
####################################################

def centroid(X):
    """Calcula o centro geométrico (centróide) de um conjunto de pontos."""
    if len(X) == 0: return np.array([0,0,0])
    return np.mean(X, axis=0)

def rmsd(P, Q):
    """Calcula o RMSD simples entre duas matrizes de mesma dimensão."""
    if P.shape != Q.shape: return 999.999
    diff = P - Q
    return np.sqrt((diff * diff).sum() / P.shape[0])

def kabsch(P, Q):
    """Implementa o Algoritmo de Kabsch para encontrar a matriz de rotação ideal."""
    C = np.dot(np.transpose(P), Q)
    V, S, W = np.linalg.svd(C)
    # Garante que a matriz de rotação não resulte em uma reflexão (mão trocada)
    if (np.linalg.det(V) * np.linalg.det(W)) < 0.0:
        V[:, -1] = -V[:, -1]
    return np.dot(V, W)

def rotate(P, Q):
    """Rotaciona a matriz P para se alinhar da melhor forma com a matriz Q."""
    U = kabsch(P, Q)
    return np.dot(P, U)

def kabsch_rmsd(P, Q):
    """Calcula o RMSD após a rotação ideal de Kabsch (requer centralização)."""
    P_c = P - centroid(P)
    Q_c = Q - centroid(Q)
    P_rot = rotate(P_c, Q_c)
    return rmsd(P_rot, Q_c)

def fit(P, Q, max_iter=1000):
    """Otimização de translação para encontrar o melhor encaixe antes da rotação final."""
    step_size = np.max(P, axis=0) - np.min(P, axis=0)
    if (step_size == 0).all(): step_size = np.array([1.0, 1.0, 1.0])
    
    threshold = step_size * 1e-9
    rmsd_best = kabsch_rmsd(P, Q)
    n_iter = 0
    P_current = P.copy()
    
    while n_iter <= max_iter:
        n_iter += 1
        improved = False
        for i in range(3):
            temp = np.zeros(3)
            temp[i] = step_size[i]
            
            # Testa movimento na direção positiva
            rmsd_new = kabsch_rmsd(P_current + temp, Q)
            if rmsd_new < rmsd_best:
                rmsd_best = rmsd_new
                P_current[:, i] += step_size[i]
                improved = True
            else:
                # Testa movimento na direção negativa
                rmsd_new = kabsch_rmsd(P_current - temp, Q)
                if rmsd_new < rmsd_best:
                    rmsd_best = rmsd_new
                    P_current[:, i] -= step_size[i]
                    improved = True
                else:
                    # Reduz o passo se nenhum movimento ajudou
                    step_size[i] /= 2
                    
        if (step_size <= threshold).all(): break
            
    if n_iter > max_iter: return None
    return rmsd_best, P_current

def get_coordinates(filename):
    """Lê um arquivo .xyz e extrai as coordenadas (ignora Hidrogênios)."""
    if not os.path.exists(filename): return None, None
    with open(filename, "r") as f:
        try:
            line = f.readline().strip()
            if not line: return None, None
            n_atoms = int(line)
        except: return None, None
        
        f.readline() # Pula a linha de comentário do .xyz
        V, atoms = [], []
        for i, line in enumerate(f):
            if i >= n_atoms: break
            parts = line.split()
            if len(parts) < 4: continue
            atom = parts[0]
            # Ignora Hidrogênios conforme pedido original do fluxo
            if atom.upper() == 'H': continue
            try:
                V.append([float(parts[1]), float(parts[2]), float(parts[3])])
                atoms.append(atom)
            except: continue
        return np.array(V), atoms

def write_coordinates(fname, atoms, V):
    """Escreve as coordenadas em um novo arquivo no formato .xyz."""
    with open(fname, "w") as output:
        output.write(f"{len(V)}\n\n")
        for i in range(len(V)):
            output.write(f"{atoms[i]:2s} {V[i, 0]:10.5f} {V[i, 1]:10.5f} {V[i, 2]:10.5f}\n")

####################################################
# BLOCO 2 - VISUALIZAÇÃO E TRATAMENTO DE IMAGENS
####################################################

# Tabelas de referência para propriedades atômicas
ATOMS_DICT = {"H":1,"He":2,"Li":3,"Be":4,"B":5,"C":6,"N":7,"O":8,"F":9,"Ne":10,"Na":11,"Mg":12,"Al":13,"Si":14,"P":15,"S":16,"Cl":17,"Ar":18,"K":19,"Ca":20,"Sc":21,"Ti":22,"V":23,"Cr":24,"Mn":25,"Fe":26,"Co":27,"Ni":28,"Cu":29,"Zn":30,"Ga":31,"Ge":32,"As":33,"Se":34,"Br":35,"Kr":36,"Rb":37,"Sr":38,"Y":39,"Zr":40,"Nb":41,"Mo":42,"Tc":43,"Ru":44,"Rh":45,"Pd":46,"Ag":47,"Cd":48,"In":49,"Sn":50,"Sb":51,"Te":52,"I":53,"Xe":54,"Cs":55,"Ba":56,"La":57,"Ce":58,"Pr":59,"Nd":60,"Pm":61,"Sm":62,"Eu":63,"Gd":64,"Tb":65,"Dy":66,"Ho":67,"Er":68,"Tm":69,"Yb":70,"Lu":71,"Hf":72,"Ta":73,"W":74,"Re":75,"Os":76,"Ir":77,"Pt":78,"Au":79,"Hg":80,"Tl":81,"Pb":82,"Bi":83,"Po":84,"At":85,"Rn":86,"Fr":87,"Ra":88,"Ac":89,"Th":90,"Pa":91,"U":92,"Np":93,"Pu":94,"Am":95,"Cm":96,"Bk":97,"Cf":98,"Es":99,"Fm":100,"Md":101,"No":102,"Lr":103}
RADII_DICT = {"H":0.42,"He":1.6,"Li":0.68,"Be":0.352,"B":0.832,"C":0.72,"N":0.68,"O":0.68,"F":0.64,"Ne":1.12,"Na":0.972,"Mg":1.1,"Al":1.352,"Si":1.2,"P":1.036,"S":1.02,"Cl":1.0,"Ar":1.568,"K":1.328,"Ca":0.992,"Sc":1.44,"Ti":1.472,"V":1.328,"Cr":1.352,"Mn":1.352,"Fe":1.34,"Co":1.328,"Ni":1.62,"Cu":1.52,"Zn":1.448,"Ga":1.22,"Ge":1.168,"As":1.208,"Se":1.22,"Br":1.208,"Kr":1.6,"Rb":1.472,"Sr":1.12,"Y":1.78,"Zr":1.56,"Nb":1.48,"Mo":1.472,"Tc":1.352,"Ru":1.4,"Rh":1.448,"Pd":1.5,"Ag":1.592,"Cd":1.688,"In":1.632,"Sn":1.46,"Sb":1.46,"Te":1.472,"I":1.4,"Xe":1.7,"Cs":1.672,"Ba":1.64,"La":1.872,"Ce":1.932,"Pr":1.992,"Nd":1.992,"Pm":1.992,"Sm":1.992,"Eu":1.992,"Gd":1.892,"Tb":1.86,"Dy":1.852,"Ho":1.84,"Er":1.828,"Tm":1.82,"Yb":1.94,"Lu":1.82,"Hf":1.868,"Ta":1.432,"W":1.368,"Re":1.352,"Os":1.368,"Ir":1.32,"Pt":1.5,"Au":1.5,"Hg":1.7,"Tl":1.552,"Pb":1.54,"Bi":1.54,"Po":1.68,"At":1.208,"Rn":1.9,"Fr":1.8,"Ra":1.432,"Ac":1.18,"Th":1.02,"Pa":0.888,"U":0.968,"Np":0.952,"Pu":0.928,"Am":0.92,"Cm":0.912,"Bk":0.9,"Cf":0.888,"Es":0.88,"Fm":0.872,"Md":0.86,"No":0.848,"Lr":0.84,"Ln":1.992}

def get_atomic_number(atom): return ATOMS_DICT.get(atom, 6)
def get_radius(atom): return RADII_DICT.get(atom, 0.72)

def set_molecule_actor(atoms, coords, color):
    """Cria um objeto (actor) do VTK para representar a molécula em 3D."""
    mol = vtk.vtkMolecule()
    mol_list = []
    # Adiciona átomos à molécula
    for a, c in zip(atoms, coords):
        mol_list.append(mol.AppendAtom(get_atomic_number(a), c[0], c[1], c[2]))
    
    # Adiciona ligações baseadas na distância entre os átomos
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            d = np.linalg.norm(coords[i] - coords[j])
            if d < (get_radius(atoms[i]) + get_radius(atoms[j])) * 1.16:
                mol.AppendBond(mol_list[i], mol_list[j], 1)
    
    # Define cores e raios para cada átomo
    m_colors = vtk.vtkDoubleArray(); m_colors.SetName("Colors"); m_colors.SetNumberOfComponents(3); m_colors.SetNumberOfTuples(len(atoms))
    m_radii = vtk.vtkFloatArray(); m_radii.SetName("radii_array"); m_radii.SetNumberOfComponents(1); m_radii.SetNumberOfTuples(len(atoms))
    for i, atom in enumerate(atoms):
        m_colors.InsertTypedTuple(i, color)
        r = get_radius(atom)
        m_radii.InsertTypedTuple(i, [r * 0.1 if get_atomic_number(atom) > 56 else r * 0.2])
    
    mol.GetVertexData().AddArray(m_radii); mol.GetVertexData().AddArray(m_colors)
    mapper = vtk.vtkMoleculeMapper(); mapper.SetInputData(mol); mapper.UseBallAndStickSettings()
    mapper.SetAtomicRadiusScaleFactor(0.20); mapper.SetBondRadius(0.08); mapper.SetMapScalars(False)
    mapper.SetInputArrayToProcess(0, 0, 0, 4, "Colors"); mapper.SetAtomicRadiusArrayName("radii_array"); mapper.SetAtomicRadiusTypeToCustomArrayRadius()
    actor = vtk.vtkActor(); actor.SetMapper(mapper); return actor

def process_image(png_fname, rmsd_val):
    """Realiza o pós-processamento da imagem: recorta margens e adiciona o valor do RMSD."""
    if not os.path.exists(png_fname): return
    try:
        img = Image.open(png_fname)
        bg = Image.new(img.mode, img.size, (255, 255, 255))
        diff = ImageChops.difference(img, bg)
        diff = ImageChops.add(diff, diff, scale=2.0, offset=-100)
        bbox = diff.getbbox()
        if bbox:
            trim = (1.05, 0.90)
            bbox = (int(bbox[0]*trim[1]), int(bbox[1]*trim[1]), int(bbox[2]*trim[0]), int(bbox[3]*trim[0]))
            img = img.crop(bbox)
        
        draw = ImageDraw.Draw(img)
        # Tenta carregar uma fonte arial, se falhar usa a padrão
        try: font = ImageFont.truetype("arial.ttf", 80)
        except: font = ImageFont.load_default()
        draw.text((50, img.size[1] * 0.9), f"RMSD: {rmsd_val:.4f}", (0, 0, 0), font=font)
        img.save(png_fname)
    except Exception as e:
        print(f"Erro ao processar imagem {png_fname}: {e}")

####################################################
# BLOCO 3 - ORQUESTRADOR E INTERFACE
####################################################

def select_items(items, title=""):
    """Auxiliar para seleção granular de itens através do terminal."""
    print(f"\n>>> {title}")
    print("-" * 60)
    for i, item in enumerate(items):
        print(f"[{i+1}] {item}")
    print("-" * 60)
    choice = input("\nEnter numbers (e.g., 1 3 5), 'a' for all, or 'n' for none: ").strip().lower()
    if choice == 'a': return set(items)
    if choice == 'n' or not choice: return set()
    try:
        indices = [int(x) - 1 for x in choice.split() if x.isdigit()]
        return set([items[i] for i in indices if 0 <= i < len(items)])
    except:
        print("Invalid selection. Proceeding with none.")
        return set()

def main():
    print("=" * 60)
    print("      RMSD ANALYSIS & MOLECULAR ALIGNMENT TOOL")
    print("=" * 60)
    
    # PASSO 1: Configuração de entrada
    print("\n[STEP 1] INPUT CONFIGURATION")
    print("1 - Use subfolders (INPUT/referencias and INPUT/estruturas)")
    print("2 - Use reference files (.ref) and structures (.xyz) in root/INPUT")
    input_choice = input("\nSelect option (1 or 2): ").strip()
    
    dir_input = "INPUT"
    dir_ref = os.path.join(dir_input, "referencias")
    dir_str = os.path.join(dir_input, "estruturas")
    
    if input_choice == '1':
        for d in [dir_ref, dir_str]:
            if not os.path.exists(d): os.makedirs(d)
        print(f"\n[PAUSE] Place REFERENCES in: {dir_ref}")
        print(f"        Place STRUCTURES in : {dir_str}")
    else:
        if not os.path.exists(dir_input): os.makedirs(dir_input)
        print("\n[PAUSE] Place .ref and .xyz files in the root or 'INPUT' folder.")

    input("\nPress ENTER once files are ready...")

    # Coleta de arquivos baseada na escolha do usuário
    if input_choice == '1':
        ref_files = glob.glob(os.path.join(dir_ref, "*.xyz"))
        str_files = glob.glob(os.path.join(dir_str, "*.xyz"))
    else:
        ref_files = glob.glob("*.ref") + glob.glob(os.path.join(dir_input, "*.ref"))
        str_files = glob.glob("*.xyz") + glob.glob(os.path.join(dir_input, "*.xyz"))

    ref_files = list(dict.fromkeys([os.path.abspath(f) for f in ref_files]))
    str_files = list(dict.fromkeys([os.path.abspath(f) for f in str_files]))
    
    # Identifica os pares de comparação usando Regex para evitar sobreposições
    comparisons_list = [] # Lista de tuplas (caminho_ref, caminho_str, rótulo)
    for r_path in ref_files:
        r_name = os.path.basename(r_path).replace(".ref", "").replace(".xyz", "")
        # Regex: deve começar com r_name seguido de separador (_ ou .) ou fim da string
        pattern = re.compile(rf"^{re.escape(r_name)}([._]|$)")
        for s_path in str_files:
            if s_path == r_path: continue
            if pattern.search(os.path.basename(s_path)):
                label = f"{os.path.basename(s_path)} vs {r_name}"
                comparisons_list.append((r_path, s_path, label))

    # PASSO 2: Estatísticas de identificação
    print("\n[STEP 2] IDENTIFIED COMPARISONS")
    if not comparisons_list:
        print("No matching pairs found. Check filenames."); return
    
    print(f"Total pairs found: {len(comparisons_list)}")
    for _, _, label in comparisons_list: print(f"  -> {label}")

    # PASSO 3: Opções de Saída Opcionais
    print("\n[STEP 3] OUTPUT OPTIONS")
    print("RMSD.txt is mandatory and will be generated for ALL pairs.")
    
    save_aligned = input("\nExport aligned structures (.xyz)? (y/n): ").strip().lower() == 'y'
    selected_xyz = set()
    if save_aligned:
        labels = [c[2] for c in comparisons_list]
        selected_xyz = select_items(labels, "SELECT STRUCTURES TO SAVE AS .XYZ")

    save_images = input("\nGenerate overlay images (.png)? (y/n): ").strip().lower() == 'y'
    selected_png = set()
    if save_images:
        if save_aligned and selected_xyz:
            if input("Use the same selection for images? (y/n): ").strip().lower() == 'y':
                selected_png = selected_xyz
            else:
                labels = [c[2] for c in comparisons_list]
                selected_png = select_items(labels, "SELECT STRUCTURES FOR IMAGES")
        else:
            labels = [c[2] for c in comparisons_list]
            selected_png = select_items(labels, "SELECT STRUCTURES FOR IMAGES")

    # Cria diretórios de saída se necessário
    out_dir = "OUTPUT"
    out_aligned = os.path.join(out_dir, "Alinhamentos")
    out_images = os.path.join(out_dir, "Imagens")
    for d in [out_dir, out_aligned if save_aligned else None, out_images if save_images else None]:
        if d and not os.path.exists(d): os.makedirs(d)

    # Inicializa o motor VTK se imagens forem solicitadas
    renderer = None; renderWin = None
    if save_images and selected_png:
        renderer = vtk.vtkRenderer(); renderer.SetBackground(1, 1, 1)
        renderWin = vtk.vtkRenderWindow(); renderWin.AddRenderer(renderer)
        renderWin.SetOffScreenRendering(1); renderWin.SetSize(1200, 1200)

    # PASSO 4: Loop principal de processamento
    print("\n[STEP 4] PROCESSING...")
    results = []
    for r_path, s_path, label in comparisons_list:
        P, atomsP = get_coordinates(s_path)
        Q, atomsQ = get_coordinates(r_path)
        
        # Validação básica de contagem de átomos
        if P is None or Q is None or len(P) != len(Q):
            print(f"  ! Skipping {label}: Atom count mismatch or file error.")
            continue

        print(f"  > Calculating RMSD: {label}")
        norm_r = rmsd(P, Q)
        kabs_r = kabsch_rmsd(P, Q)
        fit_res = fit(P.copy(), Q.copy())
        fit_r, P_fit = fit_res if fit_res else (999.999, P.copy())

        # Verifica se o usuário selecionou este par para as saídas opcionais
        do_xyz = label in selected_xyz
        do_png = label in selected_png

        if do_xyz or do_png:
            # Centralização e Alinhamento Final para exportação/imagem
            P_c = P_fit - centroid(P_fit); Q_c = Q - centroid(Q); P_final = rotate(P_c, Q_c)
            
            if do_xyz:
                print(f"    - Saving aligned XYZ...")
                # Reposiciona na origem da referência para contexto
                P_output = P_final + centroid(Q)
                out_path = os.path.join(out_aligned, f"{os.path.splitext(os.path.basename(s_path))[0]}_aligned.xyz")
                write_coordinates(out_path, atomsP, P_output)

            if do_png:
                print(f"    - Saving image...")
                act_ref = set_molecule_actor(atomsQ, Q_c, (1.0, 0.0, 0.0)) # Referência = Vermelho
                act_str = set_molecule_actor(atomsP, P_final, (0.0, 1.0, 0.0)) # Estrutura = Verde
                renderer.AddActor(act_ref); renderer.AddActor(act_str)
                renderer.ResetCamera(); renderWin.Render()
                png_path = os.path.join(out_images, f"{os.path.splitext(os.path.basename(s_path))[0]}.png")
                rl = vtk.vtkRenderLargeImage(); rl.SetInput(renderer); rl.SetMagnification(2)
                w = vtk.vtkPNGWriter(); w.SetInputConnection(rl.GetOutputPort()); w.SetFileName(png_path); w.Write()
                process_image(png_path, fit_r)
                renderer.RemoveActor(act_ref); renderer.RemoveActor(act_str)

        results.append((os.path.basename(s_path), norm_r, kabs_r, fit_r))

    # Gravação do Relatório Final de RMSD
    report_path = os.path.join(out_dir, "RMSD.txt")
    with open(report_path, "w") as f:
        f.write(f"{'Molecule':<35} {'Normal':>12} {'Kabsch':>12} {'Fitted':>12}\n" + "-" * 75 + "\n")
        for r in results: f.write(f"{r[0]:<35} {r[1]:>12.4f} {r[2]:>12.4f} {r[3]:>12.4f}\n")
    
    print(f"\n{'='*60}\nSUCCESS! Report: {report_path}\n{'='*60}")

if __name__ == "__main__":
    main()

