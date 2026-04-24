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
# BLOCK 1 - MATHEMATICAL AND GEOMETRIC CALCULATIONS
####################################################

def centroid(X):
    """Calculates the geometric center (centroid) of a set of points."""
    if len(X) == 0: return np.array([0,0,0])
    return np.mean(X, axis=0)

def rmsd(P, Q):
    """Calculates the simple RMSD between two matrices of the same dimension."""
    if P.shape != Q.shape: return 999.999
    diff = P - Q
    return np.sqrt((diff * diff).sum() / P.shape[0])

def kabsch(P, Q):
    """Implements the Kabsch Algorithm to find the ideal rotation matrix."""
    C = np.dot(np.transpose(P), Q)
    V, S, W = np.linalg.svd(C)
    # Ensures the rotation matrix does not result in a reflection
    if (np.linalg.det(V) * np.linalg.det(W)) < 0.0:
        V[:, -1] = -V[:, -1]
    return np.dot(V, W)

def rotate(P, Q):
    """Rotates matrix P to best align with matrix Q."""
    U = kabsch(P, Q)
    return np.dot(P, U)

def kabsch_rmsd(P, Q):
    """Calculates the RMSD after ideal Kabsch rotation (requires centering)."""
    P_c = P - centroid(P)
    Q_c = Q - centroid(Q)
    P_rot = rotate(P_c, Q_c)
    return rmsd(P_rot, Q_c)

def fit(P, Q, max_iter=1000):
    """Translation optimization to find the best fit before final rotation."""
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
            
            # Test move in positive direction
            rmsd_new = kabsch_rmsd(P_current + temp, Q)
            if rmsd_new < rmsd_best:
                rmsd_best = rmsd_new
                P_current[:, i] += step_size[i]
                improved = True
            else:
                # Test move in negative direction
                rmsd_new = kabsch_rmsd(P_current - temp, Q)
                if rmsd_new < rmsd_best:
                    rmsd_best = rmsd_new
                    P_current[:, i] -= step_size[i]
                    improved = True
                else:
                    # Reduce step size if no move helped
                    step_size[i] /= 2
                    
        if (step_size <= threshold).all(): break
            
    if n_iter > max_iter: return None
    return rmsd_best, P_current

def get_coordinates(filename):
    """Reads a .xyz file and extracts coordinates (ignores Hydrogen atoms)."""
    if not os.path.exists(filename): return None, None
    with open(filename, "r") as f:
        try:
            line = f.readline().strip()
            if not line: return None, None
            n_atoms = int(line)
        except: return None, None
        
        f.readline() # Skip comment line
        V, atoms = [], []
        for i, line in enumerate(f):
            if i >= n_atoms: break
            parts = line.split()
            if len(parts) < 4: continue
            atom = parts[0]
            # Ignore Hydrogen as per original workflow
            if atom.upper() == 'H': continue
            try:
                V.append([float(parts[1]), float(parts[2]), float(parts[3])])
                atoms.append(atom)
            except: continue
        return np.array(V), atoms

def write_coordinates(fname, atoms, V):
    """Writes coordinates to a new file in .xyz format."""
    with open(fname, "w") as output:
        output.write(f"{len(V)}\n\n")
        for i in range(len(V)):
            output.write(f"{atoms[i]:2s} {V[i, 0]:10.5f} {V[i, 1]:10.5f} {V[i, 2]:10.5f}\n")

####################################################
# BLOCK 2 - VISUALIZATION AND IMAGE TREATMENT
####################################################

# Reference tables for atomic properties
ATOMS_DICT = {"H":1,"He":2,"Li":3,"Be":4,"B":5,"C":6,"N":7,"O":8,"F":9,"Ne":10,"Na":11,"Mg":12,"Al":13,"Si":14,"P":15,"S":16,"Cl":17,"Ar":18,"K":19,"Ca":20,"Sc":21,"Ti":22,"V":23,"Cr":24,"Mn":25,"Fe":26,"Co":27,"Ni":28,"Cu":29,"Zn":30,"Ga":31,"Ge":32,"As":33,"Se":34,"Br":35,"Kr":36,"Rb":37,"Sr":38,"Y":39,"Zr":40,"Nb":41,"Mo":42,"Tc":43,"Ru":44,"Rh":45,"Pd":46,"Ag":47,"Cd":48,"In":49,"Sn":50,"Sb":51,"Te":52,"I":53,"Xe":54,"Cs":55,"Ba":56,"La":57,"Ce":58,"Pr":59,"Nd":60,"Pm":61,"Sm":62,"Eu":63,"Gd":64,"Tb":65,"Dy":66,"Ho":67,"Er":68,"Tm":69,"Yb":70,"Lu":71,"Hf":72,"Ta":73,"W":74,"Re":75,"Os":76,"Ir":77,"Pt":78,"Au":79,"Hg":80,"Tl":81,"Pb":82,"Bi":83,"Po":84,"At":85,"Rn":86,"Fr":87,"Ra":88,"Ac":89,"Th":90,"Pa":91,"U":92,"Np":93,"Pu":94,"Am":95,"Cm":96,"Bk":97,"Cf":98,"Es":99,"Fm":100,"Md":101,"No":102,"Lr":103}
RADII_DICT = {"H":0.42,"He":1.6,"Li":0.68,"Be":0.352,"B":0.832,"C":0.72,"N":0.68,"O":0.68,"F":0.64,"Ne":1.12,"Na":0.972,"Mg":1.1,"Al":1.352,"Si":1.2,"P":1.036,"S":1.02,"Cl":1.0,"Ar":1.568,"K":1.328,"Ca":0.992,"Sc":1.44,"Ti":1.472,"V":1.328,"Cr":1.352,"Mn":1.352,"Fe":1.34,"Co":1.328,"Ni":1.62,"Cu":1.52,"Zn":1.448,"Ga":1.22,"Ge":1.168,"As":1.208,"Se":1.22,"Br":1.208,"Kr":1.6,"Rb":1.472,"Sr":1.12,"Y":1.78,"Zr":1.56,"Nb":1.48,"Mo":1.472,"Tc":1.352,"Ru":1.4,"Rh":1.448,"Pd":1.5,"Ag":1.592,"Cd":1.688,"In":1.632,"Sn":1.46,"Sb":1.46,"Te":1.472,"I":1.4,"Xe":1.7,"Cs":1.672,"Ba":1.64,"La":1.872,"Ce":1.932,"Pr":1.992,"Nd":1.992,"Pm":1.992,"Sm":1.992,"Eu":1.992,"Gd":1.892,"Tb":1.86,"Dy":1.852,"Ho":1.84,"Er":1.828,"Tm":1.82,"Yb":1.94,"Lu":1.82,"Hf":1.868,"Ta":1.432,"W":1.368,"Re":1.352,"Os":1.368,"Ir":1.32,"Pt":1.5,"Au":1.5,"Hg":1.7,"Tl":1.552,"Pb":1.54,"Bi":1.54,"Po":1.68,"At":1.208,"Rn":1.9,"Fr":1.8,"Ra":1.432,"Ac":1.18,"Th":1.02,"Pa":0.888,"U":0.968,"Np":0.952,"Pu":0.928,"Am":0.92,"Cm":0.912,"Bk":0.9,"Cf":0.888,"Es":0.88,"Fm":0.872,"Md":0.86,"No":0.848,"Lr":0.84,"Ln":1.992}

def get_atomic_number(atom): return ATOMS_DICT.get(atom, 6)
def get_radius(atom): return RADII_DICT.get(atom, 0.72)

def set_molecule_actor(atoms, coords, color):
    """Creates a VTK actor object to represent the molecule in 3D."""
    mol = vtk.vtkMolecule()
    mol_list = []
    # Add atoms to molecule
    for a, c in zip(atoms, coords):
        mol_list.append(mol.AppendAtom(get_atomic_number(a), c[0], c[1], c[2]))
    
    # Add bonds based on distance between atoms
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            d = np.linalg.norm(coords[i] - coords[j])
            if d < (get_radius(atoms[i]) + get_radius(atoms[j])) * 1.16:
                mol.AppendBond(mol_list[i], mol_list[j], 1)
    
    # Define colors and radii for each atom
    m_colors = vtk.vtkDoubleArray(); m_colors.SetName("Colors"); m_colors.SetNumberOfComponents(3); m_colors.SetNumberOfTuples(len(atoms))
    m_radii  = vtk.vtkFloatArray(); m_radii.SetName("radii_array"); m_radii.SetNumberOfComponents(1); m_radii.SetNumberOfTuples(len(atoms))
    
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
    """Performs image post-processing: trims margins and adds RMSD value."""
    if not os.path.exists(png_fname): return
    try:
        img  = Image.open(png_fname)
        bg   = Image.new(img.mode, img.size, (255, 255, 255))
        diff = ImageChops.difference(img, bg)
        diff = ImageChops.add(diff, diff, scale=2.0, offset=-100)
        bbox = diff.getbbox()
        
        if bbox:
            trim = (1.05, 0.90)
            bbox = (int(bbox[0]*trim[1]), int(bbox[1]*trim[1]), int(bbox[2]*trim[0]), int(bbox[3]*trim[0]))
            img = img.crop(bbox)
        
        draw = ImageDraw.Draw(img)
        
        # Tries to load arial font, falls back to default
        try   : font = ImageFont.truetype("arial.ttf", 80)
        except: font = ImageFont.load_default()
        
        draw.text((50, img.size[1] * 0.9), f"RMSD: {rmsd_val:.4f}", (0, 0, 0), font=font)
        img.save(png_fname)
    except Exception as e:
        print(f"Error processing image {png_fname}: {e}")

####################################################
# BLOCK 3 - ORCHESTRATOR AND INTERFACE
####################################################

def select_items(items, title=""):
    """Helper for granular selection of items through the terminal."""
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

# Helper function for keyword mode (parses selection as 'a', 'n', or list of indices)
def parse_selection(sel_str, labels):
    if sel_str == 'a':
        return set(labels)
    if sel_str == 'n' or not sel_str:
        return set()
    # Try to interpret as comma or space separated numbers
    parts = re.split(r'[,\s]+', sel_str.strip())
    indices = []
    for p in parts:
        if p.isdigit():
            indices.append(int(p) - 1)
    return set([labels[i] for i in indices if 0 <= i < len(labels)])

def main(mode='keyword'):  # default mode = keyword
    """
    mode parameter:
        'keyword' - all answers in a single input()
        'prompt'  - original interactive questions
    """
    # ---------- 1. Path definitions (no creation yet) ----------
    base_dir = os.path.dirname(os.path.abspath(__file__))
    print("=" * 60)
    print("      RMSD ANALYSIS & MOLECULAR ALIGNMENT TOOL")
    print("=" * 60)
    print(f"[PATH] Running from: {base_dir}")

    dir_input = os.path.join(base_dir, "INPUT")
    dir_ref   = os.path.join(dir_input, "referencias")
    dir_str   = os.path.join(dir_input, "estruturas")
    dir_out   = os.path.join(base_dir, "OUTPUT")

    # ---------- 2. Reading user choices ----------
    if mode == 'prompt':
        print("\n[STEP 1] INPUT CONFIGURATION")
        print("1 - Use subfolders (Prefix matching: reference.xyz vs structure.xyz)")
        print("2 - Use root folder (Cross-match: ALL *.ref files vs ALL *.xyz files)")
        input_choice = input("\nSelect option (1 or 2) [Default = 1]: ").strip() or '1'
        
        if input_choice == '1':
            print(f"\n[PAUSE] Place REFERENCES (.xyz) in: {dir_ref}")
            print(f"        Place STRUCTURES (.xyz) in : {dir_str}")
        else:
            print(f"\n[PAUSE] Place .ref and .xyz files in: {base_dir} or {dir_input}")
            print("        Every *.ref file will be compared against every *.xyz file.")
            
        input("\nPress ENTER once files are ready to process...")
        kw_xyz_sel = None
        kw_img_sel = None
    else:  # keyword mode
        print("\n" + "=" * 80)
        print("                     SINGLE INPUT CONFIGURATION (KEYWORD MODE)")
        print("=" * 80)
        print("Instructions: Type all parameters in a single line, separated by spaces.")
        print("\nPARAMETER ORDER AND HOW TO ANSWER:")
        print("-" * 80)
        print("1. [FOLDER]   : '1' for subfolders (Prefix matching: reference.xyz vs structure.xyz)")
        print("                '2' for root folder (Cross-match: ALL *.ref files vs ALL *.xyz files)")
        print("2. [SAVE XYZ] : 'y' (yes) or 'n' (no) to export aligned structures")
        print("3. [SEL XYZ]  : (Skip if answered 'n' above) Use 'a' (all), 'n' (none)")
        print("                or indices separated by commas (e.g., 1,3,5)")
        print("4. [SAVE PNG] : 'y' (yes) or 'n' (no) to generate overlay images")
        print("5. [SEL PNG]  : (Skip if answered 'n' above) Same logic as XYZ selection ('a', 'n' or numbers)")
        print("-" * 80)
        print("RESPONSE EXAMPLES:")
        print("  > 1 y a y a      (Subfolders, save all XYZ, generate all images)")
        print("  > 1 n y 1,2,3    (Subfolders, no XYZ, generate images for pairs 1, 2 and 3)")
        print("  > 2 y 1 n        (Root folder, CROSS-MATCH pair 1, no images)")
        print("=" * 80)
        
        user_input = input("ENTER YOUR RESPONSE: ").strip().split()
        if not user_input:
            print("\n[ERROR] No response provided. Exiting.")
            return
        try:
            # Input parsing
            input_choice = user_input[0]
            if input_choice not in ('1','2'): raise ValueError
            
            idx = 1
            save_aligned_choice = user_input[idx].lower()
            idx += 1
            
            kw_xyz_sel = None
            if save_aligned_choice == 'y':
                if idx >= len(user_input): raise ValueError
                kw_xyz_sel = user_input[idx]
                idx += 1
            
            save_images_choice = user_input[idx].lower()
            idx += 1
            
            kw_img_sel = None
            if save_images_choice == 'y':
                if idx >= len(user_input): raise ValueError
                kw_img_sel = user_input[idx]
        except:
            print("\n[ERROR] Invalid or incomplete response format.")
            print("Ensure you follow the order and provide selections for each chosen 'y'.")
            return

    # ---------- 3. Conditional Directory Creation ----------
    folders_to_create = []
    if input_choice == '1':
        # Scaffolding for subfolders mode
        folders_to_create = [dir_input, dir_ref, dir_str, dir_out]
    else:
        # Only ensure output folder exists for root mode
        folders_to_create = [dir_out]

    created_folders = []
    for d in folders_to_create:
        if not os.path.exists(d):
            os.makedirs(d)
            created_folders.append(os.path.basename(d))
    
    if created_folders:
        print(f"\n[INFO] Folders created at script location: {', '.join(created_folders)}")

    # --- File collection ---
    extensions = ["*.xyz"]
    ref_extensions = ["*.ref"]
    ref_files = []
    str_files = []

    if input_choice == '1':
        # Subfolders mode: Search ONLY in the created subfolders
        for ext in extensions:
            ref_files.extend(glob.glob(os.path.join(dir_ref, ext)))
            str_files.extend(glob.glob(os.path.join(dir_str, ext)))
    else:
        # Root mode: Search in base folder AND optionally in INPUT if it happens to exist
        for rext in ref_extensions:
            ref_files.extend(glob.glob(os.path.join(base_dir, rext)))
            if os.path.exists(dir_input):
                ref_files.extend(glob.glob(os.path.join(dir_input, rext)))
        for ext in extensions:
            str_files.extend(glob.glob(os.path.join(base_dir, ext)))
            if os.path.exists(dir_input):
                str_files.extend(glob.glob(os.path.join(dir_input, ext)))

    ref_files = list(dict.fromkeys([os.path.abspath(f) for f in ref_files]))
    str_files = list(dict.fromkeys([os.path.abspath(f) for f in str_files]))

    comparisons_list = []
    for r_path in ref_files:
        r_base = os.path.basename(r_path)
        r_prefix = re.split(r'[_.]', r_base)[0]
        for s_path in str_files:
            if s_path == r_path: continue
            s_base = os.path.basename(s_path)
            s_prefix = re.split(r'[_.]', s_base)[0]
            
            # Match logic: 
            # Option 1: Prefix matching (e.g., ORIGII_ref.xyz vs ORIGII_str.xyz)
            # Option 2: Cross-match (every .ref file compared against every .xyz file)
            if input_choice == '2' or r_prefix.lower() == s_prefix.lower():
                label = f"{s_base} vs {r_base}"
                comparisons_list.append((r_path, s_path, label))

    print("\n[STEP 2] IDENTIFIED COMPARISONS")
    if not comparisons_list:
        print("No matching pairs found. Check filenames."); return
    print(f"Total pairs found: {len(comparisons_list)}")
    for _, _, label in comparisons_list: print(f"  -> {label}")

    # --- Get selections ---
    labels = [c[2] for c in comparisons_list]
    selected_xyz = set()
    selected_png = set()

    if mode == 'prompt':
        print("\n[STEP 3] OUTPUT OPTIONS")
        print("RMSD.txt is mandatory and will be generated for ALL pairs.")
        save_aligned = input("\nExport aligned structures (.xyz)? (y/n): ").strip().lower() == 'y'
        if save_aligned:
            selected_xyz = select_items(labels, "SELECT STRUCTURES TO SAVE AS .XYZ")
        save_images = input("\nGenerate overlay images (.png)? (y/n): ").strip().lower() == 'y'
        if save_images:
            if save_aligned and selected_xyz:
                if input("Use the same selection for images? (y/n): ").strip().lower() == 'y':
                    selected_png = selected_xyz
                else:
                    selected_png = select_items(labels, "SELECT STRUCTURES FOR IMAGES")
            else:
                selected_png = select_items(labels, "SELECT STRUCTURES FOR IMAGES")
    else:  # keyword mode
        save_aligned = (kw_xyz_sel is not None)
        save_images = (kw_img_sel is not None)
        if save_aligned: selected_xyz = parse_selection(kw_xyz_sel, labels)
        if save_images: selected_png = parse_selection(kw_img_sel, labels)

    out_aligned = os.path.join(dir_out, "Alinhamentos")
    out_images  = os.path.join(dir_out, "Imagens")
    for d in [dir_out, out_aligned if save_aligned else None, out_images if save_images else None]:
        if d and not os.path.exists(d): os.makedirs(d, exist_ok=True)

    renderer = None; renderWin = None
    if save_images and selected_png:
        renderer  = vtk.vtkRenderer(); renderer.SetBackground(1, 1, 1)
        renderWin = vtk.vtkRenderWindow(); renderWin.AddRenderer(renderer)
        renderWin.SetOffScreenRendering(1); renderWin.SetSize(1200, 1200)

    print("\n[STEP 4] PROCESSING...")
    results = []
    for r_path, s_path, label in comparisons_list:
        P, atomsP = get_coordinates(s_path)
        Q, atomsQ = get_coordinates(r_path)
        if P is None or Q is None or len(P) != len(Q):
            print(f"  ! Skipping {label}: Atom count mismatch or file error.")
            continue

        print(f"  > Calculating RMSD: {label}")
        norm_r  = rmsd(P, Q)
        kabs_r  = kabsch_rmsd(P, Q)
        fit_res = fit(P.copy(), Q.copy())
        fit_r, P_fit = fit_res if fit_res else (999.999, P.copy())

        do_xyz = label in selected_xyz
        do_png = label in selected_png

        if do_xyz or do_png:
            P_c = P_fit - centroid(P_fit); Q_c = Q - centroid(Q); P_final = rotate(P_c, Q_c)
            if do_xyz:
                P_output = P_final + centroid(Q)
                fname_xyz = f"{os.path.splitext(os.path.basename(s_path))[0]}_aligned.xyz"
                out_path = os.path.join(out_aligned, fname_xyz)
                print(f"    - Saving aligned XYZ: {fname_xyz}")
                write_coordinates(out_path, atomsP, P_output)
            if do_png:
                fname_png = f"{os.path.splitext(os.path.basename(s_path))[0]}.png"
                png_path = os.path.join(out_images, fname_png)
                print(f"    - Saving image: {fname_png}")
                act_ref = set_molecule_actor(atomsQ, Q_c, (1.0, 0.0, 0.0))
                act_str = set_molecule_actor(atomsP, P_final, (0.0, 1.0, 0.0))
                renderer.AddActor(act_ref); renderer.AddActor(act_str)
                renderer.ResetCamera(); renderWin.Render()
                rl = vtk.vtkRenderLargeImage(); rl.SetInput(renderer); rl.SetMagnification(2)
                w = vtk.vtkPNGWriter(); w.SetInputConnection(rl.GetOutputPort()); w.SetFileName(png_path); w.Write()
                process_image(png_path, fit_r)
                renderer.RemoveActor(act_ref); renderer.RemoveActor(act_str)

        results.append((os.path.basename(s_path), norm_r, kabs_r, fit_r))

    report_path = os.path.join(dir_out, "RMSD.txt")
    try:
        with open(report_path, "w") as f:
            f.write(f"{'Molecule':<35} {'Normal':>12} {'Kabsch':>12} {'Fitted':>12}\n" + "-" * 75 + "\n")
            for r in results:
                f.write(f"{r[0]:<35} {r[1]:>12.4f} {r[2]:>12.4f} {r[3]:>12.4f}\n")
        print(f"\n{'='*60}\nSUCCESS! Report generated at: {report_path}\n{'='*60}")
    except Exception as e:
        print(f"\n[ERROR] Failed to write report: {e}")

if __name__ == "__main__":
    main()
