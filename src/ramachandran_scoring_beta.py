import io
import math
from Bio.PDB import PDBParser, PPBuilder
from ramachandran_scoring import calculate_ramachandran_score
# --- Scoring Logic Function ---
def score_phi_psi(phi, psi, resname):
    if phi is None or psi is None:
        return 0.0

    if resname == "GLY":
        if -180 <= phi <= 180 and -180 <= psi <= 180:
            return 1.0

    elif resname == "PRO":
        if -90 <= phi <= -35 and 145 <= psi <= 180:
            return 1.0
        else:
            return 0.2

    else:
        if -150 <= phi <= -30 and -90 <= psi <= -30:
            return 1.0
        elif -180 <= phi <= -50 and 90 <= psi <= 180:
            return 0.8
        elif -180 <= phi <= 180 and -180 <= psi <= 180:
            return 0.2
    return 0.0

# --- Main scoring function for file-like object ---
def calculate_ramachandran_score(pdb_file):
    file_content = pdb_file.read().decode("utf-8")
    handle = io.StringIO(file_content)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", handle)
    model = structure[0]
    ppb = PPBuilder()

    total_score = 0.0
    count = 0

    for pp in ppb.build_peptides(model):
        phi_psi = pp.get_phi_psi_list()
        for (phi, psi), residue in zip(phi_psi, pp):
            if phi is None or psi is None:
                continue

            resname = residue.get_resname()
            score = score_phi_psi(math.degrees(phi), math.degrees(psi), resname)
            total_score += score
            count += 1

    return total_score / count if count else 0.0

# --- Optional: Run directly on a local file (for testing outside Streamlit) ---
def test_local_file():
    pdb_filename = "../data/1ubq.pdb"
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_filename)
    model = structure[0]
    ppb = PPBuilder()

    angles = []
    all_residues = []

    for pp in ppb.build_peptides(model):
        phi_psi = pp.get_phi_psi_list()
        for i, residue in enumerate(pp):
            phi, psi = phi_psi[i]
            if phi is not None and psi is not None:
                phi_deg = math.degrees(phi)
                psi_deg = math.degrees(psi)
                angles.append((phi_deg, psi_deg))
                all_residues.append(residue)

    total_score = 0.0
    for i, ((phi, psi), residue) in enumerate(zip(angles, all_residues)):
        resname = residue.get_resname()
        score = score_phi_psi(phi, psi, resname)
        total_score += score
        print(f"Residue {i+1}: ϕ = {phi:.2f}°, ψ = {psi:.2f}°, Score = {score:.2f}")

    average_score = total_score / len(angles) if angles else 0
    print(f"\nFinal Ramachandran Score: {average_score:.3f} (average over {len(angles)} residues)")

def test_uploaded_file():
    with open("../data/1ubq.pdb","rb") as f:
        score = calculate_ramachandran_score(f)
        print(f"Final Ramachandran Score (from uploaded file): {score:.3f}")

# Uncomment below to test manually
test_local_file()
test_uploaded_file()
