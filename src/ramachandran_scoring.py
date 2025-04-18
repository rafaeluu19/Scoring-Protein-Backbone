from Bio.PDB import PDBParser, PPBuilder
import io

#Loads the PDB file
#pdb_filename = "data/1ubq.pdb" #Replace with your actual PDB filename
#parser = PDBParser(QUIET=True)
#structure = parser.get_structure("protein", pdb_filename)

#Use the first model
#model = structure[0]

#Use PPBuilder to get polypeptide chains 
#ppb = PPBuilder()

#List to store all (phi, psi) pairs in degrees
#angles = []
#all_residues = []

#for pp_index, pp in enumerate(ppb.build_peptides(model)):
    #print(f"\n--- Polypeptide{pp_index + 1} ---")

    #phi_psi = pp.get_phi_psi_list()

    #for i, residue in enumerate(pp):
        #resname = residue.get_resname()
        #resnum = residue.get_id()[1]

        #phi, psi = phi_psi[i]

        #if phi is not None and psi is not None:
           # phi_deg = math.degrees(phi)
            #psi_deg = math.degrees(psi)
            #angles.append((phi_deg, psi_deg)) 
            #all_residues.append(residue)

def calculate_ramachandran_score(pdb_file):
    pdb_file.seek(0)
    file_content = pdb_file.read().decode("utf-8")
    handle = io.StringIO(file_content)
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", handle)
    model = structure[0]
    ppb = PPBuilder()

    score = 0
    count = 0
    
    for pp in ppb.build_peptides(model):
        phi_psi = pp.get_phi_psi_list()
        for (phi,psi),residue in zip(phi_psi,pp):
            if phi is None or psi is None:
                continue

            count += 1
            resname = residue.get_resname()

            if resname == "GLY":
                if -180 <= phi <= 180 and -180 <= psi <= 180:
                    score += 1.0 

            elif resname == "PRO":
                if -90 <= phi <=-35 and 145 <= psi <= 180:
                    score += 1.0
                else: 
                    score += 0.2

            else: 
                if -150 <= phi <= -30 and -90<= psi <= -30:
                    score += 1.0
                elif -180 <= phi <= -50 and 90 <= psi <= 180:
                    score += 0.8
                elif -180 <= phi <= 180 and -180 <= psi <= 180:
                    score += 0.2
                else:
                    score += 0
    return score/count if count else 0.0

#for i, ((phi, psi),residue)in enumerate(zip(angles, all_residues)):
    #resname = residue.get_resname()
    #score = score_phi_psi(phi,psi)
    #total_score += score
    #print(f"Residue {i+1}: ϕ = {phi:.2f}°, ψ = {psi:.2f}°, Score = {score:.2f}")

#Average Score
#average_score = total_score / len(angles) if angles else 0
#print(f"\nFinal Ramachandran Score: {average_score:.3f} (average over {len(angles)} residues)")