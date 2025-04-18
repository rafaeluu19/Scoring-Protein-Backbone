import streamlit as st
from ramachandran_scoring import calculate_ramachandran_score

st.title("Ramachandran Score Calculator")
st.write("Upload a PDB file to calculate the backbone score based on phi/psi angles.")

uploaded_file = st.file_uploader("Choose a PDB file", type="pdb")

if uploaded_file is not None:
    try:
        uploaded_file.seek(0)
        score = calculate_ramachandran_score(uploaded_file)
        st.success(f"✅ Ramachandran Score: {score:.3f}")
    except Exception as e:
        st.error(f"❌ Error processing file: {str(e)}")