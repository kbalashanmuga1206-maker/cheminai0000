import io
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, rdMolDescriptors
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image as PDFImage, Table, TableStyle
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import inch

st.set_page_config(page_title="Molecule Studio", page_icon="🧪", layout="wide")
st.title("🧪 Molecule Studio")
st.caption("A simple chemistry project for SMILES validation, 2D/3D visualization and PDF reporting.")

examples = {
    "Caffeine": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "Ethanol": "CCO",
    "Benzene": "c1ccccc1",
}
choice = st.selectbox("Choose an example or enter your own SMILES", ["Custom"] + list(examples.keys()))
default_smiles = examples.get(choice, "")
smiles = st.text_input("SMILES notation", value=default_smiles,
                       placeholder="Example: CCO")

if not smiles.strip():
    st.info("Enter a SMILES string to begin.")
    st.stop()

mol = Chem.MolFromSmiles(smiles.strip())
if mol is None:
    st.error("Invalid SMILES. Please check the notation and try again.")
    st.stop()

st.success("Valid SMILES notation")
canonical = Chem.MolToSmiles(mol)
formula = rdMolDescriptors.CalcMolFormula(mol)
mw = Descriptors.MolWt(mol)
logp = Descriptors.MolLogP(mol)
hbd = Descriptors.NumHDonors(mol)
hba = Descriptors.NumHAcceptors(mol)
atoms = mol.GetNumAtoms()
bonds = mol.GetNumBonds()

st.subheader("Molecular information")
c1, c2, c3, c4 = st.columns(4)
c1.metric("Molecular formula", formula)
c2.metric("Molecular weight", f"{mw:.2f} g/mol")
c3.metric("Atoms", atoms)
c4.metric("Bonds", bonds)
st.write("**Canonical SMILES:**", canonical)
st.write(f"**LogP:** {logp:.2f}  |  **H-bond donors:** {hbd}  |  **H-bond acceptors:** {hba}")

st.subheader("2D molecular structure")
img = Draw.MolToImage(mol, size=(650, 450))
img_buffer = io.BytesIO()
img.save(img_buffer, format="PNG")
st.image(img_buffer.getvalue(), caption="2D molecular structure", use_container_width=False)

st.subheader("3D molecular structure")
mol3d = Chem.AddHs(mol)
params = AllChem.ETKDGv3()
params.randomSeed = 42
status = AllChem.EmbedMolecule(mol3d, params)
if status == 0:
    try:
        AllChem.MMFFOptimizeMolecule(mol3d)
    except Exception:
        pass
    xyz = Chem.MolToMolBlock(mol3d)
    try:
        import py3Dmol
        view = py3Dmol.view(width=700, height=450)
        view.addModel(xyz, "mol")
        view.setStyle({"stick": {}, "sphere": {"scale": 0.25}})
        view.setBackgroundColor("0xFFFFFF")
        view.zoomTo()
        components.html(view._make_html(), height=470, scrolling=False)
    except Exception:
        st.warning("3D viewer could not load in this environment. The 3D structure was generated, but the interactive viewer is unavailable.")
else:
    st.warning("Could not generate a 3D conformer for this molecule.")

def make_pdf():
    output = io.BytesIO()
    doc = SimpleDocTemplate(output, pagesize=A4, rightMargin=42, leftMargin=42,
                            topMargin=42, bottomMargin=42)
    styles = getSampleStyleSheet()
    story = [
        Paragraph("Molecule Studio – Molecular Analysis Report", styles["Title"]),
        Spacer(1, 12),
        Paragraph("Project: Molecular structure and property analysis", styles["Heading2"]),
        Paragraph(f"<b>Input SMILES:</b> {smiles.strip()}", styles["BodyText"]),
        Paragraph(f"<b>Canonical SMILES:</b> {canonical}", styles["BodyText"]),
        Spacer(1, 12),
        Paragraph("Molecular properties", styles["Heading2"]),
    ]
    data = [
        ["Property", "Value"],
        ["Molecular formula", formula],
        ["Molecular weight", f"{mw:.2f} g/mol"],
        ["Number of atoms", str(atoms)],
        ["Number of bonds", str(bonds)],
        ["LogP", f"{logp:.2f}"],
        ["H-bond donors", str(hbd)],
        ["H-bond acceptors", str(hba)],
    ]
    table = Table(data, colWidths=[2.5*inch, 3.0*inch])
    table.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#243B53")),
        ("TEXTCOLOR", (0,0), (-1,0), colors.white),
        ("GRID", (0,0), (-1,-1), 0.5, colors.grey),
        ("PADDING", (0,0), (-1,-1), 7),
        ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, colors.HexColor("#F1F5F9")]),
    ]))
    story += [table, Spacer(1, 16), Paragraph("2D molecular structure", styles["Heading2"])]
    story.append(PDFImage(io.BytesIO(img_buffer.getvalue()), width=4.5*inch, height=3.1*inch))
    story += [Spacer(1, 12), Paragraph(
        "Conclusion: The entered SMILES was successfully parsed and its 2D structure and "
        "basic molecular properties were generated using RDKit. A 3D conformer was attempted "
        "for visualization. This report is intended for educational project use.",
        styles["BodyText"])]
    doc.build(story)
    return output.getvalue()

st.subheader("Final project report")
st.write("Create a PDF containing the SMILES, molecular properties and 2D structure.")
pdf = make_pdf()
st.download_button("📄 Download molecular report (PDF)", data=pdf,
                   file_name="molecule_analysis_report.pdf",
                   mime="application/pdf", type="primary")

with st.expander("About this project"):
    st.write("Built with Python, Streamlit, RDKit, py3Dmol and ReportLab. "
             "Molecular descriptors are calculated from the entered structure.")
