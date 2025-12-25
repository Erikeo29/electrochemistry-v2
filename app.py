import streamlit as st
import base64
import os
import pandas as pd

# --- Configuration de la page ---
st.set_page_config(page_title="Electrochemistry Simulation Platform", layout="wide")

# --- Styles CSS personnalisés (Look Académique/Latex) ---
custom_css = """
<style>
/* Masquer éléments Streamlit inutiles */
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
header {visibility: hidden;}
.stDeployButton {display: none;}
[data-testid="stToolbar"] {display: none;}

/* Style des onglets */
.stTabs [data-baseweb="tab-list"] {
    gap: 8px;
}
.stTabs [data-baseweb="tab"] {
    height: 50px;
    padding: 10px 24px;
    background-color: #f0f2f6;
    border-radius: 8px 8px 0 0;
    font-weight: 600;
    font-size: 16px;
}
.stTabs [aria-selected="true"] {
    background-color: #004b87;
    color: white;
}

/* Bouton retour en haut */
.back-to-top {
    position: fixed;
    bottom: 30px;
    right: 30px;
    z-index: 9999;
    background-color: #004b87;
    color: white;
    border: none;
    border-radius: 50%;
    width: 50px;
    height: 50px;
    font-size: 24px;
    cursor: pointer;
    box-shadow: 0 4px 8px rgba(0,0,0,0.3);
    display: flex;
    align-items: center;
    justify-content: center;
    text-decoration: none;
}
.back-to-top:hover {
    background-color: #003366;
    transform: scale(1.1);
}

/* Style scientifique (Latex-like) */
.main {
    font-family: 'serif';
}
h1, h2, h3 {
    color: #004b87;
}
</style>

<!-- Bouton retour en haut -->
<a href="#top" class="back-to-top" title="Retour en haut">&#8679;</a>
<div id="top"></div>
"""
st.markdown(custom_css, unsafe_allow_html=True)

# --- Chemins ---
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
DOC_PATH = os.path.join(ROOT_DIR, "docs")
ASSETS_PATH = os.path.join(ROOT_DIR, "assets")

# --- Fonctions Utilitaires ---
def load_file_content(path):
    try:
        with open(path, 'r', encoding='utf-8') as f: return f.read()
    except Exception: return f"Document en cours de rédaction : {os.path.basename(path)}"

def display_smart_markdown(content):
    if "```python" in content:
        parts = content.split("```python")
        for i, part in enumerate(parts):
            if i > 0:
                if "```" in part:
                    code, text = part.split("```", 1)
                    st.code(code.strip(), language='python')
                    if text.strip(): st.markdown(text)
                else:
                    st.code(part.strip(), language='python')
            elif part.strip():
                st.markdown(part)
    elif "```cpp" in content:
        parts = content.split("```cpp")
        for i, part in enumerate(parts):
            if i > 0:
                if "```" in part:
                    code, text = part.split("```", 1)
                    st.code(code.strip(), language='cpp')
                    if text.strip(): st.markdown(text)
                else:
                    st.code(part.strip(), language='cpp')
            elif part.strip():
                st.markdown(part)
    else:
        st.markdown(content)

# --- Barre Latérale ---
st.sidebar.title("Electrochemistry Simulation")

# Section Accueil
st.sidebar.subheader("Général")
main_nav = st.sidebar.radio("Navigation principale", ["Accueil"])

# Section Cyclic Voltammetry
st.sidebar.markdown("---")
st.sidebar.subheader("Cyclic Voltammetry (CV)")
cv_nav = st.sidebar.selectbox("Modules CV", 
    ["Introduction", "Comparaison des modèles", "Python (Firedrake)", "OpenFOAM", "Conclusion", "Équations clés", "Lexique", "Un peu d'histoire"],
    index=None, placeholder="Sélectionnez un module...")

# Section Electroplating
st.sidebar.markdown("---")
st.sidebar.subheader("Electroplating")
plating_nav = st.sidebar.selectbox("Modules Plating", 
    ["Introduction", "Python (Antigravity)", "Conclusion", "Équations clés", "Lexique", "Un peu d'histoire"],
    index=None, placeholder="Sélectionnez un module...")

st.sidebar.markdown("---")
st.sidebar.markdown("**Version 1.0.0**")

# --- Logique de Navigation ---
if main_nav == "Accueil" and cv_nav is None and plating_nav is None:
    st.title("Plateforme de Simulation Électrochimique")
    st.markdown(load_file_content(os.path.join(DOC_PATH, "accueil/accueil.md")))
    
    col1, col2 = st.columns(2)
    with col1:
        st.info("### Cyclic Voltammetry")
        res_img = os.path.join(ASSETS_PATH, "cv/png/cv_result_example.png")
        if os.path.exists(res_img):
            st.image(res_img, use_container_width=True)
        st.write("Modélisation du transport de masse et de la cinétique aux électrodes par éléments finis.")
        
    with col2:
        st.success("### Electroplating")
        res_img = os.path.join(ASSETS_PATH, "plating/results/plating_result_example.png")
        if os.path.exists(res_img):
            st.image(res_img, use_container_width=True)
        st.write("Simulation de dépôt électrolytique et distribution de courant secondaire.")

elif cv_nav:
    st.title(f"CV : {cv_nav}")
    if cv_nav == "Introduction":
        display_smart_markdown(load_file_content(os.path.join(DOC_PATH, "intro/intro_cv.md")))
    elif cv_nav == "Comparaison des modèles":
        display_smart_markdown(load_file_content(os.path.join(DOC_PATH, "intro/comparaison_cv.md")))
    elif cv_nav == "Python (Firedrake)":
        tab1, tab2 = st.tabs(["Physique & Résultats", "Code Source"])
        with tab1:
            display_smart_markdown(load_file_content(os.path.join(DOC_PATH, "physics/cv_firedrake.md")))
            res_img = os.path.join(ASSETS_PATH, "cv/png/cv_result_example.png")
            if os.path.exists(res_img):
                st.image(res_img, caption="Exemple de Voltammogramme (Firedrake)", use_container_width=True)
        with tab2:
            display_smart_markdown(load_file_content(os.path.join(DOC_PATH, "code/cv_firedrake_code.md")))
    elif cv_nav == "OpenFOAM":
        tab1, tab2 = st.tabs(["Physique & Résultats", "Configuration"])
        with tab1:
            display_smart_markdown(load_file_content(os.path.join(DOC_PATH, "physics/cv_openfoam.md")))
        with tab2:
            display_smart_markdown(load_file_content(os.path.join(DOC_PATH, "code/cv_openfoam_code.md")))
    elif cv_nav == "Conclusion":
        display_smart_markdown(load_file_content(os.path.join(DOC_PATH, "conclusion/cv_conclusion.md")))
    elif cv_nav == "Équations clés":
        display_smart_markdown(load_file_content(os.path.join(DOC_PATH, "equations/cv_equations.md")))
    elif cv_nav == "Lexique":
        display_smart_markdown(load_file_content(os.path.join(DOC_PATH, "lexique/cv_lexique.md")))
    elif cv_nav == "Un peu d'histoire":
        display_smart_markdown(load_file_content(os.path.join(DOC_PATH, "histoire/cv_histoire.md")))

elif plating_nav:
    st.title(f"Electroplating : {plating_nav}")
    if plating_nav == "Introduction":
        display_smart_markdown(load_file_content(os.path.join(DOC_PATH, "intro/intro_plating.md")))
    elif plating_nav == "Python (Antigravity)":
        tab1, tab2 = st.tabs(["Physique & Résultats", "Code Source"])
        with tab1:
            display_smart_markdown(load_file_content(os.path.join(DOC_PATH, "physics/plating_antigravity.md")))
            res_img = os.path.join(ASSETS_PATH, "plating/results/plating_result_example.png")
            if os.path.exists(res_img):
                st.image(res_img, caption="Distribution d'épaisseur du dépôt (Antigravity)", use_container_width=True)
        with tab2:
            display_smart_markdown(load_file_content(os.path.join(DOC_PATH, "code/plating_antigravity_code.md")))
    elif plating_nav == "Conclusion":
        display_smart_markdown(load_file_content(os.path.join(DOC_PATH, "conclusion/plating_conclusion.md")))
    elif plating_nav == "Équations clés":
        display_smart_markdown(load_file_content(os.path.join(DOC_PATH, "equations/plating_equations.md")))
    elif plating_nav == "Lexique":
        display_smart_markdown(load_file_content(os.path.join(DOC_PATH, "lexique/plating_lexique.md")))
    elif plating_nav == "Un peu d'histoire":
        display_smart_markdown(load_file_content(os.path.join(DOC_PATH, "histoire/plating_histoire.md")))