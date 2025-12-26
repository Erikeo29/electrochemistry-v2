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
    bottom: 50%;
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
"")
st.markdown(custom_css, unsafe_allow_html=True)

# --- Chemins ---
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
DOC_PATH = os.path.join(ROOT_DIR, "docs")
ASSETS_PATH = os.path.join(ROOT_DIR, "assets")

# --- Dictionnaire de Traduction UI ---
TRANSLATIONS = {
    "fr": {
        "title": "Plateforme de Simulation Électrochimique",
        "sidebar_title": "Simulation Électrochimique",
        "gen_header": "Général",
        "gen_home": "Accueil",
        "cv_header": "Voltamétrie Cyclique (CV)",
        "cv_modules": [
            "Introduction", "Comparaison des modèles", "Python (Firedrake)", 
            "OpenFOAM", "Conclusion", "Équations clés", "Lexique", "Un peu d'histoire"
        ],
        "plating_header": "Électrodéposition",
        "plating_modules": [
            "Introduction", "Python (Antigravity)", "Conclusion", 
            "Équations clés", "Lexique", "Un peu d'histoire"
        ],
        "version_info": """**Version 1.0.1**\nDec 2025\n*EQU*\n\n**Nouveautés :**\n- Documentation code enrichie (imports, détails)\n- Support bilingue FR/EN\n\n**Précédemment (1.0.0) :**\n- Comparaison FEM vs FVM\n- Simulation Galvanostatique""",
        "tabs_cv_python": ["Physique", "Code", "Exemples GIF", "Exemples PNG"],
        "tabs_cv_openfoam": ["Physique", "Configuration", "Résultats"],
        "tabs_plating": ["Physique", "Code", "Exemples GIF", "Exemples PNG"],
        "card_cv_title": "### Voltamétrie Cyclique",
        "card_cv_text": "Modélisation du transport de masse et de la cinétique aux électrodes par éléments finis.",
        "card_plating_title": "### Électrodéposition",
        "card_plating_text": "Simulation de dépôt électrolytique et distribution de courant secondaire.",
        "gif_coming_soon": "Visualisation dynamique (Gifs) - À venir",
        "no_gif": "Aucune animation GIF disponible pour le moment.",
        "png_results": "Résultats Graphiques",
        "png_thickness": "Cartes d'épaisseur",
        "openfoam_soon": "Résultats OpenFOAM - À venir"
    },
    "en": {
        "title": "Electrochemical Simulation Platform",
        "sidebar_title": "Electrochemical Simulation",
        "gen_header": "General",
        "gen_home": "Home",
        "cv_header": "Cyclic Voltammetry (CV)",
        "cv_modules": [
            "Introduction", "Technical Comparison", "Python (Firedrake)", 
            "OpenFOAM", "Conclusion", "Key Equations", "Glossary", "A Bit of History"
        ],
        "plating_header": "Electroplating",
        "plating_modules": [
            "Introduction", "Python (Antigravity)", "Conclusion", 
            "Key Equations", "Glossary", "A Bit of History"
        ],
        "version_info": """**Version 1.0.1**\nDec 2025\n*EQU*\n\n**New Features:**\n- Enriched code documentation\n- Bilingual support FR/EN\n\n**Previously (1.0.0):**\n- FEM vs FVM Comparison\n- Galvanostatic Simulation""",
        "tabs_cv_python": ["Physics", "Code", "GIF Examples", "PNG Examples"],
        "tabs_cv_openfoam": ["Physics", "Configuration", "Results"],
        "tabs_plating": ["Physics", "Code", "GIF Examples", "PNG Examples"],
        "card_cv_title": "### Cyclic Voltammetry",
        "card_cv_text": "Modeling of mass transport and electrode kinetics using finite elements.",
        "card_plating_title": "### Electroplating",
        "card_plating_text": "Simulation of electrolytic deposition and secondary current distribution.",
        "gif_coming_soon": "Dynamic Visualization (Gifs) - Coming Soon",
        "no_gif": "No GIF animation available at the moment.",
        "png_results": "Graphical Results",
        "png_thickness": "Thickness Maps",
        "openfoam_soon": "OpenFOAM Results - Coming Soon"
    }
}

# --- Fonctions Utilitaires ---
def get_language():
    if 'lang' not in st.session_state:
        st.session_state.lang = 'fr'
    return st.session_state.lang

def t(key):
    lang = get_language()
    return TRANSLATIONS[lang].get(key, key)

def load_file_content(relative_path):
    lang = get_language()
    full_path = os.path.join(DOC_PATH, lang, relative_path)
    try:
        with open(full_path, 'r', encoding='utf-8') as f: return f.read()
    except Exception: return f"Document not found / Document non trouvé : {os.path.join(lang, relative_path)}"

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

def search_images(base_path, extensions=['.png', '.jpg', '.jpeg']):
    images = []
    if not os.path.exists(base_path): return images
    for root, dirs, files in os.walk(base_path):
        for file in files:
            if any(file.lower().endswith(ext) for ext in extensions):
                images.append(os.path.join(root, file))
    return images

# --- Barre Latérale ---

# Sélecteur de langue
col_l1, col_l2 = st.sidebar.columns(2)
lang_selection = st.sidebar.radio(
    "Language",
    ["🇫🇷 FR", "🇬🇧 EN"],
    horizontal=True,
    label_visibility="collapsed"
)

if "FR" in lang_selection:
    st.session_state.lang = "fr"
else:
    st.session_state.lang = "en"

st.sidebar.title(t("sidebar_title"))

# Callbacks
def on_change_gen():
    st.session_state.nav_cv = None
    st.session_state.nav_plating = None

def on_change_cv():
    st.session_state.nav_gen = None
    st.session_state.nav_plating = None

def on_change_plating():
    st.session_state.nav_gen = None
    st.session_state.nav_cv = None

# Init states
if 'nav_gen' not in st.session_state: st.session_state.nav_gen = t("gen_home")
if 'nav_cv' not in st.session_state: st.session_state.nav_cv = None
if 'nav_plating' not in st.session_state: st.session_state.nav_plating = None

# Navigation
st.sidebar.subheader(t("gen_header"))
main_nav = st.sidebar.radio(
    "Nav Gen", 
    [t("gen_home")],
    key="nav_gen", 
    on_change=on_change_gen,
    label_visibility="collapsed"
)

st.sidebar.markdown("---")
st.sidebar.subheader(t("cv_header"))
cv_nav = st.sidebar.radio(
    "Nav CV", 
    t("cv_modules"),
    key="nav_cv",
    index=None,
    on_change=on_change_cv,
    label_visibility="collapsed"
)

st.sidebar.markdown("---")
st.sidebar.subheader(t("plating_header"))
plating_nav = st.sidebar.radio(
    "Nav Plating", 
    t("plating_modules"),
    key="nav_plating",
    index=None,
    on_change=on_change_plating,
    label_visibility="collapsed"
)

st.sidebar.markdown("---")
st.sidebar.markdown(t("version_info"))

# --- Mapping Modules (FR/EN correspondence) ---
# Helper to check selection regardless of language
def is_sel(selection, key_fr, key_en):
    if not selection: return False
    return selection == key_fr or selection == key_en

modules_cv = TRANSLATIONS["fr"]["cv_modules"]
modules_cv_en = TRANSLATIONS["en"]["cv_modules"]
modules_pl = TRANSLATIONS["fr"]["plating_modules"]
modules_pl_en = TRANSLATIONS["en"]["plating_modules"]

# --- Contenu Principal ---

# ACCUEIL
if main_nav == t("gen_home") and cv_nav is None and plating_nav is None:
    st.title(t("title"))
    st.markdown(load_file_content("accueil/accueil.md"))
    
    col1, col2 = st.columns(2)
    with col1:
        st.info(t("card_cv_title"))
        res_img = os.path.join(ASSETS_PATH, "cv/png/cv_result_example.png")
        if os.path.exists(res_img):
            st.image(res_img, use_container_width=True)
        st.write(t("card_cv_text"))
        
    with col2:
        st.success(t("card_plating_title"))
        res_img = os.path.join(ASSETS_PATH, "plating/results/plating_result_example.png")
        if os.path.exists(res_img):
            st.image(res_img, use_container_width=True)
        st.write(t("card_plating_text"))

# CV PAGES
elif cv_nav:
    st.title(f"CV : {cv_nav}")
    
    # Identify index to map FR logic
    try:
        idx = modules_cv.index(cv_nav)
    except ValueError:
        idx = modules_cv_en.index(cv_nav)
    
    # 0: Intro
    if idx == 0:
        display_smart_markdown(load_file_content("intro/intro_cv.md"))
    # 1: Comparaison
    elif idx == 1:
        display_smart_markdown(load_file_content("intro/comparaison_cv.md"))
    # 2: Python (Firedrake)
    elif idx == 2:
        tabs = st.tabs(t("tabs_cv_python"))
        with tabs[0]:
            display_smart_markdown(load_file_content("physics/cv_firedrake.md"))
        with tabs[1]:
            display_smart_markdown(load_file_content("code/cv_firedrake_code.md"))
        with tabs[2]:
            st.info(t("gif_coming_soon"))
            gifs = search_images(os.path.join(ASSETS_PATH, "cv/gif"), ['.gif'])
            if gifs:
                for gif in gifs:
                    st.image(gif, caption=os.path.basename(gif), use_container_width=True)
            else:
                st.warning(t("no_gif"))
        with tabs[3]:
            st.subheader(t("png_results"))
            res_img = os.path.join(ASSETS_PATH, "cv/png/cv_result_example.png")
            if os.path.exists(res_img):
                st.image(res_img, caption="Voltammogramme Firedrake", use_container_width=True)
    # 3: OpenFOAM
    elif idx == 3:
        tabs = st.tabs(t("tabs_cv_openfoam"))
        with tabs[0]:
            display_smart_markdown(load_file_content("physics/cv_openfoam.md"))
        with tabs[1]:
            display_smart_markdown(load_file_content("code/cv_openfoam_code.md"))
        with tabs[2]:
            st.info(t("openfoam_soon"))
    # 4: Conclusion
    elif idx == 4:
        display_smart_markdown(load_file_content("conclusion/cv_conclusion.md"))
    # 5: Equations
    elif idx == 5:
        display_smart_markdown(load_file_content("equations/cv_equations.md"))
    # 6: Lexique
    elif idx == 6:
        display_smart_markdown(load_file_content("lexique/cv_lexique.md"))
    # 7: Histoire
    elif idx == 7:
        display_smart_markdown(load_file_content("histoire/cv_histoire.md"))

# PLATING PAGES
elif plating_nav:
    st.title(f"Plating : {plating_nav}")
    
    try:
        idx = modules_pl.index(plating_nav)
    except ValueError:
        idx = modules_pl_en.index(plating_nav)

    # 0: Intro
    if idx == 0:
        display_smart_markdown(load_file_content("intro/intro_plating.md"))
    # 1: Python (Antigravity)
    elif idx == 1:
        tabs = st.tabs(t("tabs_plating"))
        with tabs[0]:
            display_smart_markdown(load_file_content("physics/plating_antigravity.md"))
        with tabs[1]:
            display_smart_markdown(load_file_content("code/plating_antigravity_code.md")
        with tabs[2]:
             st.info(t("gif_coming_soon"))
             gifs = search_images(os.path.join(ASSETS_PATH, "plating/results"), ['.gif'])
             if gifs:
                for gif in gifs:
                    st.image(gif, caption=os.path.basename(gif), use_container_width=True)
             else:
                st.warning(t("no_gif"))
        with tabs[3]:
            st.subheader(t("png_thickness"))
            res_img = os.path.join(ASSETS_PATH, "plating/results/plating_result_example.png")
            if os.path.exists(res_img):
                st.image(res_img, caption="Thickness Distribution (Antigravity)", use_container_width=True)
    # 2: Conclusion
    elif idx == 2:
        display_smart_markdown(load_file_content("conclusion/plating_conclusion.md"))
    # 3: Equations
    elif idx == 3:
        display_smart_markdown(load_file_content("equations/plating_equations.md"))
    # 4: Lexique
    elif idx == 4:
        display_smart_markdown(load_file_content("lexique/plating_lexique.md"))
    # 5: Histoire
    elif idx == 5:
        display_smart_markdown(load_file_content("histoire/plating_histoire.md"))