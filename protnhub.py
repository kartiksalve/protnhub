import streamlit as st
import requests
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mpl
import openai

# ------------- Configuration -------------
st.set_page_config(page_title="ProtHub", layout="wide")
openai.api_key = st.secrets["openai"]["api_key"]

# ------------- Constants -----------------
STRING_API_URL = "https://string-db.org/api"
STRING_OUTPUT_FORMAT = "json"
STRING_METHOD = "network"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

# ------------- Helper Functions -----------------
def map_sequence_to_uniprot(input_text):
    lines = input_text.strip().splitlines()
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    if not sequence:
        return None
    params = {"query": f'sequence:"{sequence}"', "format": "json", "size": 1}
    r = requests.get(UNIPROT_SEARCH_URL, params=params)
    if r.status_code == 200 and r.json().get("results"):
        return r.json()["results"][0]["primaryAccession"]
    return None

def get_string_interactions(uniprot_id, species=9606, min_score=0.4):
    params = {
        "identifiers": uniprot_id,
        "species": species,
        "caller_identity": "streamlit_app",
        "required_score": int(min_score * 1000)
    }
    r = requests.post(f"{STRING_API_URL}/{STRING_OUTPUT_FORMAT}/{STRING_METHOD}", data=params)
    return r.json() if r.status_code == 200 else None

def build_network(data):
    G = nx.DiGraph()
    for i in data:
        G.add_edge(i["preferredName_A"], i["preferredName_B"], weight=i["score"])
    return G

def find_hub_genes(G, top_n=5):
    degrees = dict(G.degree())
    sorted_genes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)
    return [gene for gene, _ in sorted_genes[:top_n]]

def visualize_network(G, hub_genes):
    seed = 13648
    pos = nx.spring_layout(G, seed=seed)
    node_sizes = [400 + 100 * G.degree(n) for n in G.nodes()]
    M = G.number_of_edges()
    edge_colors = range(2, M + 2)
    edge_alphas = [(5 + i) / (M + 4) for i in range(M)]
    cmap = plt.cm.plasma

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_title("Interaction Network", fontsize=14)
    nodes = nx.draw_networkx_nodes(
        G, pos, ax=ax, node_size=node_sizes,
        node_color=["gold" if n in hub_genes else "skyblue" for n in G.nodes()]
    )
    edges = nx.draw_networkx_edges(
        G, pos, ax=ax, edge_color=edge_colors, edge_cmap=cmap,
        arrowstyle="->", arrowsize=10, width=2
    )
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=10, font_color="black")
    for i, e in enumerate(edges):
        e.set_alpha(edge_alphas[i])
    pc = mpl.collections.PatchCollection(edges, cmap=cmap)
    pc.set_array(edge_colors)
    ax.set_axis_off()
    plt.colorbar(pc, ax=ax)
    st.pyplot(fig)

def explain_hub_genes(hub_genes):
    if not hub_genes:
        return "No hub genes found."
    prompt = (
        f"Explain the biological significance and potential functions of these hub genes in a protein-protein interaction network:\n"
        f"{', '.join(hub_genes)}"
    )
    try:
        response = openai.ChatCompletion.create(
            model="gpt-4",
            messages=[
                {"role": "system", "content": "You are a helpful bioinformatics assistant."},
                {"role": "user", "content": prompt}
            ],
            max_tokens=500,
            temperature=0.7
        )
        return response["choices"][0]["message"]["content"].strip()
    except Exception as e:
        return f"GeneChat error: {e}"

# ------------- Streamlit UI -----------------
st.title("🔬 ProtHub: Protein Hub Explorer + GeneChat")

input_type = st.radio("Input Type", ["Protein Name / UniProt ID", "Raw Sequence"])
user_input = st.text_area("Paste your input here:")

species_dict = {
    "Human (Homo sapiens)": 9606,
    "Mouse (Mus musculus)": 10090,
    "Rat (Rattus norvegicus)": 10116,
    "Zebrafish (Danio rerio)": 7955,
    "Fruit fly (Drosophila melanogaster)": 7227,
    "Custom (enter manually)": None,
}
selected_species = st.selectbox("Choose species", list(species_dict.keys()))
species = st.number_input("Enter NCBI Taxonomy ID:", value=9606) if selected_species == "Custom (enter manually)" else species_dict[selected_species]
score_threshold = st.slider("Minimum Interaction Score", 0.0, 1.0, 0.4, 0.05)

if st.button("Analyze Network"):
    with st.spinner("Fetching and building network..."):
        if input_type == "Raw Sequence":
            uniprot_id = map_sequence_to_uniprot(user_input.strip())
            if not uniprot_id:
                st.error("Could not map the sequence to a UniProt ID.")
                st.stop()
            else:
                st.success(f"Mapped to UniProt ID: {uniprot_id}")
        else:
            uniprot_id = user_input.strip()

        data = get_string_interactions(uniprot_id, species, score_threshold)
        if not data:
            st.error("No interaction data found.")
        else:
            G = build_network(data)
            hub_genes = find_hub_genes(G)
            st.subheader("🔗 Top Hub Genes")
            st.success(", ".join(hub_genes))
            st.subheader("📊 Network Visualization")
            visualize_network(G, hub_genes)
            st.subheader("🧠 GeneChat: Explain These Hub Genes")
            st.info(explain_hub_genes(hub_genes))