import streamlit as st

st.set_page_config(page_title="Mol* Demo — STICkit", layout="wide")
st.title("Mol* (WebGL) Viewer Demo")

pdb_id = st.text_input("PDB ID", "1CRN")
height = st.number_input("Viewer height (px)", 600, 1600, 800, step=50)

# Embed Mol* via CDN using HTML component
molstar_html = f'''
<div id="app" style="width:100%; height:{height}px;"></div>
<script type="module">
  import {{ Viewer }} from "https://unpkg.com/molstar@latest/build/viewer/molstar.js";
  const viewer = new Viewer("app", {{ layoutIsExpanded: true, layoutShowControls: true }});
  viewer.loadPdb(`https://files.rcsb.org/view/{pdb_id}.pdb`);
</script>
'''
st.components.v1.html(molstar_html, height=height + 40, scrolling=False)
