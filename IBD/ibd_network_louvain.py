#!/usr/bin/env python3 
import pandas as pd 
import networkx as nx 
import matplotlib.pyplot as plt 
import os 
import community as community_louvain 
 
# ========== CONFIGURATION ========== 
chr = "3R"  # Chromosome to analyze 
ibd_file = f"{chr}_phased_data.ibd.gz" 
meta_file = "sample_data.csv" 
output_dir = "./IBD_results" 
min_cm_length = 0.8  # Minimum IBD length in cM 
color_by = 'population'  # Options: 'population' or 'community' 
 
# ========== LOAD IBD AND METADATA ========== 
columns = ['ID1', 'HAP1', 'ID2', 'HAP2', 'CHR', 'START', 'END', 'LOD', 'CM_LEN'] 
ibd_df = pd.read_csv(ibd_file, sep='\t', names=columns, compression='gzip') 
ibd_df = ibd_df[ibd_df['CM_LEN'] >= min_cm_length] 
 
# Sum IBD lengths between each individual pair 
ibd_pairs = ibd_df.groupby(['ID1', 'ID2'])['CM_LEN'].sum().reset_index() 
 
# Load population metadata 
meta_df = pd.read_csv(meta_file) 
meta_dict = dict(zip(meta_df['pop'], meta_df['taxon'])) 
 
# ========== BUILD GRAPH ========== 
G = nx.Graph() 
for _, row in ibd_pairs.iterrows(): 
    G.add_edge(row['ID1'], row['ID2'], weight=row['CM_LEN']) 
 
# Assign metadata 
for node in G.nodes(): 
    G.nodes[node]['population'] = meta_dict.get(node, 'Unknown') 
 
# ========== COLORING SCHEME ========== 
if color_by == 'population': 
    pop_set = sorted(set(nx.get_node_attributes(G, 'population').values())) 
    #cmap = plt.cm.get_cmap("tab10", len(pop_set)) 
    #color_map = {pop: cmap(i) for i, pop in enumerate(pop_set)} 
    color_map = { 
        "gambiae": "blue", 
        "arabiensis": "green", 
        "coluzzii": "orange", 
        "melas": "red", 
        "bissau": "purple", 
    } 
    
    # Custom legend labels with formatting
    legend_label_map = {
        "gambiae": "An. gambiae s.s",
        "arabiensis": "An. arabiensis",
        "coluzzii": "An. coluzzii", 
        "melas": "An. melas",
        "bissau": "Bissau",
    }
    
    node_colors = [color_map[G.nodes[node]['population']] for node in G.nodes()] 
    legend_items = [(legend_label_map.get(pop, pop), color_map[pop], pop) for pop in pop_set] 
    legend_title = "Population" 
    filename = "ibd_network_population.png" 
 
elif color_by == 'community': 
    partition = community_louvain.best_partition(G) 
    communities = sorted(set(partition.values())) 
    cmap = plt.cm.get_cmap("tab20", len(communities)) 
    node_colors = [cmap(partition[node]) for node in G.nodes()] 
    legend_items = [(f"Community {comm}", cmap(comm)) for comm in communities] 
    legend_title = "Louvain Communities" 
    filename = "ibd_network_community.png" 
 
# ========== NODE/EDGE PROPERTIES ========== 
degrees = dict(G.degree()) 
node_sizes = [degrees[node]*10 for node in G.nodes()] 
edge_widths = [G[u][v]['weight'] for u, v in G.edges()] 
 
# Label top-degree nodes 
top_nodes = sorted(G.degree, key=lambda x: x[1], reverse=True)[:5] 
all_labels = {node: node for node, _ in top_nodes}

# Split labels for different styling
# Make the top 2 nodes italic, the rest bold
top_2_nodes = {node: node for node, _ in top_nodes[:2]}  # Top 2 for italic
remaining_nodes = {node: node for node, _ in top_nodes[2:]}  # Rest for bold
 
# ========== PLOT ========== 
os.makedirs(output_dir, exist_ok=True) 
plt.figure(figsize=(14, 12)) 
pos = nx.spring_layout(G, seed=42, k=0.15) 
 
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.9) 
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.4) 


# Legend with custom formatting
legend_handles = []
legend_labels = []

for display_label, color, original_pop in legend_items: 
    handle = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, 
                       markersize=10, linestyle='')
    legend_handles.append(handle)
    legend_labels.append(display_label)

# Create legend with larger, bold text
legend = plt.legend(legend_handles, legend_labels, title=legend_title, loc='upper left', 
                   fontsize=22, title_fontsize=21, bbox_to_anchor=(1, 1))

# Make legend title bold
plt.setp(legend.get_title(), fontweight='bold')

# Apply italic formatting to specific labels
for text, (display_label, color, original_pop) in zip(legend.get_texts(), legend_items):
    text.set_fontweight('bold')
    if original_pop != 'bissau':
        text.set_fontstyle('italic') 
 
plt.title(rf"Identity by Descent Network of $\boldsymbol{{Anopheles\ gambiae}}$ complex on Chromosome {chr}", 
          fontsize=22, fontweight='bold', pad=20)
plt.axis('off') 
plt.tight_layout() 
plt.savefig(f"{output_dir}/{chr}_{filename}", dpi=300) 
plt.close()