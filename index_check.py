#!/usr/bin/env python
# coding: utf-8

# In[1]:

import numpy as np 
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from io import BytesIO
from itertools import combinations

st.title("Index Color Balance")

st.subheader("XLEAP SBS reagents on the NextSeq 1000/2000 and NovaSeq X/X Plus")
st.write("T = Green")
st.write("C = Blue + Green")
st.write("A = Blue")
st.write("G = Dark (no label)")

st.subheader("\n**Checking for:**")
st.write("- Both signals are present in both channels for every cycle (it is acceptable to have only the Green channel from T or C if needed)")
st.write("- Avoid having only a signal from the Blue channel from A or A+G in any cycle")
st.write("- Avoid no signal cycles (only G present)")
st.write("- Either of the first two cycles must start with one base other than G")
st.write("- Ensure indexes are sufficiently diverse to prevent misassignment")
st.write("- Check for over-represented nucleotide biases in I7/I5 pairings")

def generate_template(): 
    template_df = pd.DataFrame ({"Sample Name": [], "I7 ID": [],"I7 Sequence": [],"I5 ID": [],"I5 Sequence": []})
    output = BytesIO()
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        template_df.to_excel(writer, index=False, sheet_name="Index Template")
    output.seek(0)
    return output

st.download_button(
    label="Download Excel Template",
    data=generate_template(),
    file_name="Index_Template.xlsx",
    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

uploaded_file = st.file_uploader("Upload your filled-in template (Excel file)", type=["xlsx"])

def hamming_distance(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def check_index_diversity(indexes, index_names, min_distance=3):
    too_close_pairs = []
    
    # Ensure Hamming distance is checked separately for I7 and I5
    half = len(indexes) // 2  # Assuming equal numbers of I7 and I5 indexes
    i7_indexes = indexes[:half]
    i5_indexes = indexes[half:]
    i7_names = index_names[:half]
    i5_names = index_names[half:]
    
    # Check I7 separately
    for (idx1, seq1), (idx2, seq2) in combinations(enumerate(i7_indexes), 2):
        dist = hamming_distance(seq1, seq2)
        if dist < min_distance:
            too_close_pairs.append((i7_names[idx1], i7_names[idx2], dist))  # Use correct names
    
    # Check I5 separately
    for (idx1, seq1), (idx2, seq2) in combinations(enumerate(i5_indexes), 2):
        dist = hamming_distance(seq1, seq2)
        if dist < min_distance:
            too_close_pairs.append((i5_names[idx1], i5_names[idx2], dist))  # Use correct names
    
    return too_close_pairs

def check_nucleotide_bias(i7_df, i5_df, threshold=60):
    problematic_cycles = []
    for cycle in range(len(i7_df)):
        for base in "ATCG":
            i7_freq = i7_df.iloc[cycle][base]
            i5_freq = i5_df.iloc[cycle][base]
            
            if i7_freq > threshold and i5_freq > threshold:
                problematic_cycles.append((cycle + 1, f"Overrepresented nucleotide '{base}' in both I7 and I5 (> {threshold}%)"))
    return problematic_cycles

def check_color_balance(indexes, index_names):
    index_length = len(indexes[0])
    
    if any(len(idx) != index_length for idx in indexes):
        st.error("All index sequences must be the same length.")
        return None, []

    index_matrix = np.array([list(idx.upper()) for idx in indexes])
    nucleotide_counts = {nuc: np.count_nonzero(index_matrix == nuc, axis=0) for nuc in "ATCG"}
    df = pd.DataFrame(nucleotide_counts)
    df = df.div(df.sum(axis=1), axis=0) * 100  
    
    # Shift cycle positions to start from 1
    df.index = df.index + 1
    
    # Detect potential sequencing issues based on NovaSeq X+ guidelines
    problematic_cycles = []
    for cycle, row in df.iterrows():
        if row['G'] >= 90:  # Cycle-wide check: if all indexes have G
            problematic_cycles.append((cycle, "Dark Cycle: Only G detected (No signal)"))
        elif row['A'] + row['G'] >= 90:  # Cycle-wide check: if all indexes only contain A or G
            problematic_cycles.append((cycle, "Potential Issue: Only A + G detected (Blue channel only)"))
    
    return df, problematic_cycles

if uploaded_file:
    try:
        df = pd.read_excel(uploaded_file, sheet_name="Index Template", dtype=str)
        df.columns = df.columns.str.strip() 
        
        if "I7 Sequence" in df.columns and "I5 Sequence" in df.columns:
            df.columns = df.columns.str.strip()  # Ensure no trailing spaces in headers
            df["I7 Sequence"] = df["I7 Sequence"].str.strip()
            df["I5 Sequence"] = df["I5 Sequence"].str.strip()
        
            i7_indexes = df["I7 Sequence"].dropna().astype(str).tolist()
            i5_indexes = df["I5 Sequence"].dropna().astype(str).tolist()
            i7_ids = df["I7 ID"].astype(str).tolist()
            i5_ids = df["I5 ID"].astype(str).tolist()
            
            st.subheader("I7 Index")
            color_balance_df_i7, problematic_cycles_i7 = check_color_balance(i7_indexes, i7_ids)
            if color_balance_df_i7 is not None:
                st.dataframe(color_balance_df_i7)
                
                fig, ax = plt.subplots()
                color_balance_df_i7.plot(kind="bar", stacked=True, ax=ax, color=["blue", "green", "cyan", "gray"])
                ax.set_xlabel("Cycle Position")
                ax.set_ylabel("Nucleotide %")
                ax.set_title("I7 Nucleotide % Per Cycle")
                ax.set_xticks(range(len(color_balance_df_i7)))
                ax.set_xticklabels(range(1, len(color_balance_df_i7) + 1))
                st.pyplot(fig)
                
                if problematic_cycles_i7:
                    st.warning("⚠️ Potential sequencing issues detected in I7 Index:")
                    for cycle, issue in problematic_cycles_i7:
                        st.write(f"Cycle {cycle}: {issue}")
            
            st.subheader("I5 Index")
            color_balance_df_i5, problematic_cycles_i5 = check_color_balance(i5_indexes, i5_ids)
            if color_balance_df_i5 is not None:
                st.dataframe(color_balance_df_i5)
                
                fig, ax = plt.subplots()
                color_balance_df_i5.plot(kind="bar", stacked=True, ax=ax, color=["blue", "green", "cyan", "gray"])
                ax.set_xlabel("Cycle Position")
                ax.set_ylabel("Nucleotide %")
                ax.set_title("I5 Nucleotide % Per Cycle")
                ax.set_xticks(range(len(color_balance_df_i5)))
                ax.set_xticklabels(range(1, len(color_balance_df_i5) + 1))
                st.pyplot(fig)
                
                if problematic_cycles_i5:
                    st.warning("⚠️ Potential sequencing issues detected in I5 Index:")
                    for cycle, issue in problematic_cycles_i5:
                        st.write(f"Cycle {cycle}: {issue}")
            
            # Check for overrepresented nucleotide bias in I7/I5 pairs
            bias_issues = check_nucleotide_bias(color_balance_df_i7, color_balance_df_i5)
            if bias_issues:
                st.warning("⚠️ Overrepresented Nucleotide Bias Detected:")
                for cycle, issue in bias_issues:
                    st.write(f"Cycle {cycle}: {issue}")
            
            # Check index diversity
            diversity_issues = check_index_diversity(i7_indexes + i5_indexes, i7_ids + i5_ids)
            if diversity_issues:
                st.warning("⚠️ Index diversity warning: Some index pairs are too similar:")
                for idx1, idx2, distance in diversity_issues:
                    st.write(f"Indexes {idx1} and {idx2} have a Hamming distance of {distance} (too close)")
    
    except Exception as e:
        st.error(f"Error processing file: {e}")


# In[ ]:




