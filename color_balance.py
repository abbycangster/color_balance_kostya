#!/usr/bin/env python
# coding: utf-8

# In[1]:

import numpy as np 
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from io import BytesIO

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

def check_color_balance(indexes):
    index_length = len(indexes[0])
    
    if any(len(idx) != index_length for idx in indexes):
        st.error("All index sequences must be the same length.")
        return None

    index_matrix = np.array([list(idx.upper()) for idx in indexes])
    nucleotide_counts = {nuc: np.count_nonzero(index_matrix == nuc, axis=0) for nuc in "ATCG"}
    df = pd.DataFrame(nucleotide_counts)
    df = df.div(df.sum(axis=1), axis=0) * 100  
    
    # Shift cycle positions to start from 1
    df.index = df.index + 1
    
    # Detect potential sequencing issues based on Illumina guidelines
    problematic_cycles = []
    for cycle, row in df.iterrows():
        if row['G'] == 100:
            problematic_cycles.append((cycle, "Dark Cycle: Only G detected (No signal)"))
        elif row['A'] + row['G'] == 100:
            problematic_cycles.append((cycle, "Potential Issue: Only A + G detected (Blue channel only)"))
    
    # Check that either of the first two cycles contains a base other than G
    if df.iloc[0]['G'] == 100 and df.iloc[1]['G'] == 100:
        problematic_cycles.append(("1-2", "Warning: First two cycles contain only G (No signal at start)"))
    
    return df, problematic_cycles

if uploaded_file:
    try:
        df = pd.read_excel(uploaded_file, sheet_name="Index Template")
        
        if "I7 Sequence" in df.columns and "I5 Sequence" in df.columns:
            i7_indexes = df["I7 Sequence"].dropna().astype(str).tolist()
            i5_indexes = df["I5 Sequence"].dropna().astype(str).tolist()
            
            st.subheader("I7 Index")
            color_balance_df_i7, problematic_cycles_i7 = check_color_balance(i7_indexes)
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
            color_balance_df_i5, problematic_cycles_i5 = check_color_balance(i5_indexes)
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
            
    except Exception as e:
        st.error(f"Error processing file: {e}")

# In[ ]:




