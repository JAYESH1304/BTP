import os
os.system("pip install streamlit pandas xlsxwriter openpyxl pymongo matplotlib seaborn")

import streamlit as st
import pandas as pd
import xlsxwriter
from io import BytesIO
from collections import defaultdict
import hashlib
import matplotlib.pyplot as plt
import seaborn as sns

# MongoDB Setup
try:
    from pymongo import MongoClient
    client = MongoClient("mongodb+srv://dhruvmangroliya:Eussmh5MbCBIkLJ6@cluster0.rrnbxfw.mongodb.net/BTP_DB?retryWrites=true&w=majority")
    db = client['BTP_DB']
    results_collection = db['protein_results']
except:
    results_collection = None

# Utility Functions
def is_homo_repeat(s):
    return all(c == s[0] for c in s)

def hash_sequence(sequence):
    return hashlib.md5(sequence.encode()).hexdigest()

@st.cache_data(show_spinner=False)
def fragment_protein_sequence(sequence, max_length=1000):
    return [sequence[i:i+max_length] for i in range(0, len(sequence), max_length)]

def find_homorepeats(protein):
    n = len(protein)
    freq = defaultdict(int)
    i = 0
    while i < n:
        curr = protein[i]
        repeat = ""
        while i < n and curr == protein[i]:
            repeat += protein[i]
            i += 1
        if len(repeat) > 1:
            freq[repeat] += 1
    return freq

def find_hetero_amino_acid_repeats(sequence):
    repeat_counts = defaultdict(int)
    for length in range(2, len(sequence) + 1):
        for i in range(len(sequence) - length + 1):
            substring = sequence[i:i+length]
            repeat_counts[substring] += 1
    return {k: v for k, v in repeat_counts.items() if v > 1}

def check_boundary_repeats(fragments, final_repeats, overlap=50):
    for i in range(len(fragments) - 1):
        left_overlap = fragments[i][-overlap:]
        right_overlap = fragments[i + 1][:overlap]
        overlap_region = left_overlap + right_overlap
        boundary_repeats = find_hetero_amino_acid_repeats(overlap_region)
        for substring, count in boundary_repeats.items():
            if any(aa in left_overlap for aa in substring) and any(aa in right_overlap for aa in substring):
                final_repeats[substring] += count
    return final_repeats

def find_new_boundary_repeats(fragments, final_repeats, overlap=50):
    new_repeats = defaultdict(int)
    for i in range(len(fragments) - 1):
        left_overlap = fragments[i][-overlap:]
        right_overlap = fragments[i + 1][:overlap]
        overlap_region = left_overlap + right_overlap
        boundary_repeats = find_hetero_amino_acid_repeats(overlap_region)
        for substring, count in boundary_repeats.items():
            if any(aa in left_overlap for aa in substring) and any(aa in right_overlap for aa in substring):
                if substring not in final_repeats:
                    new_repeats[substring] += count
    return new_repeats

def get_or_process_sequence(sequence, analysis_type, overlap=50):
    if results_collection is None:
        return {}
    hash_input = f"{sequence}_{analysis_type}"
    sequence_hash = hash_sequence(hash_input)
    cached = results_collection.find_one({"_id": sequence_hash})
    if cached:
        return cached["repeats"]

    fragments = fragment_protein_sequence(sequence)
    final_repeats = defaultdict(int)

    if analysis_type == "Hetero":
        for fragment in fragments:
            fragment_repeats = find_hetero_amino_acid_repeats(fragment)
            for k, v in fragment_repeats.items():
                final_repeats[k] += v
        final_repeats = check_boundary_repeats(fragments, final_repeats, overlap)
        new_repeats = find_new_boundary_repeats(fragments, final_repeats, overlap)
        for k, v in new_repeats.items():
            final_repeats[k] += v
        final_repeats = {k: v for k, v in final_repeats.items() if not is_homo_repeat(k)}

    elif analysis_type == "Homo":
        final_repeats = find_homorepeats(sequence)

    elif analysis_type == "Both":
        hetero_repeats = defaultdict(int)
        for fragment in fragments:
            fragment_repeats = find_hetero_amino_acid_repeats(fragment)
            for k, v in fragment_repeats.items():
                hetero_repeats[k] += v
        hetero_repeats = check_boundary_repeats(fragments, hetero_repeats, overlap)
        new_repeats = find_new_boundary_repeats(fragments, hetero_repeats, overlap)
        for k, v in new_repeats.items():
            hetero_repeats[k] += v
        hetero_repeats = {k: v for k, v in hetero_repeats.items() if not is_homo_repeat(k)}
        homo_repeats = find_homorepeats(sequence)
        final_repeats = homo_repeats.copy()
        for k, v in hetero_repeats.items():
            final_repeats[k] += v

    results_collection.insert_one({
        "_id": sequence_hash,
        "sequence": sequence,
        "analysis_type": analysis_type,
        "repeats": dict(final_repeats)
    })
    return final_repeats

def process_excel(excel_data, analysis_type):
    repeats = set()
    sequence_data = []
    count = 0
    for sheet_name in excel_data.sheet_names:
        df = excel_data.parse(sheet_name)
        if len(df.columns) < 3:
            st.error(f"Error: The sheet '{sheet_name}' must have at least three columns: ID, Protein Name, Sequence")
            return None, None
        for _, row in df.iterrows():
            entry_id = str(row[0])
            protein_name = str(row[1])
            sequence = str(row[2]).replace('"', '').replace(' ', '').strip()
            if not sequence:
                continue
            count += 1
            freq = get_or_process_sequence(sequence, analysis_type)
            sequence_data.append((entry_id, protein_name, freq))
            repeats.update(freq.keys())
    st.toast(f"{count} sequences processed.")
    return repeats, sequence_data

def create_excel(sequences_data, repeats, filenames):
    output = BytesIO()
    workbook = xlsxwriter.Workbook(output, {'in_memory': True})
    for file_index, file_data in enumerate(sequences_data):
        filename = filenames[file_index]
        worksheet = workbook.add_worksheet(filename[:31])
        worksheet.write(0, 0, "Entry")
        worksheet.write(0, 1, "Protein Name")
        col = 2
        for repeat in sorted(repeats):
            worksheet.write(0, col, repeat)
            col += 1
        row = 1
        for entry_id, protein_name, freq in file_data:
            worksheet.write(row, 0, entry_id)
            worksheet.write(row, 1, protein_name)
            col = 2
            for repeat in sorted(repeats):
                worksheet.write(row, col, freq.get(repeat, 0))
                col += 1
            row += 1
    workbook.close()
    output.seek(0)
    return output

# Streamlit UI
st.set_page_config(page_title="Protein Tool", layout="wide")
st.title("🧬 Protein Analysis Toolkit by SCBL, IITG")

app_choice = st.radio("Choose an option", ["🔁 Protein Repeat Finder", "📊 Protein Comparator", "🧪 Amino Acid Percentage Analyzer"])

if app_choice == "🔁 Protein Repeat Finder":
    analysis_type = st.radio("Select analysis type:", ["Homo", "Hetero", "Both"], index=2)
    uploaded_files = st.file_uploader("Upload Excel files", accept_multiple_files=True, type=["xlsx"])

    if 'all_sequences_data' not in st.session_state:
        st.session_state.all_sequences_data = []
        st.session_state.all_repeats = set()
        st.session_state.filenames = []
        st.session_state.excel_file = None

    if uploaded_files and st.button("Process Files"):
        st.session_state.all_repeats = set()
        st.session_state.all_sequences_data = []
        st.session_state.filenames = []
        for file in uploaded_files:
            excel_data = pd.ExcelFile(file)
            repeats, sequence_data = process_excel(excel_data, analysis_type)
            if repeats is not None:
                st.session_state.all_repeats.update(repeats)
                st.session_state.all_sequences_data.append(sequence_data)
                st.session_state.filenames.append(file.name)
        if st.session_state.all_sequences_data:
            st.toast(f"Processed {len(uploaded_files)} file(s) successfully.")
            st.session_state.excel_file = create_excel(
                st.session_state.all_sequences_data,
                st.session_state.all_repeats,
                st.session_state.filenames
            )

    if st.session_state.excel_file:
        st.download_button(
            label="Download Excel file",
            data=st.session_state.excel_file,
            file_name="Protein_Repeats_Analysis.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )

    # Display results table and repeat cluster visualization
    if st.checkbox("Show Results Table"):
        rows = []
        for file_index, file_data in enumerate(st.session_state.all_sequences_data):
            filename = st.session_state.filenames[file_index]
            for entry_id, protein_name, freq in file_data:
                row = {"Filename": filename, "Entry": entry_id, "Protein Name": protein_name}
                row.update({repeat: freq.get(repeat, 0) for repeat in sorted(st.session_state.all_repeats)})
                rows.append(row)
        result_df = pd.DataFrame(rows)
        st.dataframe(result_df)

        # Repeat Cluster Visualization
        repeat_counts = defaultdict(int)
        for seq_data in st.session_state.all_sequences_data:
            for _, _, freq_dict in seq_data:
                for repeat, count in freq_dict.items():
                    repeat_counts[repeat] += count

        if repeat_counts:
            sorted_repeats = sorted(repeat_counts.items(), key=lambda x: x[1], reverse=True)
            top_n = st.slider("Select number of top repeats to visualize", min_value=5, max_value=50, value=20)
            top_repeats = sorted_repeats[:top_n]
            repeats, counts = zip(*top_repeats)

            plt.figure(figsize=(12, 6))
            sns.barplot(x=list(repeats), y=list(counts), palette="viridis")
            plt.xticks(rotation=45, ha='right')
            plt.xlabel("Repeats")
            plt.ylabel("Total Frequency")
            plt.title("Top Repeat Clusters Across All Sequences")
            st.pyplot(plt.gcf())
        else:
            st.warning("No repeat data available to visualize. Please upload files first.")



elif app_choice == "📊 Protein Comparator":
    st.write("Upload two Excel files with protein data to compare repeat frequencies.")

    file1 = st.file_uploader("Upload First Excel File", type=["xlsx"], key="comp1")
    file2 = st.file_uploader("Upload Second Excel File", type=["xlsx"], key="comp2")

    if file1 and file2:
        df1 = pd.read_excel(file1)
        df2 = pd.read_excel(file2)

        df1.columns = df1.columns.astype(str)
        df2.columns = df2.columns.astype(str)

        id_col = df1.columns[0]
        name_col = df1.columns[1]
        repeat_columns = df1.columns[2:]

        diff_data = []
        for i in range(min(len(df1), len(df2))):
            row1 = df1.iloc[i]
            row2 = df2.iloc[i]
            diff_row = {"Entry": row1[id_col], "Protein Name": row1[name_col]}
            for repeat in repeat_columns:
                val1 = row1.get(repeat, 0)
                val2 = row2.get(repeat, 0)
                change = ((val2 - val1) / val1 * 100) if val1 != 0 else (100 if val2 > 0 else 0)
                diff_row[repeat] = change
            diff_data.append(diff_row)

        result_df = pd.DataFrame(diff_data)
        percent_cols = result_df.select_dtypes(include='number').columns
        st.dataframe(result_df.style.format({col: "{:.2f}%" for col in percent_cols}))

        def to_excel_with_colors(df):
            output = BytesIO()
            workbook = xlsxwriter.Workbook(output, {'in_memory': True})
            worksheet = workbook.add_worksheet('Comparison')

            green_format = workbook.add_format({'font_color': 'green'})
            red_format = workbook.add_format({'font_color': 'red'})
            header_format = workbook.add_format({'bold': True, 'bg_color': '#D7E4BC'})

            for col_num, col_name in enumerate(df.columns):
                worksheet.write(0, col_num, col_name, header_format)

            for row_num, row in enumerate(df.itertuples(index=False), start=1):
                for col_num, value in enumerate(row):
                    if col_num < 2:
                        worksheet.write(row_num, col_num, value)
                    else:
                        fmt = green_format if value > 0 else red_format if value < 0 else None
                        worksheet.write(row_num, col_num, f"{value:.2f}%", fmt)

            workbook.close()
            output.seek(0)
            return output

        excel_file = to_excel_with_colors(result_df)

        st.download_button(
            label="Download Colored Comparison Excel",
            data=excel_file,
            file_name="comparison_result_colored.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )

elif app_choice == "🧪 Amino Acid Percentage Analyzer":
    import matplotlib.pyplot as plt  # Needed for pie chart

    AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

    uploaded_file = st.file_uploader("Upload Excel file (with Entry, Protein Name, Sequence)", type=["xlsx"])

    if uploaded_file and st.button("Analyze File"):
        df = pd.read_excel(uploaded_file)

        if len(df.columns) < 3:
            st.error("The file must have at least three columns: Entry, Protein Name, Sequence")
        else:
            entry_col = df.columns[0]
            name_col = df.columns[1]
            seq_col = df.columns[2]

            from collections import Counter
            all_counts = Counter()
            all_length = 0
            result_rows = []

            for _, row in df.iterrows():
                entry = str(row[entry_col])
                name = str(row[name_col])
                sequence = str(row[seq_col]).replace(" ", "").replace("\"", "").strip().upper()
                sequence = ''.join(filter(lambda c: c in AMINO_ACIDS, sequence))
                length = len(sequence)

                if length == 0:
                    continue

                count = Counter(sequence)
                all_counts.update(count)
                all_length += length
                percentage = {aa: round(count[aa] / length * 100, 2) for aa in AMINO_ACIDS}
                result_rows.append({"Entry": entry, "Protein Name": name, **percentage})

            overall_percentage = {aa: round(all_counts[aa] / all_length * 100, 2) for aa in AMINO_ACIDS}
            overall_row = {"Entry": "OVERALL", "Protein Name": "ALL SEQUENCES", **overall_percentage}
            df_result = pd.concat([pd.DataFrame([overall_row]), pd.DataFrame(result_rows)], ignore_index=True)

            st.dataframe(df_result)

            # 🔵 Pie Chart
            st.subheader("🧁 Overall Amino Acid Composition (Pie Chart)")
            fig, ax = plt.subplots(figsize=(9, 9))
            labels = list(overall_percentage.keys())
            sizes = list(overall_percentage.values())
            filtered = [(label, size) for label, size in zip(labels, sizes) if size > 0]

            if filtered:
                labels, sizes = zip(*filtered)
                ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, counterclock=False)
                ax.axis('equal')
                st.pyplot(fig)
            else:
                st.info("No valid amino acids found to display in pie chart.")

            # Excel Export
            def to_excel(df):
                output = BytesIO()
                workbook = xlsxwriter.Workbook(output, {'in_memory': True})
                worksheet = workbook.add_worksheet("Amino Acid %")
                header_format = workbook.add_format({'bold': True, 'bg_color': '#CDEDF6'})
                for col_num, col_name in enumerate(df.columns):
                    worksheet.write(0, col_num, col_name, header_format)
                for row_num, row in enumerate(df.itertuples(index=False), start=1):
                    for col_num, value in enumerate(row):
                        worksheet.write(row_num, col_num, value)
                workbook.close()
                output.seek(0)
                return output

            excel_file = to_excel(df_result)

            st.download_button(
                label="Download Excel Report",
                data=excel_file,
                file_name="amino_acid_percentage.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
