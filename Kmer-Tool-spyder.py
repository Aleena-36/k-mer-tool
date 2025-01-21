# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 15:17:51 2025

@author: Mohib
"""

import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from collections import Counter
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd  # For data export
import numpy as np


class KmerAnalysisTool:
    def __init__(self, root):
        self.root = root
        self.root.title("K-mer Analysis Tool")
        self.root.geometry("900x700")
        self.root.configure(bg="#e8f1f2")

        # Title Label
        title_label = tk.Label(root, text="K-mer Analysis Tool", font=("Arial", 24, "bold"), bg="#007BFF", fg="white")
        title_label.pack(pady=10, fill=tk.X)

        # Input Frame
        input_frame = tk.Frame(root, bg="#e8f1f2", pady=10)
        input_frame.pack(fill=tk.X, padx=20)

        tk.Label(input_frame, text="Enter DNA Sequence:", font=("Arial", 12), bg="#e8f1f2").grid(row=0, column=0, sticky="w", pady=5)
        self.sequence_entry = tk.Text(input_frame, height=5, width=60, font=("Arial", 10), relief=tk.GROOVE, borderwidth=2)
        self.sequence_entry.grid(row=1, column=0, columnspan=2, pady=5)

        upload_button = ttk.Button(input_frame, text="Upload FASTA File", command=self.upload_fasta)
        upload_button.grid(row=2, column=0, pady=10, sticky="w")

        tk.Label(input_frame, text="K-mer Size:", font=("Arial", 12), bg="#e8f1f2").grid(row=3, column=0, sticky="w", pady=5)
        self.kmer_size_entry = ttk.Entry(input_frame, width=10)
        self.kmer_size_entry.grid(row=3, column=1, sticky="w")

        # Buttons Frame
        button_frame = tk.Frame(root, bg="#e8f1f2")
        button_frame.pack(pady=10)

        ttk.Button(button_frame, text="Perform K-mer Analysis", command=self.perform_analysis).grid(row=0, column=0, padx=10, pady=5)
        ttk.Button(button_frame, text="Detect AMR Genes", command=self.detect_amr_genes).grid(row=0, column=1, padx=10, pady=5)
        ttk.Button(button_frame, text="Generate Phylogenetic Tree", command=self.generate_phylogenetic_tree).grid(row=0, column=2, padx=10, pady=5)
        ttk.Button(button_frame, text="Generate Heatmap", command=self.generate_heatmap).grid(row=0, column=3, padx=10, pady=5)
        ttk.Button(button_frame, text="Export Data", command=self.export_data).grid(row=0, column=4, padx=10, pady=5)

        # Output Frame
        output_frame = tk.Frame(root, bg="#e8f1f2")
        output_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)

        tk.Label(output_frame, text="Output:", font=("Arial", 12, "bold"), bg="#e8f1f2").pack(anchor="w", pady=5)
        self.output_text = tk.Text(output_frame, height=10, font=("Arial", 10), relief=tk.GROOVE, borderwidth=2)
        self.output_text.pack(fill=tk.BOTH, expand=True)

    def upload_fasta(self):
        filepath = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")])
        if filepath:
            with open(filepath, "r") as file:
                lines = file.readlines()
                sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
                self.sequence_entry.delete("1.0", tk.END)
                self.sequence_entry.insert("1.0", sequence)

    def perform_analysis(self):
        sequence = self.sequence_entry.get("1.0", tk.END).strip().upper()
        try:
            k = int(self.kmer_size_entry.get())
            if k <= 0 or k > len(sequence):
                raise ValueError
        except ValueError:
            messagebox.showerror("Error", "Invalid k-mer size.")
            return

        self.kmer_counts = self.count_kmers(sequence, k)

        self.output_text.delete("1.0", tk.END)
        self.output_text.insert("1.0", f"K-mer Counts (k={k}):\n")
        for kmer, count in self.kmer_counts.items():
            self.output_text.insert(tk.END, f"{kmer}: {count}\n")

    def count_kmers(self, sequence, k):
        return Counter(sequence[i:i + k] for i in range(len(sequence) - k + 1))

    def generate_heatmap(self):
        if not hasattr(self, 'kmer_counts'):
            messagebox.showerror("Error", "No k-mer data to generate heatmap.")
            return

        kmer_list = list(self.kmer_counts.keys())
        frequencies = list(self.kmer_counts.values())
        heatmap_data = np.array([frequencies])
        sns.heatmap(heatmap_data, annot=True, fmt="d", cmap="coolwarm", xticklabels=kmer_list, yticklabels=["Frequency"])
        plt.title("K-mer Frequencies")
        plt.xlabel("K-mers")
        plt.ylabel("Frequency")
        plt.show()

    def detect_amr_genes(self):
        sequence = self.sequence_entry.get("1.0", tk.END).strip().upper()

        amr_genes = ["bla", "tet", "erm", "sul", "aac"]
        amr_kmers = ["AACCGGTT", "TTGGCCAA"]

        detected_genes = [gene for gene in amr_genes if gene.upper() in sequence]
        detected_kmers = [kmer for kmer in amr_kmers if kmer.upper() in sequence]

        self.output_text.delete("1.0", tk.END)

        if detected_genes or detected_kmers:
            self.output_text.insert(tk.END, "Detected AMR Genes and K-mers:\n")
            if detected_genes:
                self.output_text.insert(tk.END, f"Genes: {', '.join(detected_genes)}\n")
            if detected_kmers:
                self.output_text.insert(tk.END, f"K-mers: {', '.join(detected_kmers)}\n")
        else:
            self.output_text.insert(tk.END, "No AMR genes or k-mers detected.\n")

    def generate_phylogenetic_tree(self):
        try:
            alignment = MultipleSeqAlignment([
                SeqRecord(Seq("ATCGTACGAT"), id="Seq1"),
                SeqRecord(Seq("ATCGTTCGAT"), id="Seq2"),
                SeqRecord(Seq("ATCGTACCAT"), id="Seq3"),
            ])
            calculator = DistanceCalculator("identity")
            distance_matrix = calculator.get_distance(alignment)
            constructor = DistanceTreeConstructor()
            tree = constructor.upgma(distance_matrix)
            Phylo.draw(tree)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate phylogenetic tree: {str(e)}")

    def export_data(self):
        if not hasattr(self, 'kmer_counts') or not self.kmer_counts:
            messagebox.showerror("Error", "No analysis data to export.")
            return

        data = {"K-mer": list(self.kmer_counts.keys()), "Count": list(self.kmer_counts.values())}
        df = pd.DataFrame(data)

        save_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")])
        if save_path:
            df.to_csv(save_path, index=False)
            messagebox.showinfo("Success", f"Data exported successfully to {save_path}")


if __name__ == "__main__":
    root = tk.Tk()
    app = KmerAnalysisTool(root)
    root.mainloop()
