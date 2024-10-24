function probe_design_main_Github(FASTA_dir, out_dir)

[summary_table] = probe_design_targetseq_Github(FASTA_dir, out_dir);

max_probe_number = 40;
probe_design_hairpin_Github(FASTA_dir, summary_table, out_dir, max_probe_number)