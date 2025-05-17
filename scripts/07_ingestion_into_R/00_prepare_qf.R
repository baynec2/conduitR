################################################################################
# Preparing QFeatures object from matrcies
################################################################################
# Inputs
annotation = snakemake@input[["annotation"]]
precursor_matrix = snakemake@input[["precursor_matrix"]]
peptide_matrix = snakemake@input[["peptide_matrix"]]
protein_group_matrix = snakemake@input[["protein_group_matrix"]]
superkingdom_matrix = snakemake@input[["superkingdom_matrix"]]
kingdom_matrix= snakemake@input[["kingdom_matrix"]]
phylum_matrix= snakemake@input[["phylum_matrix"]]
class_matrix= snakemake@input[["class_matrix"]]
order_matrix= snakemake@input[["order_matrix"]]
family_matrix= snakemake@input[["family_matrix"]]
genus_matrix= snakemake@input[["genus_matrix"]]
species_matrix= snakemake@input[["species_matrix"]]
go_matrix= snakemake@input[["go_matrix"]]
go_taxa_matrix = snakemake@input[["go_taxa_matrix"]]
subcellular_locations_matrix = snakemake@input[["subcellular_locations_matrix"]]

# Output
qf = snakemake@output[["qf"]]


# Testing
#annotation = "user_input/sample_annotation.txt"
#precursor_matrix="output/05_output_files/precursor_matrix.tsv"
#peptide_matrix = "output/05_output_files/peptide_matrix.tsv"
#protein_group_matrix="output/05_output_files/protein_group_matrix.tsv"
#superkingdom_matrix="output/05_output_files/superkingdom_matrix.tsv"
#kingdom_matrix="output/05_output_files/kingdom_matrix.tsv"
#phylum_matrix="output/05_output_files/phylum_matrix.tsv"
#class_matrix="output/05_output_files/class_matrix.tsv"
#order_matrix="output/05_output_files/order_matrix.tsv"
#family_matrix="output/05_output_files/family_matrix.tsv"
#genus_matrix="output/05_output_files/genus_matrix.tsv"
#species_matrix="output/05_output_files/species_matrix.tsv"
#go_matrix = "output/05_output_files/go_matrix.tsv"
#go_taxa_matrix = "output/05_output_files/go_taxa_matrix.tsv"
#subcellular_locations_matrix = "output/05_output_files/subcellular_locations_matrix.tsv"

#qf = "output/05_output_files/qf.rds"
#Defining vector of filepaths
vector_of_matrix_fps = c(precursor_matrix,peptide_matrix,protein_group_matrix,
                         superkingdom_matrix,kingdom_matrix,phylum_matrix,
                         class_matrix,order_matrix,family_matrix,genus_matrix,
                         species_matrix,go_matrix,go_taxa_matrix,subcellular_locations_matrix)



# Loading QFeatures object
QF = conduitR::prepare_qfeature(sample_annotation_fp = annotation,
                                         vector_of_matrix_fps)


saveRDS(QF,qf)
