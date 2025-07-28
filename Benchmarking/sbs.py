from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerAssignment import Analyzer as Analyze

genInstall.install('GRCh38')
genInstall.install('GRCh37')


Analyze.cosmic_fit("/Users/liqun/scMutraceProject/Data/EGA_AML/SBS/SBS/", "/Users/liqun/scMutraceProject/Data/EGA_AML/SBS/SBS/", input_type="vcf", context_type="96",
                    collapse_to_SBS96=True, cosmic_version=3.4, exome=False,
                    genome_build="GRCh37", signature_database=None,
                    exclude_signature_subgroups=None, export_probabilities=True,
                    export_probabilities_per_mutation=False, make_plots=True,
                    sample_reconstruction_plots=True, verbose=True)

Analyze.cosmic_fit("/Users/liqun/scMutraceProject/Data/EGA_AML/SBS/SBSSkin/", "/Users/liqun/scMutraceProject/Data/EGA_AML/SBS/SBSSkin/", input_type="vcf", context_type="96",
                    collapse_to_SBS96=True, cosmic_version=3.4, exome=False,
                    genome_build="GRCh37", signature_database=None,
                    exclude_signature_subgroups=None, export_probabilities=True,
                    export_probabilities_per_mutation=False, make_plots=True,
                    sample_reconstruction_plots=True, verbose=False)


Analyze.cosmic_fit("/Users/liqun/scMutraceProject/Data/EGA_AML/SBS/SBS/Mutational_Matrix_SBS.txt", "/Users/liqun/scMutraceProject/Data/EGA_AML/SBS/SBS/", input_type="matrix", context_type="96",
                    collapse_to_SBS96=True, cosmic_version=3.4, exome=False,
                    genome_build="GRCh38", signature_database=None,
                    exclude_signature_subgroups=None, export_probabilities=False,
                    export_probabilities_per_mutation=False, make_plots=False,
                    sample_reconstruction_plots=True, verbose=False)