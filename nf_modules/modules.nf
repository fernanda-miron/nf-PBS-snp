#!/usr/bin/env nextflow

/*Modulos de prueba para  nextflow */
results_dir = "./test/results"
intermediates_dir = "./test/results/intermediates"

process fst_calculation {

publishDir "${results_dir}/fst_results_pop1_pop2/", mode:"copy"

	input:
	file vcf
	file pop_1
	file pop_2

	output:
	path "*.fst"

	"""
	vcftools --vcf ${vcf} \
					 --weir-fst-pop ${pop_1} \
					 --weir-fst-pop ${pop_2} \
					 --out pop1pop2
	"""
}

process fst_calculation_2 {

publishDir "${results_dir}/fst_results_pop1_popout/", mode:"copy"

	input:
	file vcf
	file pop_1
	file pop_out

	output:
	path "*.fst"

	"""
	vcftools --vcf ${vcf} \
					 --weir-fst-pop ${pop_1} \
					 --weir-fst-pop ${pop_out} \
					 --out pop1popout
	"""
}

process fst_calculation_3 {

publishDir "${results_dir}/fst_results_pop2_popout/", mode:"copy"

	input:
	file vcf
	file pop_2
	file pop_out

	output:
	path "*.fst"

	"""
	vcftools --vcf ${vcf} \
					 --weir-fst-pop ${pop_2} \
					 --weir-fst-pop ${pop_out} \
					 --out pop2popout
	"""
}

process pbs_by_snp {

	publishDir "${results_dir}/pbs_by_snp/",mode:"copy"

	input:
	file p4
	file pbs_calculator

	output:
	file "pbs*"

	"""
	mkdir dir_my_fst
	mv ${p4} dir_my_fst
	Rscript --vanilla pbs_calculator.R ./dir_my_fst "pbs_by_snp.png" "pbs.tsv"
	"""
}
