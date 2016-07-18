TARGETS=`{find -L data/ -type f -name '*.sif' \
	| sed -e 's#data/#results/#g' \
		-e 's#$#.nullmodel#g' \
}

NPROC=1

intraintertest:V: $TARGETS


results/indexes/%_chrom.index:	data/%.sif
	mkdir -p `dirname $target`
	Rscript index_chromosomes_in_sif.R \
		$prereq \
		-o $target

results/%.sif_log_plot.pdf	results/%.sif_plot_pdf	results/%.sif_zoom_plot.pdf	results/%.sif.stats	results/%.sif.nullmodel: \
data/%.sif	results/indexes/%_chrom.index
	DIR=`dirname $target`
	mkdir -p $DIR
	Rscript intra_inter_comparison_script.R \
		$prereq \
		--plots \
		-o $DIR
