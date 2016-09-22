TARGETS=`{./targets 2>/dev/null}

NPROC=1

intraintertest:V: $TARGETS


results/indexes/%_chrom.index:	data/%.sif
	mkdir -p `dirname $target`
	Rscript index_chromosomes_in_sif.R \
		$prereq \
		-d gene_chromosome_dictionary_july2016.txt \
		-o $target

results/%.sif_log_plot.pdf	results/%.sif_plot_pdf	results/%.sif_zoom_plot.pdf	results/%.sif.stats	results/%.sif.nullmodel: \
data/%.sif	results/indexes/%_chrom.index
	DIR=`dirname $target`
	mkdir -p $DIR
	Rscript intra_inter_comparison_script.R \
		$prereq \
		--plots \
		-o $DIR

init:V:
	mkdir -p data results

clean:V:
	rm -r results
