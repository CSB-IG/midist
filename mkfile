TARGETS=`{./targets 2>/dev/null}

NPROC=1

intraintertest:V: $TARGETS

results/%.valid.chrom_info:	data/%sif
	./chrom_info \
		$prereq \
		--vanilla results/$stem.vanilla.chrom_info \
		--output results/$stem.valid.chrom_info

results/indexes/%_chrom.index:D:	data/%.sif	results/%.valid.chrom_info
	mkdir -p `dirname $target`
	cat data/$stem.sif \
	| ./index_chromosomes_in_sif \
		results/$stem.valid.chrom_info \
	> $target

results/%.sif_log_plot.pdf	\
results/%.sif_plot_pdf	\
results/%.sif_zoom_plot.pdf	\
results/%.sif.stats	\
results/%.sif.nullmodel: \
data/%.sif	\
results/indexes/%_chrom.index
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
