TARGETS=`{./targets 2>/dev/null}

intraintertest:V: $TARGETS

data/%.sif:D:	data/%.sif.bz2
	bzip2 -k -c -d `readlink -f $prereq` > $target

results/%.valid.chrom_info:	data/%.sif
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
		-o $DIR &&
	{
		test -f data/$stem.sif.bz2 && rm data/$stem.sif || true
	}

init:V:
	mkdir -p data results

clean:V:
	rm -r results

split_mi:V:	`{./targets_split_mi}

results/mi_values/%.intra \
results/mi_values/%.inter:	data/%.sif	results/indexes/%_chrom.index
	mkdir -p `dirname $target`
	./intra_inter_mi \
		$prereq \
		--output-intra 'results/mi_values/'$stem'.intra' \
		--output-inter 'results/mi_values/'$stem'.inter'
