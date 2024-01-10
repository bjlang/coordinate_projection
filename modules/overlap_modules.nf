process get_overlap {
	input:
	path overlap_file

	output:
	stdout

	script:
	def suffix = params.method == 'pslMap' ? (params.swapMap?'_PMs_coord_merged':'_PMns_coord_merged') : '_LOap_coord_merged'
	"""
		echo -en $overlap_file\$"\\t" | sed 's/overlap_//g' | sed 's/ATAC_//g' | sed 's/_Peaks//g' | sed 's/_to_/\\t/g' | sed 's/${suffix}//g' | sed 's/.coverage//g'; 
		bc <<< "scale=5; 1 - \$(cut -f5 $overlap_file | sort | uniq -c | grep ' 0' | xargs | cut -f1 -d' ') / \$(wc -l $overlap_file | cut -f1 -d' ')"
	"""
}

process get_base_overlap {
	input:
	path overlap_file

	output:
	stdout

	script:
	def suffix = params.method == 'pslMap' ? (params.swapMap?'_PMs_coord_merged':'_PMns_coord_merged') : '_LOap_coord_merged'
	"""
		echo -en $overlap_file\$"\\t" | sed 's/overlap_//g' | sed 's/ATAC_//g' | sed 's/_Peaks//g' | sed 's/_to_/\\t/g' | sed 's/${suffix}//g' | sed 's/.coverage//g'; 
		# bc <<< "scale=5; 1 - \$( grep -v chrUn $overlap_file | awk '{s+=\$6; t+=\$7} END {print s" / "t}' )"
		bc <<< "scale=5; \$( awk '{s+=\$6; t+=\$7} END {print s" / "t}' $overlap_file )"
	"""
}

process get_region_overlap {
	input:
	path overlap_file

	output:
	stdout

	script:
	def suffix = params.method == 'pslMap' ? (params.swapMap?'_PMs_coord_merged':'_PMns_coord_merged') : '_LOap_coord_merged'
	"""
		echo -en $overlap_file\$"\\t" | sed 's/overlap_//g' | sed 's/ATAC_//g' | sed 's/_Peaks//g' | sed 's/_to_/\\t/g' | sed 's/${suffix}//g' | sed 's/.coverage//g'; 
		bc <<< "scale=5; \$( awk '{if (\$5 > 0) print \$0}' $overlap_file | cut -f4 | sort | uniq | wc -l) / \$(cut -f4 $overlap_file | sort | uniq | wc -l)"
	"""
}

process get_genome_overlap {
	input:
	path overlap_file

	output:
	stdout

	script:
	def suffix = params.method == 'pslMap' ? (params.swapMap?'_PMs_coord_merged':'_PMns_coord_merged') : '_LOap_coord_merged'
	"""
		echo -en $overlap_file\$"\\t" | sed 's/ATAC_//g' | sed 's/_Peaks//g' | sed 's/${suffix}//g' | sed 's/.bed.coverage//g'; 
		bc <<< "scale=5; 1 - \$(grep -P 'genome\\t0' $overlap_file | cut -f5)"
	"""
}