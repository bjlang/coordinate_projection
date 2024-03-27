nextflow.enable.dsl=2


process	liftOver {
	container params.liftOver_container
	label 'cluster'
	storeDir "${storeDir}"
	tag { "${bed}" }

	
	input:
	tuple val(storeDir), path(bed)
	path chain

	output:
	tuple val(storeDir), path (params.suffix ? "${bed.baseName}"+params.suffix+".bed" : "${bed.baseName}_LOap_coord.bed"), emit: newCoords
	val storeDir, emit: storeDir
	path ("${bed.baseName}_ap_unmapped.bed")

	script:
    def args = task.ext.args ?: ''
	def out  = params.suffix ? "${bed.baseName}"+params.suffix+".bed" : "${bed.baseName}_LOap_coord.bed"
	"""
		liftOver $args -bedPlus=3 -tab ${bed} $chain ${out} ${bed.baseName}_ap_unmapped.bed
	"""
}

process	pslMap {
	container params.pslMap_container
	label 'cluster'
	storeDir "${storeDir}"
	tag { "${bed}" }
	label 'error_retry'

	
	input:
	tuple val(storeDir), path(psl)
	path chain

	output:
	tuple val(storeDir), path (params.suffix ? "${psl.baseName}"+params.suffix+".psl" : (params.swapMap ? "${psl.baseName}_PMs_coord.psl" : "${psl.baseName}_PMns_coord.psl")), emit: newCoords
	val storeDir, emit: storeDir

	script:
	def args = params.swapMap ? '-chainMapFile -swapMap' : '-chainMapFile'
	def out  = params.suffix ? "${psl.baseName}"+params.suffix+".psl" : (params.swapMap ? "${psl.baseName}_PMs_coord.psl" : "${psl.baseName}_PMns_coord.psl")
	"""
		pslMap $args -mapInfo=${psl.baseName}_mapInfo.file $psl $chain ${out}
	"""
}

process bedToPsl {
	container params.bedtopsl_container

	input:
	tuple val(meta), path(bed)
	path genome

	output:
	tuple val(meta), path("${bed.baseName}.psl")

	script:
	"""
		cut -f1,2,3,4 $bed > bed4file.bed
		bedToPsl $genome bed4file.bed ${bed.baseName}.psl
	"""
}

process pslToBed {
	container params.psltobed_container

	input:
	tuple val(meta), path(psl)

	output:
	tuple val(meta), path("${psl.baseName}.bed")

	script:
	"""
		psl2bed --split < $psl > ${psl.baseName}.bed
	"""
}

process mergeBed {
	container params.bedTools_container
	storeDir "${storeDir}"
	label 'error_retry'

	input:
	tuple val(storeDir), path(bed)

	output:
	path ("${bed.baseName}_merged.bed")

	script:
	"""
		bedtools merge -c 4 -o collapse -i <(sort -k1,1 -k2,2n ${bed}) > ${bed.baseName}_merged.bed
	"""
}

process intersectBed {
	container params.bedTools_container
	storeDir "${storeDir}"
	//publishDir "${storeDir}", mode:"copy"
	label 'error_retry'

	input:
	tuple val(storeDir), path(bed_a), val(bed_b), val(bed_out)

	output:
	path bed_out
	
	script:
	"""
		bedtools intersect -a ${bed_a} -b ${bed_b} > ${bed_out}
	"""
}

process multiIntersectBed {
	container params.bedTools_container
	storeDir "${storeDir}"
	//publishDir "${storeDir}", mode:"copy"
	label 'error_retry'

	input:
	tuple val(storeDir), val(bed_list), val(num), val(bed_out)

	output:
	path bed_out
	
	script:
	"""
		if [[ $num -eq 1 ]]
		then
			cp ${bed_list} ${bed_out}
		else
			#bedtools multiinter -i ${bed_list} | awk '{if (\$4==${num}) print \$0}' > ${bed_out}
			bedtools multiinter -i ${bed_list} | awk '{if (\$4==${num}) print \$1,\$2,\$3,\$1"_"\$2"_"\$3}' OFS='\t' > ${bed_out}
		fi
	"""
}

process genomeCoverageBed {
	container params.bedTools_container
	tag { "${bed}" }
	
	input:
	path bed
	path ref

	output:
	path ("${bed}.coverage")

	script:
	"""
		sort -k1,1 -k2,2n ${bed} > ${bed.baseName}.sorted.bed
		#sort -k1,1 -k2,2n ${bed} | grep -v chrUn > ${bed.baseName}.sorted.bed
		bedtools genomecov -i ${bed.baseName}.sorted.bed -g ${ref} > ${bed}.coverage
	"""
}

process coverageBed{
	container params.bedTools_container
	storeDir "${params.cache}"
	//publishDir "${params.cache}", mode:"copy"
	label 'error_retry'

	input:
	each path(ref_bed)
	tuple( val(query_name), file("query.bed") )


	output:
	path ("overlap_${ref_bed.baseName}_to_${query_name}.coverage")

	script:
	"""	
		# from v2.24 on a and b are swapped
		bedtools coverage -b <(sort -k1,1 -k2,2n ${ref_bed} | cut -f1,2,3,4) -a <(sort -k1,1 -k2,2n query.bed | cut -f1,2,3,4) > overlap_${ref_bed.baseName}_to_${query_name}.coverage
	"""
}

process coverageBed_rev{
	container params.bedTools_container
	storeDir "${params.rev_cache}"
	//publishDir "${params.cache}", mode:"copy"
	label 'error_retry'

	input:
	each path(ref_bed)
	tuple( val(query_name), file("query.bed") )


	output:
	path ("overlap_${ref_bed.baseName}_to_${query_name}.coverage")

	script:
	"""	
		# from v2.24 on a and b are swapped
		bedtools coverage -a <(sort -k1,1 -k2,2n ${ref_bed} | cut -f1,2,3,4) -b <(sort -k1,1 -k2,2n query.bed | cut -f1,2,3,4) > overlap_${ref_bed.baseName}_to_${query_name}.coverage
	"""
}

process jaccard{
	container params.bedTools_container

	input:
	tuple val(ref_name) , path(ref_bed) , val(query_name), file("query.bed") 


	output:
	stdout

	script:
	"""	
		sort -k1,1 -k2,2n $ref_bed | cut -f1,2,3,4 > sorted_$ref_bed
		sort -k1,1 -k2,2n query.bed | cut -f1,2,3,4 > sorted_query.bed
		echo -en $ref_name\$"\\t"$query_name\$"\\t";
		bedtools jaccard -b sorted_${ref_bed} -a sorted_query.bed | cut -f3 | tail -n1
	"""
}

include { get_overlap 		 as get_overlap_sens 		} from './modules/overlap_modules'
include { get_overlap 		 as get_overlap_prec 		} from './modules/overlap_modules'
include { get_base_overlap	 as get_base_overlap_sens	} from './modules/overlap_modules'
include { get_base_overlap	 as get_base_overlap_prec	} from './modules/overlap_modules'
include { get_region_overlap as get_region_overlap_sens	} from './modules/overlap_modules'
include { get_region_overlap as get_region_overlap_prec	} from './modules/overlap_modules'
include { get_genome_overlap 							} from './modules/overlap_modules'



refTissues = Channel.fromPath(params.referenceTissues + "*/", type: 'dir', relative: true)
		    .map {t -> tuple("$t", file("$params.referenceTissues/${t}/"))}
qryTissues = Channel.fromPath(params.queryTissues + "*/", type: 'dir', relative: true) 
		    .map {t -> tuple("$t", file("$params.queryTissues/${t}/"))}
		    

def filterFiles( f_list ) {
	return f_list.findAll { !it.toString().contains("_coord.bed") && !it.toString().contains("unmapped.bed") && !it.toString().contains("merged.bed") && !it.toString().contains("cons.bed") && !it.toString().contains("Intrsct.bed") }
}


def matrixOutput( data_jaccard, data_list_prec , data_list_sens , genome_data_input) {
	println "Precision"
	genomeData = genome_data_input.toSpreadMap()
	prev = data_list_prec[0][0]
	data_list_prec.each { if(it[0] == prev) {System.out.print "\t${it[1]}"} else {} }
	System.out.print "\tgenome\n$prev\t"
	data_list_prec.each { if(it[0] == prev) {System.out.print "${it[2]}\t"} else {System.out.print "${genomeData.get(prev)}\n${it[0]}\t${it[2]}\t"; prev = it[0]} }
	System.out.print "${genomeData.get(prev)}\ngenome"
	data_list_prec.each { if(it[0] == prev) {System.out.print "\t${genomeData.get(it[1]) ?: genomeData.get(it[1].split('_')[1])}"} else {} }
	println ""
	println "Sensitivity"
	prev = data_list_sens[0][0]
	data_list_sens.each { if(it[0] == prev) {System.out.print "\t${it[1]}"} else {} }
	System.out.print "\tgenome\n$prev\t"
	data_list_sens.each { if(it[0] == prev) {System.out.print "${it[2]}\t"} else {System.out.print "${genomeData.get(prev)}\n${it[0]}\t${it[2]}\t"; prev = it[0]} }
	System.out.print "${genomeData.get(prev)}\ngenome"
	data_list_sens.each { if(it[0] == prev) {System.out.print "\t${genomeData.get(it[1]) ?: genomeData.get(it[1].split('_')[1])}"} else {} }
	println ""
	println "Jaccard"
	prev = data_jaccard[0][0]
	data_jaccard.each { if(it[0] == prev) {System.out.print "\t${it[1]}"} else {} }
	System.out.print "\tgenome\n$prev\t"
	data_jaccard.each { if(it[0] == prev) {System.out.print "${it[2]}\t"} else {System.out.print "${genomeData.get(prev)}\n${it[0]}\t${it[2]}\t"; prev = it[0]} }
	System.out.print "${genomeData.get(prev)}\ngenome"
	data_jaccard.each { if(it[0] == prev) {System.out.print "\t${genomeData.get(it[1]) ?: genomeData.get(it[1].split('_')[1])}"} else {} }
	println ""
}


workflow {
	tissues = refTissues.join(qryTissues)

	ch_refRegionFiles = tissues.flatMap { t -> filterFiles( file(t[1] + "/*.bed") ) }.unique()
	ch_qryRegionFiles = tissues.flatMap { t -> filterFiles( file(t[2] + "/*.bed") ) }


	// create channel of shape [<parent folder>, <space separated file list>, <list length>, <file name>] as input for multiIntersectBed
	tissues	.map{ t->tuple(	t[1],
						filterFiles(file(t[1] + "/*.bed")).inject(""){ a,b -> a+" "+b }, 
						filterFiles(file(t[1] + "/*.bed")).size(),
						t[0]+"_refIntrsct.bed" ) }
			.mix( tissues
				  .map{ t->tuple( t[2],
						filterFiles(file(t[2] + "/*.bed")).inject(""){ a,b -> a+" "+b }, 
						filterFiles(file(t[2] + "/*.bed")).size(),
						t[0]+"_qryIntrsct.bed" ) } )
			.filter{it[1] != null}
			.set {combined_tissues_ch}
	multiIntersectBed(combined_tissues_ch) 
			| branch {
	  			ref:	it.parent.parent == file(params.referenceTissues)
	  			qry:	it.parent.parent == file(params.queryTissues)
	  			other:	true  }
			| set {ch_intersectedFiles} 

	if ( params.full_intersection ) {
		ch_intersected_qryRegionFiles = tissues.flatMap { t -> filterFiles( file(t[2] + "/*.bed") ) }
											   .mix(ch_intersectedFiles.qry)
	} else {
		ch_intersected_qryRegionFiles = ch_intersectedFiles.qry
	}

	if (params.intersect_conservation) {
		ch_intersected_qryRegionFiles
			   .map { t -> tuple(t.parent, t, params.queryConservation,  t.baseName+"_cons.bed") }
			   | intersectBed
		ch_qryRegionFiles = intersectBed.out 
	} else {
		ch_qryRegionFiles = ch_intersected_qryRegionFiles
	}


    if (params.method == 'pslMap') {
		qryRegion_Psl  = bedToPsl( ch_qryRegionFiles.map{ t -> tuple(t.parent, t)}, params.queryReferenceGenome ) 
		pslMap( qryRegion_Psl , params.swapMap ? params.swapped_chain : params.chain ) 
		lifted_qry_bed = pslToBed( pslMap.out.newCoords )
	} else {
		liftOver( ch_qryRegionFiles.map{ t -> tuple(t.parent, t)} , params.chain )
		lifted_qry_bed = liftOver.out.newCoords
	}
	ch_lifted_qry_merged = mergeBed(lifted_qry_bed )

	if ( params.full_intersection ) {
		cov_file 	 = coverageBed( ch_refRegionFiles.mix(ch_intersectedFiles.ref) , 
		 							ch_lifted_qry_merged.map{ t -> tuple(t.baseName.contains("_qryIntrsct") ? t.baseName : t.parent.baseName+"_"+t.baseName, t) } )
		cov_file_rev = coverageBed_rev( ch_refRegionFiles.mix(ch_intersectedFiles.ref) , 
				 						ch_lifted_qry_merged.map{ t -> tuple(t.baseName.contains("_qryIntrsct") ? t.baseName : t.parent.baseName+"_"+t.baseName, t) } )
		jaccard( ch_refRegionFiles.mix(ch_intersectedFiles.ref).map{ t -> tuple(t.baseName.contains("_refIntrsct") ? t.baseName : t.parent.baseName+"_"+t.baseName, t) }
					.combine( ch_lifted_qry_merged.map{ t -> tuple(t.baseName.contains("_qryIntrsct") ? t.baseName : t.parent.baseName+"_"+t.baseName, t) } ) )
	} else {
		cov_file 	 = coverageBed( ch_intersectedFiles.ref ,
				 					ch_lifted_qry_merged.map{ t -> tuple(t.baseName.contains("_qryIntrsct") ? t.baseName : t.parent.baseName+"_"+t.baseName, t) } )
		cov_file_rev = coverageBed_rev( ch_intersectedFiles.ref , 
				 						ch_lifted_qry_merged.map{ t -> tuple(t.baseName.contains("_qryIntrsct") ? t.baseName : t.parent.baseName+"_"+t.baseName, t) } )
		jaccard( ch_intersectedFiles.ref.map{ t -> tuple(t.baseName.contains("_refIntrsct") ? t.baseName : t.parent.baseName+"_"+t.baseName, t) }
					.combine( ch_lifted_qry_merged.map{ t -> tuple(t.baseName.contains("_qryIntrsct") ? t.baseName : t.parent.baseName+"_"+t.baseName, t) } ) )
	}

			 
	peak_overlap = Channel.empty()
    if (params.compute_overlap == 'base') {
		peak_overlap 	 = get_base_overlap_prec(cov_file)
		peak_overlap_rev = get_base_overlap_sens(cov_file_rev)
	} else if (params.compute_overlap == 'region') {
		peak_overlap 	 = get_region_overlap_prec(cov_file)
		peak_overlap_rev = get_region_overlap_sens(cov_file_rev)
	} else {
		peak_overlap 	 = get_overlap_prec(cov_file)
		peak_overlap_rev = get_overlap_sens(cov_file_rev)
	}

	if ( params.full_intersection ) {
		genomeCoverageBed( ch_lifted_qry_merged.mix( ch_refRegionFiles ).mix(ch_intersectedFiles.ref) , params.targetReferenceGenome )
	} else {
		genomeCoverageBed( ch_lifted_qry_merged.mix(ch_intersectedFiles.ref) , params.targetReferenceGenome )	
	}
	genome_overlap = get_genome_overlap(genomeCoverageBed.out)
		.splitCsv(sep:'\t')
		//.groupTuple(size:1)
		.flatten()
		.toList()
		.map { it -> [0, it]}
		.set{ch_genome_data}

	peak_overlap
		.splitCsv(sep:'\t')
		.toSortedList({ x,y -> x[0] <=> y[0] ?: x[1] <=> y[1] })
		.map { it -> [0, it]}
		.set{ch_prec_data}

	peak_overlap_rev
		.splitCsv(sep:'\t')
		.toSortedList({ x,y -> x[0] <=> y[0] ?: x[1] <=> y[1] })
		.map { it -> [0, it]}
		.set{ch_sens_data} 

	jaccard.out
		.splitCsv(sep:'\t')
		.toSortedList({ x,y -> x[0] <=> y[0] ?: x[1] <=> y[1] })
		.map{ it -> [0, it]}
		.set{ch_jaccard_data} 

	ch_jaccard_data.join(ch_prec_data, by:0)
				.join(ch_sens_data, by:0)
				.join(ch_genome_data, by:0)
				.map{it -> matrixOutput(it[1], it[2], it[3], it[4])}
}




