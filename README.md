# GDAN AWG SV Consensus
Adapted from: https://bitbucket.org/weischenfeldt/pcawg_sv_merge/src/docker/

## Inputs:
BEDPE format **must** have the following columns. The code will error very early on without this exact format:
```
chrom1	start1	end1	chrom2	start2	end2	sv_id	tumreads	strand1	strand2	svclass	svmethod
```

- Please keep the column header "chrom1" labeled as such -- this substring is specifically searched for later in the code.
- Do NOT include 'chr' in your chromosome columns. We've tried adding error handling to bypass this, but best use would be to not include this prefix.

## Dependencies
bedtools (for pairToPair) >= v.2.28.0
Python dependencies (using python3): bgzip, pandas, pyvcf, matplotlib, networkx, numpy, os, re, sys

## Running the code
Sample test command: 
```
scripts/optimizing_sv_consensus.bash -a HCM-BROD-0002-C71-85A -o ./test_output_HCM-BROAD-0002 -i "./test_bedpe/HCM-BROD-0002-C71-85A_mskcc.bedpe ./test_bedpe/HCM-BROD-0002-C71-85A_nygc.bedpe ./test_bedpe/HCM-BROD-0002-C71-85A_purple.bedpe" -n "mskcc nygc washu"
```

## Outputs:
Consensus SV BEDPE (aliquot_id.somatic.sv.bedpe) with following headers:
```
chrom1	start1	end1	chrom2	start2	end2	sv_id	pe_support	strand1	strand2	sv_class	svmethod	center
```


-----------

# Notes for us, will need to clean this up later:

## To-Do:
- Some additional error checking with generate_bedpe_new (or perhaps the clique generation beforehand?). The # consensus SVs seem very different from numbers SURVIVOR is generating for these samples, though it is unclear exactly why yet. Reached out to Seongmin to check on parameters
- Should we make filtering steps optional (add a flag for this)? I think we will need to perform liftover of blacklist_files for this

## Error checking I have tried to introduce:
- Input bedpe files should have ^I (tabs) between different fields and no ^M (carriage return) at end of line, otherwise code will fail (error detected 11/14 when testing CTSP-ACY8) (11/14)
- Notes emphasizing need for specific header formatting since this caused numerous issues and groups have not been following this (11/15)
- Extra 'chr' removal steps to handle some of the above (11/15)

## Directories / files we can probably remove(?):
- ./data subdirectories: 161127.sv_merge, brass, delly, dranger, snowman
- ./data files: pcawg_release_mar2016.tsv, PCAWG-QC_Summary-of-Measures.tsv
- ./orig_test_bedpe_copy
- ./plots
- ./scripts files: pcawg_multiple_sv_merging_heatmap.R, pcawg6_make_readdepth_plot.sh, pcawg6_sv_merge_annotate.py, remove_splicing_type_svs_yl.sh, run_sv_merge.py
- ./test_bedpe
- ./test_output_*
- ./temp_file
- ./tmpfile
- ./unresolved.txt
