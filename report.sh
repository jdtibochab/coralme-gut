find ./me-models/ -maxdepth 2 -name "*.xlsx" | cut -d / -f 3 | sort -u > step1.txt
find ./me-models/ -maxdepth 2 -name "MEModel*step2*" | cut -d / -f 3 | sort -u > step2.txt
find ./me-models/ -maxdepth 2 -name "MEModel*step3*" | cut -d / -f 3 | sort -u > step3.txt
find ./me-models/ -maxdepth 2 -name "MEBuilder*log" | xargs grep -l -e "Gene overlap between M-model and Genbank : 0%" | cut -d / -f 3 | sort -u > no_overlap.txt
find ./me-models/ -maxdepth 2 -name "curation_notes.txt" | xargs grep -l -e "Some tRNAs could be missing in the files" | cut -d / -f 3 | sort -u > trna_issues.txt

find ./me-models/ -maxdepth 2 -name "METroub*log" | xargs grep -l -e "TS_generic_" | cut -d / -f 3 | sort -u > trna_ts_issues.txt
find ./me-models/ -maxdepth 2 -name "MEModel*BIO*" | cut -d / -f 3 | sort -u > biomass_constrained.txt

find ./me-models/ -maxdepth 2 -name "curation_notes.txt" | xargs grep -l -e "Could not identify RNA polymerase" | cut -d / -f 3 | sort -u > rna_polymerase_issues.txt
find ./me-models/ -name "protein_mod.txt" -exec grep -l "_mod_zn2(1)_mod_mg2(2).*RNA_Polymerase" {} \; | cut -d / -f 3 | sort -u > zinc_in_polymerase.txt

find ./cases/fluxes/base/ -name "*.csv" -exec grep -H "sink_.heme" {} \; | cut -d / -f 5 | cut -d . -f 1 | sort -u > heme_sink.txt