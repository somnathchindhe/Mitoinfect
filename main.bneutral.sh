#!/bin/bash
# This script is being written by Mr. Somnath Chindhe from abopenlab, BIC, CSIR-IMTECH, on 13-10-2023, for Mitoinfect project.
# Contact: somnathchindhe@imtech.res.in or somnathchindhe.bioinfo@gmail.com
#Usage: ./main.sh ensembl_id_List.txt
# Save pal2pal.pl to binary (Uncomment these lines if needed)
# sudo mv pal2pal.pl /usr/local/bin/
# chmod +x /usr/local/bin/pal2pal.pl

mod2="B-Neutral"
cod1="F61"
cod2="F3x4"

echo -e "Ensembl_Gene_ID\tChlorocebus_sabaeus(w)\tNomascus_leucogenys(w)\tGorilla_gorilla(w)\tPan_troglodytes(w)\tPan_paniscus(w)\tHomo_sapiens(w)\tPongo_abelii(w)\tLikelihood(lnL)\tparameters(np)\tntime\tkappa(ts/tv)" > "Mitochondrail_Neutral_ratio_${cod1}_model_dnds.tsv"
echo -e "Ensembl_Gene_ID\tChlorocebus_sabaeus(w)\tNomascus_leucogenys(w)\tGorilla_gorilla(w)\tPan_troglodytes(w)\tPan_paniscus(w)\tHomo_sapiens(w)\tPongo_abelii(w)\tLikelihood(lnL)\tparameters(np)\tntime\tkappa(ts/tv)" > "Mitochondrail_Neutral_ratio_${cod2}_model_dnds.tsv"

# Assign model variables
#special script for branch neutral model


# 1. Processing the cds and protein file for alignment and ML
 #id="ENSG00000000419"  # Example ID
# The following loop iterates through a list of IDs, but it's currently commented out. You can uncomment it to run the script for each ID in the ensemble_id.txt file.
for id in $(cat $1); do

    # Run PRANK for alignment
    # prank -d=protein/"$id"__protein.fa -f=fasta -o="$id".aln +F -protein
   # clustalo -i "protein/${id}__protein.fa" -t Protein --out "${id}.pep.aln.best.fas" --output-order input-order
    # Run pal2nal.pl for alignment between protein and cds
  #  ./pal2nal.pl "${id}.pep.aln.best.fas" "cds/${id}__cds.fa" -output fasta -nomismatch -nogap > "${id}.cds.aln"
    # You have an option here, you can run pal2nal.pl as above, or use the commented line with # to use the same line you provided, which is currently commented out.
    
    # Replace the header of the alignment with species names
   # species=("Nomascus_leucogenys" "Pongo_abelii" "Gorilla_gorilla" "Homo_sapiens" "Pan_troglodytes" "Pan_paniscus" "Chlorocebus_sabaeus")
    
    #for spe in "${species[@]}"; do
      #  header=$(grep "$spe" "${id}.cds.aln")
       # new=">$spe"
      #  sed -i "s/$header/$new/g" "${id}.cds.aln"
   # done


mod2="B-Neutral"
cod1="F61"
cod2="F3x4"
codon1="0" #F61 #1
codon2="2" # F3x4 #2
omega1="1" #fix omega value to 1 i.e model assume omega one and calculated its likelyhood
omega0="0"
# Remove summary tsv files to prevent being overwritten
#################rm "ono_ratio_${cod1}_model_dnds.tsv" "ono_ratio_${cod2}_model_dnds.tsv" "free_ratio_${cod1}_model_dnds.tsv" "free_ratio_${cod2}_model_dnds.tsv" "two_ratio_${cod1}_model_dnds.tsv" "two_ratio_${cod2}_model_dnds.tsv"
cp -r new.demo.ctl "${id}_${mod2}_${cod1}.ctl"
cp -r new.demo.ctl "${id}_${mod2}_${cod2}.ctl"

# Create table files which store all dN/dS values and model parameters


sed -i "s/iiii/${id}.cds.aln/g" "${id}_${mod2}_${cod1}.ctl"
    sed -i "s/oooo/${id}_${mod2}_${cod1}/g" "${id}_${mod2}_${cod1}.ctl"
    sed -i "s/mmmm/2/g"  "${id}_${mod2}_${cod1}.ctl"
    sed -i "s/cccc/${codon1}/g"  "${id}_${mod2}_${cod1}.ctl" #########PLEASE MANUALY CHECK 
    sed -i "s/wwww/${omega1}/g"  "${id}_${mod2}_${cod1}.ctl"
    codeml "${id}_${mod2}_${cod1}.ctl"
    awk '/w ratios as labels for TreeView:/ {getline; print}' "${id}_${mod2}_${cod1}.out" > "${id}_${mod2}_${cod1}_raw_dnds.tre"


cs21=$(grep -o 'Chlorocebus_sabaeus #[0-9.]\+' "${id}_${mod2}_${cod1}_raw_dnds.tre" | cut -d '#' -f2)
nl21=$(grep -o 'Nomascus_leucogenys #[0-9.]\+' "${id}_${mod2}_${cod1}_raw_dnds.tre" | cut -d '#' -f2)
gg21=$(grep -o 'Gorilla_gorilla #[0-9.]\+' "${id}_${mod2}_${cod1}_raw_dnds.tre" | cut -d '#' -f2)
pt21=$(grep -o 'Pan_troglodytes #[0-9.]\+' "${id}_${mod2}_${cod1}_raw_dnds.tre" | cut -d '#' -f2)
pp21=$(grep -o 'Pan_paniscus #[0-9.]\+' "${id}_${mod2}_${cod1}_raw_dnds.tre" | cut -d '#' -f2)
hs21=$(grep -o 'Homo_sapiens #[0-9.]\+' "${id}_${mod2}_${cod1}_raw_dnds.tre" | cut -d '#' -f2)
pa21=$(grep -o 'Pongo_abelii #[0-9.]\+' "${id}_${mod2}_${cod1}_raw_dnds.tre" | cut -d '#' -f2)
np21=$(grep "lnL"  ${id}_${mod2}_${cod1}.out |cut -f2 -d ":" | sed 's/ //g'| sed 's/np//g')
ntime21=$(grep "lnL"  ${id}_${mod2}_${cod1}.out  |cut -f3 -d ":" | sed 's/ //g'| sed 's/)//g')
lnL21=$(grep "lnL"  ${id}_${mod2}_${cod1}.out |sed 's/..*\:\ *//' | sed 's/\ ..*//')
kappa21=$(grep "ts/tv" ${id}_${mod2}_${cod1}.out |cut -f2 -d "=" |sed 's/ //g')
echo -e "$id\t$cs21\t$nl21\t$gg21\t$pt21\t$pp21\t$hs21\t$pa21\t$lnL21\t$np21\t$ntime21\t$kappa21" >> "Mitochondrail-Neutral_ratio_${cod1}_model_dnds.tsv"

################################22. two ratio	& codon freqeuncy model="F3x4"############################################################################

sed -i "s/iiii/${id}.cds.aln/g" "${id}_${mod2}_${cod2}.ctl"
    sed -i "s/oooo/${id}_${mod2}_${cod2}/g" "${id}_${mod2}_${cod2}.ctl"
    sed -i "s/mmmm/2/g"  "${id}_${mod2}_${cod2}.ctl"
    sed -i "s/cccc/${codon2}/g"  "${id}_${mod2}_${cod2}.ctl" #########PLEASE MANUALY CHECK 
     sed -i "s/wwww/${omega1}/g"  "${id}_${mod2}_${cod2}.ctl"
    codeml "${id}_${mod2}_${cod2}.ctl"
    echo $mod2
    awk '/w ratios as labels for TreeView:/ {getline; print}' "${id}_${mod2}_${cod2}.out" > "${id}_${mod2}_${cod2}_raw_dnds.tre"


cs22=$(grep -o 'Chlorocebus_sabaeus #[0-9.]\+' "${id}_${mod2}_${cod2}_raw_dnds.tre" | cut -d '#' -f2)
nl22=$(grep -o 'Nomascus_leucogenys #[0-9.]\+' "${id}_${mod2}_${cod2}_raw_dnds.tre" | cut -d '#' -f2)
gg22=$(grep -o 'Gorilla_gorilla #[0-9.]\+' "${id}_${mod2}_${cod2}_raw_dnds.tre" | cut -d '#' -f2)
pt22=$(grep -o 'Pan_troglodytes #[0-9.]\+' "${id}_${mod2}_${cod2}_raw_dnds.tre" | cut -d '#' -f2)
pp22=$(grep -o 'Pan_paniscus #[0-9.]\+' "${id}_${mod2}_${cod2}_raw_dnds.tre" | cut -d '#' -f2)
hs22=$(grep -o 'Homo_sapiens #[0-9.]\+' "${id}_${mod2}_${cod2}_raw_dnds.tre" | cut -d '#' -f2)
pa22=$(grep -o 'Pongo_abelii #[0-9.]\+' "${id}_${mod2}_${cod2}_raw_dnds.tre" | cut -d '#' -f2)
np22=$(grep "lnL"  ${id}_${mod2}_${cod2}.out |cut -f2 -d ":" | sed 's/ //g'| sed 's/np//g')
ntime22=$(grep "lnL"  ${id}_${mod2}_${cod2}.out  |cut -f3 -d ":" | sed 's/ //g'| sed 's/)//g')
lnL22=$(grep "lnL"  ${id}_${mod2}_${cod2}.out |sed 's/..*\:\ *//' | sed 's/\ ..*//')
kappa22=$(grep "ts/tv" ${id}_${mod2}_${cod2}.out |cut -f2 -d "=" |sed 's/ //g')
echo -e "$id\t$cs22\t$nl22\t$gg22\t$pt22\t$pp22\t$hs22\t$pa22\t$lnL22\t$np22\t$ntime22\t$kappa22" >> "Mitochondrail-Neutral_ratio_${cod2}_model_dnds.tsv"


directories=( "codeml_models" "codeml_out" "codeml_out/Neutral_ratio/F61" "codeml_out/Neutral_ratio/F3X4")

# Loop to create directories only if they don't exist
for dir in "${directories[@]}"; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
 
   fi
done


# Now move respective files to destination folders
#mv "$id".pep.aln.best.fas protein_aln_prank/
#mv "$id".cds.aln cds_aln_pan2pal/
mv $id**ctl codeml_models/

cp -r "${id}_${mod2}_${cod1}.out" "codeml_out/Neutral_ratio/F61/"
cp -r "${id}_${mod2}_${cod1}_raw_dnds.tre" "codeml_out/Neutral_ratio/F61/"
cp -r "${id}_${mod2}_${cod2}.out""codeml_out/Neutral_ratio/F3X4"
cp -r "${id}_${mod2}_${cod2}_raw_dnds.tre" "codeml_out/Neutral_ratio/F3X4"

done
