#!/bin/bash
# This script is being written by Mr. Somnath Chindhe from abopenlab, BIC, CSIR-IMTECH, on 13-10-2023, for Mitoinfect project.
# Contact: somnathchindhe@imtech.res.in or somnathchindhe.bioinfo@gmail.com
#Usage: ./main.sh ensembl_id_List.txt
# Save pal2pal.pl to binary (Uncomment these lines if needed)
# sudo mv pal2pal.pl /usr/local/bin/
# chmod +x /usr/local/bin/pal2pal.pl

# Assign model variables
mod0="0"
mod1="1"
mod2="2"
cod1="F61"
cod2="F3x4"

# Remove summary tsv files to prevent being overwritten
#################rm "ono_ratio_${cod1}_model_dnds.tsv" "ono_ratio_${cod2}_model_dnds.tsv" "free_ratio_${cod1}_model_dnds.tsv" "free_ratio_${cod2}_model_dnds.tsv" "two_ratio_${cod1}_model_dnds.tsv" "two_ratio_${cod2}_model_dnds.tsv"

# Create table files which store all dN/dS values and model parameters
echo -e "Ensembl_Gene_ID\tChlorocebus_sabaeus(w)\tNomascus_leucogenys(w)\tGorilla_gorilla(w)\tPan_troglodytes(w)\tPan_paniscus(w)\tHomo_sapiens(w)\tPongo_abelii(w)\tLikelihood(lnL)\tparameters(np)\tntime\tkappa(ts/tv)" > "free_ratio_${cod1}_model_dnds.tsv"
echo -e "Ensembl_Gene_ID\tChlorocebus_sabaeus(w)\tNomascus_leucogenys(w)\tGorilla_gorilla(w)\tPan_troglodytes(w)\tPan_paniscus(w)\tHomo_sapiens(w)\tPongo_abelii(w)\tLikelihood(lnL)\tparameters(np)\tntime\tkappa(ts/tv)" > "free_ratio_${cod2}_model_dnds.tsv"
echo -e "Ensembl_Gene_ID\tChlorocebus_sabaeus(w)\tNomascus_leucogenys(w)\tGorilla_gorilla(w)\tPan_troglodytes(w)\tPan_paniscus(w)\tHomo_sapiens(w)\tPongo_abelii(w)\tLikelihood(lnL)\tparameters(np)\tntime\tkappa(ts/tv)" > "two_ratio_${cod1}_model_dnds.tsv"
echo -e "Ensembl_Gene_ID\tChlorocebus_sabaeus(w)\tNomascus_leucogenys(w)\tGorilla_gorilla(w)\tPan_troglodytes(w)\tPan_paniscus(w)\tHomo_sapiens(w)\tPongo_abelii(w)\tLikelihood(lnL)\tparameters(np)\tntime\tkappa(ts/tv)" > "two_ratio_${cod2}_model_dnds.tsv"
echo -e "Ensembl_Gene_ID\tOmega(w)\tLikelihood(lnL)\tparameters(np)\tntime\tkappa(ts/tv)" > "ono_ratio_${cod1}_model_dnds.tsv"
echo -e "Ensembl_Gene_ID\tOmega(w)\tLikelihood(lnL)\tparameters(np)\tntime\tkappa(ts/tv)" > "ono_ratio_${cod2}_model_dnds.tsv"

# Generate ensemble_id.txt containing all the ensemble gene IDs which are further processed for alignment and dn/ds calculation
# ls cds/*cds.fa | cut -f2 -d "/" | cut -f1 -d "_" > ensemble_id.txt 

# 1. Processing the cds and protein file for alignment and ML
 #id="ENSG00000000419"  # Example ID
# The following loop iterates through a list of IDs, but it's currently commented out. You can uncomment it to run the script for each ID in the ensemble_id.txt file.
for id in $(cat $1); do
    # Run PRANK for alignment
    # prank -d=protein/"$id"__protein.fa -f=fasta -o="$id".aln +F -protein
    clustalo -i "protein/${id}__protein.fa" -t Protein --out "${id}.pep.aln.best.fas" --output-order input-order
    # Run pal2nal.pl for alignment between protein and cds
    ./pal2nal.pl "${id}.pep.aln.best.fas" "cds/${id}__cds.fa" -output fasta -nomismatch -nogap > "${id}.cds.aln"
    # You have an option here, you can run pal2nal.pl as above, or use the commented line with # to use the same line you provided, which is currently commented out.
    
    # Replace the header of the alignment with species names
    species=("Nomascus_leucogenys" "Pongo_abelii" "Gorilla_gorilla" "Homo_sapiens" "Pan_troglodytes" "Pan_paniscus" "Chlorocebus_sabaeus")
    
    for spe in "${species[@]}"; do
        header=$(grep "$spe" "${id}.cds.aln")
        new=">$spe"
        sed -i "s/$header/$new/g" "${id}.cds.aln"
    done
    
    
    mod0="0" #M0 model
mod1="1" #free ratio model
mod2="2" #Two ratio model
cod1="F61" #  0 IN CTL FILE
cod2="F3x4" # 2 in CTL FILE
codon1="0" #F61 #1
codon2="2" # F3x4 #2

##############
cp -r demo.ctl "${id}_${mod0}_${cod1}.ctl"
cp -r demo.ctl "${id}_${mod1}_${cod1}.ctl"
cp -r demo.ctl "${id}_${mod2}_${cod1}.ctl"
cp -r demo.ctl "${id}_${mod0}_${cod2}.ctl"
cp -r demo.ctl "${id}_${mod1}_${cod2}.ctl"
cp -r demo.ctl "${id}_${mod2}_${cod2}.ctl"

################################################***11.free model , cod1="FX61"###########################################################################################################
#control file processing for free model , cod1="FX61" an run codeml model

sed -i "s/iiii/${id}.cds.aln/g" "${id}_${mod1}_${cod1}.ctl"
sed -i "s/oooo/${id}_${mod1}_${cod1}/g" "${id}_${mod1}_${cod1}.ctl"
sed -i "s/mmmm/${mod1}/g"  "${id}_${mod1}_${cod1}.ctl"
sed -i "s/cccc/${codon1}/g"  "${id}_${mod1}_${cod1}.ctl"
codeml "${id}_${mod1}_${cod1}.ctl"

awk '/w ratios as labels for TreeView:/ {getline; print}' "${id}_${mod1}_${cod1}.out" > "${id}_${mod1}_${cod1}_raw_dnds.tre"

#to perse result into table 
#ENSG00000000419.out.out:(Chlorocebus_sabaeus #0.0789804 , Nomascus_leucogenys #0.0001 , ((Gorilla_gorilla #6.08045 , ((Pan_troglodytes #8.58569 , Pan_paniscus #28.4795 ) #8.39211 , Homo_sapiens #8.74115 ) #6.68481 ) #0.0001 , Pongo_abelii #1.92827 ) #0.0001 )
cs11=$(grep -o 'Chlorocebus_sabaeus #[0-9.]\+' "${id}_${mod1}_${cod1}_raw_dnds.tre" | cut -d '#' -f2)
nl11=$(grep -o 'Nomascus_leucogenys #[0-9.]\+' "${id}_${mod1}_${cod1}_raw_dnds.tre" | cut -d '#' -f2)
gg11=$(grep -o 'Gorilla_gorilla #[0-9.]\+' "${id}_${mod1}_${cod1}_raw_dnds.tre" | cut -d '#' -f2)
pt11=$(grep -o 'Pan_troglodytes #[0-9.]\+' "${id}_${mod1}_${cod1}_raw_dnds.tre" | cut -d '#' -f2)
pp11=$(grep -o 'Pan_paniscus #[0-9.]\+' "${id}_${mod1}_${cod1}_raw_dnds.tre" | cut -d '#' -f2)
hs11=$(grep -o 'Homo_sapiens #[0-9.]\+' "${id}_${mod1}_${cod1}_raw_dnds.tre" | cut -d '#' -f2)
pa11=$(grep -o 'Pongo_abelii #[0-9.]\+' "${id}_${mod1}_${cod1}_raw_dnds.tre" | cut -d '#' -f2)
np11=$(grep "lnL"  ${id}_${mod1}_${cod1}.out |cut -f2 -d ":" | sed 's/ //g'| sed 's/np//g')
ntime11=$(grep "lnL"  ${id}_${mod1}_${cod1}.out  |cut -f3 -d ":" | sed 's/ //g'| sed 's/)//g')
lnL11=$(grep "lnL"  ${id}_${mod1}_${cod1}.out |sed 's/..*\:\ *//' | sed 's/\ ..*//')
kappa11=$(grep "ts/tv" ${id}_${mod1}_${cod1}.out |cut -f2 -d "=" |sed 's/ //g')
#kappa11=$(grep "ts/tv" ${id}_${mod1}_${cod1}.out |cut -f2 -d "=" |sed 's/ //g')

echo -e "$id\t$cs11\t$nl11\t$gg11\t$pt11\t$pp11\t$hs11\t$pa11\t$lnL11\t$np11\t$ntime11\t$kappa11" >> "free_ratio_${cod1}_model_dnds.tsv" 

#create table file containing dnds a for all branches and lnl, np n ntime value

################################################**1.free model  && 2.codon freqeuncy model="F3x4"##############################################################################
sed -i "s/iiii/${id}.cds.aln/g" "${id}_${mod1}_${cod2}.ctl"
    sed -i "s/oooo/${id}_${mod1}_${cod2}/g" "${id}_${mod1}_${cod2}.ctl"
    sed -i "s/mmmm/${mod1}/g"  "${id}_${mod1}_${cod2}.ctl"
    sed -i "s/cccc/${codon2}/g"  "${id}_${mod1}_${cod2}.ctl" #########PLEASE MANUALY CHECK
    codeml "${id}_${mod1}_${cod2}.ctl"
    awk '/w ratios as labels for TreeView:/ {getline; print}' "${id}_${mod1}_${cod2}.out" > "${id}_${mod1}_${cod2}_raw_dnds.tre"

    # To parse results into a table
    # ENSG00000000419.out.out:(Chlorocebus_sabaeus #0.0789804 , Nomascus_leucogenys #0.0001 , ((Gorilla_gorilla #6.08045 , ((Pan_troglodytes #8.58569 , Pan_paniscus #28.4795 ) #8.39211 , Homo_sapiens #8.74115 ) #6.68481 ) #0.0001 , Pongo_abelii #1.92827 ) #0.0001 )
cs12=$(grep -o 'Chlorocebus_sabaeus #[0-9.]\+' "${id}_${mod1}_${cod2}_raw_dnds.tre" | cut -d '#' -f2)
nl12=$(grep -o 'Nomascus_leucogenys #[0-9.]\+' "${id}_${mod1}_${cod2}_raw_dnds.tre" | cut -d '#' -f2)
gg12=$(grep -o 'Gorilla_gorilla #[0-9.]\+' "${id}_${mod1}_${cod2}_raw_dnds.tre" | cut -d '#' -f2)
pt12=$(grep -o 'Pan_troglodytes #[0-9.]\+' "${id}_${mod1}_${cod2}_raw_dnds.tre" | cut -d '#' -f2)
pp12=$(grep -o 'Pan_paniscus #[0-9.]\+' "${id}_${mod1}_${cod2}_raw_dnds.tre" | cut -d '#' -f2)
hs12=$(grep -o 'Homo_sapiens #[0-9.]\+' "${id}_${mod1}_${cod2}_raw_dnds.tre" | cut -d '#' -f2)
pa12=$(grep -o 'Pongo_abelii #[0-9.]\+' "${id}_${mod1}_${cod2}_raw_dnds.tre" | cut -d '#' -f2)
np12=$(grep "lnL"  ${id}_${mod1}_${cod2}.out |cut -f2 -d ":" | sed 's/ //g'| sed 's/np//g')
ntime12=$(grep "lnL"  ${id}_${mod1}_${cod2}.out  |cut -f3 -d ":" | sed 's/ //g'| sed 's/)//g')
lnL12=$(grep "lnL"  ${id}_${mod1}_${cod2}.out |sed 's/..*\:\ *//' | sed 's/\ ..*//')
kappa12=$(grep "ts/tv" ${id}_${mod1}_${cod2}.out |cut -f2 -d "=" |sed 's/ //g')
echo -e "$id\t$cs12\t$nl12\t$gg12\t$pt12\t$pp12\t$hs12\t$pa12\t$lnL12\t$np12\t$ntime12\t$kappa12" >> "free_ratio_${cod2}_model_dnds.tsv"  # create a table file containing dN/dS and other information

################################21. two ratio	& codon freqeuncy model="FX61############################################################################
sed -i "s/iiii/${id}.cds.aln/g" "${id}_${mod2}_${cod1}.ctl"
    sed -i "s/oooo/${id}_${mod2}_${cod1}/g" "${id}_${mod2}_${cod1}.ctl"
    sed -i "s/mmmm/${mod2}/g"  "${id}_${mod2}_${cod1}.ctl"
    sed -i "s/cccc/${codon1}/g"  "${id}_${mod2}_${cod1}.ctl" #########PLEASE MANUALY CHECK 
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
echo -e "$id\t$cs21\t$nl21\t$gg21\t$pt21\t$pp21\t$hs21\t$pa21\t$lnL21\t$np21\t$ntime21\t$kappa21" >> "two_ratio_${cod1}_model_dnds.tsv"

################################22. two ratio	& codon freqeuncy model="F3x4"############################################################################

sed -i "s/iiii/${id}.cds.aln/g" "${id}_${mod2}_${cod2}.ctl"
    sed -i "s/oooo/${id}_${mod2}_${cod2}/g" "${id}_${mod2}_${cod2}.ctl"
    sed -i "s/mmmm/${mod2}/g"  "${id}_${mod2}_${cod2}.ctl"
    sed -i "s/cccc/${codon2}/g"  "${id}_${mod2}_${cod2}.ctl" #########PLEASE MANUALY CHECK 
    codeml "${id}_${mod2}_${cod2}.ctl"
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
echo -e "$id\t$cs22\t$nl22\t$gg22\t$pt22\t$pp22\t$hs22\t$pa22\t$lnL22\t$np22\t$ntime22\t$kappa22" >> "two_ratio_${cod2}_model_dnds.tsv"

################################01.MO ratio	& codon freqeuncy model="FX61"#########################################################################################################

sed -i "s/iiii/${id}.cds.aln/g" "${id}_${mod0}_${cod1}.ctl"
    sed -i "s/oooo/${id}_${mod0}_${cod1}/g" "${id}_${mod0}_${cod1}.ctl"
    sed -i "s/mmmm/${mod0}/g"  "${id}_${mod0}_${cod1}.ctl"
    sed -i "s/cccc/${codon1}/g"  "${id}_${mod0}_${cod1}.ctl" #########PLEASE MANUALY CHECK 
    codeml "${id}_${mod0}_${cod1}.ctl"


np01=$(grep "lnL" ${id}_${mod0}_${cod1}.out  |cut -f2 -d ":" | sed 's/ //g' | sed 's/np//g')
ntime01=$(grep "lnL" ${id}_${mod0}_${cod1}.out  |cut -f3 -d ":" | sed 's/ //g'| sed 's/)//g')
lnL01=$(grep "lnL" ${id}_${mod0}_${cod1}.out |sed 's/..*\:\ *//' | sed 's/\ ..*//')
kappa01=$(grep "ts/tv" ${id}_${mod0}_${cod1}.out |cut -f2 -d "=" |sed 's/ //g')
omega01=$(grep "omega (dN/dS)"  ${id}_${mod0}_${cod1}.out |cut -f2 -d "=" |sed 's/ //g')
#np01=$(grep "lnL"  ${id}_${mod0}_${cod2}.out |cut -f2 -d ":" | sed 's/ //g'| sed 's/np//g')

echo -e "$id\t$omega01\t$lnL01\t$np01\t$ntime01\t$kappa01" >> "ono_ratio_${cod1}_model_dnds.tsv" 

################################02.MO rati1	& codon freqeuncy model="F3X4"#################################################################################################
 sed -i "s/iiii/${id}.cds.aln/g" "${id}_${mod0}_${cod2}.ctl"
    sed -i "s/oooo/${id}_${mod0}_${cod2}/g" "${id}_${mod0}_${cod2}.ctl"
    sed -i "s/mmmm/${mod0}/g"  "${id}_${mod0}_${cod2}.ctl"
    sed -i "s/cccc/${codon2}/g"  "${id}_${mod0}_${cod2}.ctl" #########PLEASE MANUALY CHECK 
    codeml "${id}_${mod0}_${cod2}.ctl"
 



ntime02=$(grep "lnL" ${id}_${mod0}_${cod2}.out  |cut -f3 -d ":" | sed 's/ //g'| sed 's/)//g')
lnL02=$(grep "lnL" ${id}_${mod0}_${cod2}.out |sed 's/..*\:\ *//' | sed 's/\ ..*//')
kappa02=$(grep "ts/tv" ${id}_${mod0}_${cod2}.out |cut -f2 -d "=" |sed 's/ //g')
#omega02=$(grep "omega (dN/dS)"  ${id}_${mod0}_${cod2}.out |cut -f2 -d "=" |sed 's/ //g')
omega02=$(grep "omega (dN/dS)"  ${id}_${mod0}_${cod2}.out |cut -f2 -d "=" |sed 's/ //g')
#np02=$(grep "lnL"  ${id}_${mod0}_${cod2}.out |cut -f2 -d ":" | sed 's/ //g'| sed 's/np//g')
np02=$(grep "lnL" ${id}_${mod0}_${cod2}.out  |cut -f2 -d ":" | sed 's/ //g' | sed 's/np//g')
#echo -e "$id\t$omega02\t$lnL02\t$np02\t$ntime02\t$kappa02" >> "ono_ratio_${cod2}_model_dnds.tsv"
echo -e "$id\t$omega02\t$lnL02\t$np02\t$ntime02\t$kappa02" >> "ono_ratio_${cod2}_model_dnds.tsv"


#############

directories=("protein_aln_prank" "cds_aln_pan2pal" "codeml_models" "codeml_out" "codeml_out/free_ratio/F61" "codeml_out/free_ratio/F3X4" "codeml_out/M0/F61" "codeml_out/M0/F3X4" "codeml_out/two_ratio/F3X4" "codeml_out/two_ratio/F61")

# Loop to create directories only if they don't exist
for dir in "${directories[@]}"; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
 
   fi
done


# Now move respective files to destination folders
mv "$id".pep.aln.best.fas protein_aln_prank/
mv "$id".cds.aln cds_aln_pan2pal/
mv $id**ctl codeml_models/
cp -r "${id}_${mod0}_${cod1}.out" codeml_out/M0/F61/
cp -r "${id}_${mod1}_${cod1}.out" codeml_out/free_ratio/F61/
cp -r "${id}_${mod1}_${cod1}_raw_dnds.tre" codeml_out/free_ratio/F61/
cp -r "${id}_${mod2}_${cod1}.out" codeml_out/two_ratio/F61/
cp -r "${id}_${mod2}_${cod1}_raw_dnds.tre" codeml_out/two_ratio/F61/
cp -r "${id}_${mod0}_${cod2}.out" codeml_out/M0/F3X4/
cp -r "${id}_${mod1}_${cod2}.out" codeml_out/free_ratio/F3X4/
cp -r "${id}_${mod1}_${cod2}_raw_dnds.tre" codeml_out/free_ratio/F3X4/
cp -r "${id}_${mod2}_${cod2}.out" codeml_out/two_ratio/F3X4/
cp -r "${id}_${mod2}_${cod2}_raw_dnds.tre" codeml_out/two_ratio/F3X4/

done
