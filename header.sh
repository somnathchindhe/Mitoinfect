find sequences/ -type f -name "*fas" | while read -r file; do
  while read -r spe; do
    header=$(grep "$spe" "$file")
    new=$(echo ">$spe")
    sed -i "s/$header/$new/g" "$file"
  done < "/home/somnath/MitoHPI/Data/ortho_3975/7_primate_data-analysis/sequence_Mitolink_2898_1to_1_ortholog_7_primate/sequences/new_cds/species_name"
done
