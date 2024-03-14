for i in `cut -f2 Combined_mitoinfect_two_free_ratio_psg_enseble_uniprot.tsv`
do
echo $i
grep $i *_dnds.tsv
done
