PROJECTDIR="/storage/praha1/home/davidhoksza/projects/apoholo-trj-comp/"

module load pymol

# apo-*, holo-* and outpout params passed as variable
python ${PROJECTDIR}apo-holo-trj-comp.py \
    --apo-gro $apo_gro \
    --apo-trj $apo_trj \
    --apo-ndx $apo_ndx \
    --apo-ix-range $apo_ix_range \
    --holo-gro $holo_gro \
    --holo-trj $holo_trj \
    --holo-ndx $holo_ndx \
    --holo-ix-range $holo_ix_range \
    --blosum $blosum \
    --output $output