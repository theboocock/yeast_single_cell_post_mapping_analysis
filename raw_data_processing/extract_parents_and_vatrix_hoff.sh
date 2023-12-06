#$ -l h_data=4G,highp -pe shared 16
#$ -l time=336:00:00
#$ -cwd
#$ -M theboocock@ucla.edu
#$ -m a
#$ -e logs/
#$ -o logs/

IN_ROW=$(sed -n ${SGE_TASK_ID}p sge_files/vartrix.txt)
CELLRANGER_FOLDER=$(echo $IN_ROW | awk '{print $1}')
PARENTS=$(echo $IN_ROW | awk '{print $2}')
echo $CELLRANGER_FOLDER
echo $PARENTS
### use single cell environmennt

python extract_parents_and_vatrix.py --cellranger-outdir $CELLRANGER_FOLDER --crosses $PARENTS --vcf /u/home/s/smilefre/project-kruglyak/sceqtl/2021/ref/parents.nostar.vcf.gz
