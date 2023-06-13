#!/bin/bash
#
#SBATCH --job-name=idr
#SBATCH --output=slurm-fin_pk-%j.o.out
#SBATCH --error=slurm-fin_pk-%j.e.out
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8


### Se crea una carpeta temporal, que se borrara al final del programa, para almacenar el input y output: SCRATCH
SCRATCH=/scratch/cvicgar$SLURM_JOB_ID

mkdir -p $SCRATCH || exit $?

### Guardo en run una carpeta por punto, con los .bed de las replicas. Se recoge el nombre de los archivos y los archivos mismos.
### Ademas se recoge el working directory para grabar los datos.

wd=$PWD
folder=$1
name_folder=${folder:0:(-1)}
file1=$folder/*1_nucfree.bed
file2=$folder/*2_nucfree.bed

### Cuando se llama al programa se le da el tamanio del genoma como argumento. 

gsize=$2

### Se copian los datos requeridos en el archivo temporal SCRATCH

cp $file1 $SCRATCH/1.bed
cp $file2 $SCRATCH/2.bed

### Se crea la carpeta results en SCRATCH. 

cd $SCRATCH
mkdir results
cd results

### Peak calling en las replicas y en el pool de replicas: IDR conservative peaks --> high confidence peaks, that represent
### events across true biological replciates and that account for true biological and technical noise.

#activate the virtual env for python 2 which has MACS2


echo 'Starting proper replicate peak calling'
mkdir CONS_PEAKS
macs2 callpeak -f BED -t ../1.bed --nomodel --extsize 100 --shift -45  -p 0.02 -g ${gsize} -n ${name_folder}_rep1 --outdir ./CONS_PEAKS &
macs2 callpeak -f BED -t ../2.bed --nomodel --extsize 100 --shift -45  -p 0.02 -g ${gsize} -n ${name_folder}_rep2 --outdir ./CONS_PEAKS &
macs2 callpeak -f BED -t ../1.bed ../2.bed --nomodel --extsize 100 --shift -45  -p 0.02 -g ${gsize} -n ${name_folder}_pooled_rep --outdir ./CONS_PEAKS &
wait %1 %2 %3
echo 'Done with proper replicate peak calling'

### Se filtran los 5E5 mejores picos (p-value en 8 columna) --> nr indica numerico y orden reverso

echo 'Extracting top 500000 peaks'
sort -k8,8nr ./CONS_PEAKS/${name_folder}_rep1_peaks.narrowPeak | head -n 500000 > ./CONS_PEAKS/${name_folder}_rep1_toppeaks.bed &
sort -k8,8nr ./CONS_PEAKS/${name_folder}_rep2_peaks.narrowPeak | head -n 500000 > ./CONS_PEAKS/${name_folder}_rep2_toppeaks.bed &
sort -k8,8nr ./CONS_PEAKS/${name_folder}_pooled_rep_peaks.narrowPeak | head -n 500000 > ./CONS_PEAKS/${name_folder}_pool_toppeaks.bed &
wait %1 %2 %3
echo 'Finished extracting top 500000 peaks'

### Un segundo set de picos son los denominados picos optimos --> high-confidence peaks, that represent
### reproducible events and that accounting for read sampling noise. The optimal set is more sensitive, 
### especially when one of the replicates has lower data quality than the other.
### En este caso se han de comparar pseudoreplicas. Para cada replica biologica, se obtienen 2 pseudoreplicas
### para lo cual se barajan los datos de dicha replica y se divide en dos archivos con el mismo numero de lineas.
### Los datos que se barajan son los originales, los nucfree.bed, que contiene coordenadas de reads, pero estan de una en una,
### no agrupadas con su dato de frecuencia (como hace MACS2). Por eso es posible barajar y dividir por la mitad: si una zona
### tiene un pico de verdad, al dividir las reads en dos, si es suficientemente fuerte, los dos ficheros tendran suficientes 
### reads en esa zona como para dar lugar a un pico.

## Creacion de pseudoreplicas

echo 'Creating pseudoreplicates'
mkdir OPT_PEAKS

nlines1=$( cat ../1.bed | wc -l ) 
halfnlines1=$(( (nlines1 + 1) / 2 )) 
cat ../1.bed | shuf | split -d -l ${halfnlines1} - ./OPT_PEAKS/${name_folder}_rep1_ps  ## -d indica que se pondran sufijos numericos

nlines2=$( cat ../2.bed | wc -l ) 
halfnlines2=$(( (nlines2 + 1) / 2 )) 
cat ../2.bed | shuf | split -d -l ${halfnlines2} - ./OPT_PEAKS/${name_folder}_rep2_ps

cat ../1.bed ../2.bed | shuf > ./pooled.bed
nlinesp=$( cat ./pooled.bed | wc -l ) 
halfnlinesp=$(( (nlinesp + 1) / 2 )) 
cat ./pooled.bed | shuf | split -d -l ${halfnlinesp} - ./OPT_PEAKS/${name_folder}_pool_ps 

echo 'Done creating pseudoreplicates'

## Peak calling de las pseudoreplicas

echo 'Starting pseudoreplicate peak calling'
macs2 callpeak -f BED -t ./OPT_PEAKS/${name_folder}_rep1_ps00 --nomodel --extsize 100 --shift -45  -p 0.02 -g ${gsize} -n ${name_folder}_rep1_ps00 --outdir ./OPT_PEAKS &
macs2 callpeak -f BED -t ./OPT_PEAKS/${name_folder}_rep1_ps01 --nomodel --extsize 100 --shift -45  -p 0.02 -g ${gsize} -n ${name_folder}_rep1_ps01 --outdir ./OPT_PEAKS &

macs2 callpeak -f BED -t ./OPT_PEAKS/${name_folder}_rep2_ps00 --nomodel --extsize 100 --shift -45  -p 0.02 -g ${gsize} -n ${name_folder}_rep2_ps00 --outdir ./OPT_PEAKS &
macs2 callpeak -f BED -t ./OPT_PEAKS/${name_folder}_rep2_ps01 --nomodel --extsize 100 --shift -45  -p 0.02 -g ${gsize} -n ${name_folder}_rep2_ps01 --outdir ./OPT_PEAKS &

macs2 callpeak -f BED -t ./OPT_PEAKS/${name_folder}_pool_ps00 --nomodel --extsize 100 --shift -45  -p 0.02 -g ${gsize} -n ${name_folder}_pool_ps00 --outdir ./OPT_PEAKS &
macs2 callpeak -f BED -t ./OPT_PEAKS/${name_folder}_pool_ps01 --nomodel --extsize 100 --shift -45  -p 0.02 -g ${gsize} -n ${name_folder}_pool_ps01 --outdir ./OPT_PEAKS &

wait %1 %2 %3 %4 %5 %6
echo 'Done with pseudoreplicate peak calling'

echo 'Extracting top 500000 peaks'
sort -k8,8nr ./OPT_PEAKS/${name_folder}_rep1_ps00_peaks.narrowPeak | head -n 500000 > ./OPT_PEAKS/${name_folder}_rep1_ps00_toppeaks.bed &
sort -k8,8nr ./OPT_PEAKS/${name_folder}_rep1_ps01_peaks.narrowPeak | head -n 500000 > ./OPT_PEAKS/${name_folder}_rep1_ps01_toppeaks.bed &

sort -k8,8nr ./OPT_PEAKS/${name_folder}_rep2_ps00_peaks.narrowPeak | head -n 500000 > ./OPT_PEAKS/${name_folder}_rep2_ps00_toppeaks.bed &
sort -k8,8nr ./OPT_PEAKS/${name_folder}_rep2_ps01_peaks.narrowPeak | head -n 500000 > ./OPT_PEAKS/${name_folder}_rep2_ps01_toppeaks.bed &

sort -k8,8nr ./OPT_PEAKS/${name_folder}_pool_ps00_peaks.narrowPeak | head -n 500000 > ./OPT_PEAKS/${name_folder}_pool_ps00_toppeaks.bed &
sort -k8,8nr ./OPT_PEAKS/${name_folder}_pool_ps01_peaks.narrowPeak | head -n 500000 > ./OPT_PEAKS/${name_folder}_pool_ps01_toppeaks.bed &

wait %1 %2 %3 %4 %5 %6
echo 'Finished extracting top 500000 peaks'

### IDR corre en py3 asi que hay que cambiar de entorno. 
### Se realiza el analisis de las replicas originales y de cada par de pseudoreplicas.

#de-activate the virtual env for python 2 which has MACS2
#deactivate
#activate the virtual env for python 3 which has IDR
#source /home/ska/tools/py3_venv/bin/activate

echo 'Starting IDR analysis'
idr --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.1 -s ./CONS_PEAKS/${name_folder}_rep1_toppeaks.bed ./CONS_PEAKS/${name_folder}_rep2_toppeaks.bed -p ./CONS_PEAKS/${name_folder}_pool_toppeaks.bed --verbose -l ./CONS_PEAKS/idr_consrv.log -o ./CONS_PEAKS/${name_folder}_rep1_rep2.txt --plot &
idr --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.1 -s ./OPT_PEAKS/${name_folder}_rep1_ps00_toppeaks.bed ./OPT_PEAKS/${name_folder}_rep1_ps01_toppeaks.bed -p ./CONS_PEAKS/${name_folder}_rep1_toppeaks.bed --verbose -l ./OPT_PEAKS/idr_opt_rep1.log -o ./OPT_PEAKS/${name_folder}_rep1_ps.txt --plot &
idr --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.1 -s ./OPT_PEAKS/${name_folder}_rep2_ps00_toppeaks.bed ./OPT_PEAKS/${name_folder}_rep2_ps01_toppeaks.bed -p ./CONS_PEAKS/${name_folder}_rep2_toppeaks.bed --verbose -l ./OPT_PEAKS/idr_opt_rep2.log -o ./OPT_PEAKS/${name_folder}_rep2_ps.txt --plot &
idr --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.1 -s ./OPT_PEAKS/${name_folder}_pool_ps00_toppeaks.bed ./OPT_PEAKS/${name_folder}_pool_ps01_toppeaks.bed -p ./CONS_PEAKS/${name_folder}_pool_toppeaks.bed --verbose -l ./OPT_PEAKS/idr_opt_pool.log -o ./OPT_PEAKS/${name_folder}_pool_ps.txt --plot &
wait %1 %2 %3 %4
echo 'Done IDR analysis'

#deactivate

### Se filtran los picos por un valor de IDR global de 0.01 (columna 12 >= 1 puesto que -log10(0.1)) para las replicas originales.
### Informacion sobre el output en https://github.com/nboley/idr

awk 'BEGIN{OFS="\t"} $12>='"1"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ./CONS_PEAKS/${name_folder}_rep1_rep2.txt | sort -k7n,7n  > ./CONS_PEAKS/${name_folder}_idrConsPeaks.bed &
awk 'BEGIN{OFS="\t"} $12>='"1"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ./OPT_PEAKS/${name_folder}_rep1_ps.txt | bedtools intersect -v -a - -b ./CONS_PEAKS/${name_folder}_idrConsPeaks.bed -f 0.55 | sort -k7n,7n  > ./OPT_PEAKS/${name_folder}_idrOptPeaks_intsrep1.bed &
awk 'BEGIN{OFS="\t"} $12>='"1"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ./OPT_PEAKS/${name_folder}_rep1_ps.txt | sort -k7n,7n  > ./OPT_PEAKS/${name_folder}_idrOptPeaks_rep1.bed &
awk 'BEGIN{OFS="\t"} $12>='"1"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ./OPT_PEAKS/${name_folder}_rep2_ps.txt | bedtools intersect -v -a - -b ./CONS_PEAKS/${name_folder}_idrConsPeaks.bed -f 0.55 | sort -k7n,7n  > ./OPT_PEAKS/${name_folder}_idrOptPeaks_intsrep2.bed &
awk 'BEGIN{OFS="\t"} $12>='"1"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ./OPT_PEAKS/${name_folder}_rep2_ps.txt | sort -k7n,7n  > ./OPT_PEAKS/${name_folder}_idrOptPeaks_rep2.bed &
awk 'BEGIN{OFS="\t"} $12>='"1"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ./OPT_PEAKS/${name_folder}_pool_ps.txt | bedtools intersect -v -a - -b ./CONS_PEAKS/${name_folder}_idrConsPeaks.bed -f 0.55 | sort -k7n,7n  > ./OPT_PEAKS/${name_folder}_idrOptPeaks_intspool.bed &
awk 'BEGIN{OFS="\t"} $12>='"1"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ./OPT_PEAKS/${name_folder}_pool_ps.txt | sort -k7n,7n  > ./OPT_PEAKS/${name_folder}_idrOptPeaks_pool.bed &
wait %1 %2 %3 %4 %5 %6 %7

### Se obtienen estadisticas sobre numero de picos.

conspeaks=($(wc -l ./CONS_PEAKS/${name_folder}_rep1_rep2.txt))
conspeaks_filt=($(wc -l ./CONS_PEAKS/${name_folder}_idrConsPeaks.bed))

optpeaks1=($(wc -l ./OPT_PEAKS/${name_folder}_rep1_ps.txt))
optpeaks1_filt_int=($(wc -l ./OPT_PEAKS/${name_folder}_idrOptPeaks_intsrep1.bed))
optpeaks1_filt=($(wc -l ./OPT_PEAKS/${name_folder}_idrOptPeaks_rep1.bed))

optpeaks2=($(wc -l ./OPT_PEAKS/${name_folder}_rep2_ps.txt))
optpeaks2_filt_int=($(wc -l ./OPT_PEAKS/${name_folder}_idrOptPeaks_intsrep2.bed))
optpeaks2_filt=($(wc -l ./OPT_PEAKS/${name_folder}_idrOptPeaks_rep2.bed))

optpeaksp=($(wc -l ./OPT_PEAKS/${name_folder}_pool_ps.txt))
optpeaksp_filt_int=($(wc -l ./OPT_PEAKS/${name_folder}_idrOptPeaks_intspool.bed))
optpeaksp_filt=($(wc -l ./OPT_PEAKS/${name_folder}_idrOptPeaks_pool.bed))

echo 'Writing statistics'
echo "Number of conservative peaks: $conspeaks" > idr_peaks_summary.txt
echo "Number of filtered conservative peaks (GlobalIDR >=1) (Nt): $conspeaks_filt" >> idr_peaks_summary.txt
echo "Number of optimal peaks for rep1: $optpeaks1" >> idr_peaks_summary.txt
echo "Number of filtered optimal peaks for rep1 (GlobalIDR >=1) (N1): $optpeaks1_filt" >> idr_peaks_summary.txt
echo "Number of intersected filtered optimal peaks for rep1 (GlobalIDR >=1): $optpeaks1_filt_int" >> idr_peaks_summary.txt
echo "Number of optimal peaks for rep2: $optpeaks2" >> idr_peaks_summary.txt
echo "Number of filtered optimal peaks for rep2 (GlobalIDR >=1) (N2): $optpeaks2_filt" >> idr_peaks_summary.txt
echo "Number of intersected filtered optimal peaks for rep2 (GlobalIDR >=1): $optpeaks2_filt_int" >> idr_peaks_summary.txt
echo "Number of optimal peaks for pooled data: $optpeaksp" >> idr_peaks_summary.txt
echo "Number of filtered optimal peaks for pooled data (GlobalIDR >=1) (Np): $optpeaksp_filt" >> idr_peaks_summary.txt
echo "Number of intersected filtered optimal peaks for pooled data (GlobalIDR >=1): $optpeaksp_filt_int" >> idr_peaks_summary.txt
max_rr=$(($optpeaksp_filt > $conspeaks_filt ? $optpeaksp_filt : $conspeaks_filt))
min_rr=$(($optpeaksp_filt < $conspeaks_filt ? $optpeaksp_filt : $conspeaks_filt))
echo "Rescue ratio = max (Np, Nt) / min (Np, Nt) :" >> idr_peaks_summary.txt
awk "BEGIN {print (($max_rr/$min_rr))}" >> idr_peaks_summary.txt
max_sc=$(($optpeaks1_filt > $optpeaks2_filt ? $optpeaks1_filt : $optpeaks2_filt))
min_sc=$(($optpeaks1_filt < $optpeaks2_filt ? $optpeaks1_filt : $optpeaks2_filt))
echo "Self consistency ratio = max (N1, N2) / min (N1, N2) = $sc_ratio" >> idr_peaks_summary.txt
awk "BEGIN {print (($max_sc/$min_sc))}" >> idr_peaks_summary.txt

mkdir ./TRACKS
mv ./CONS_PEAKS/${name_folder}_idrConsPeaks.bed ./TRACKS
mv ./OPT_PEAKS/${name_folder}_idrOptPeaks* ./TRACKS

echo 'DONE'

### Se guardan los datos en las carpetas correspondientes y se elimina la carpeta temporal

dest=${wd}/${folder}
cd ..
mv results/* $dest

rm -rf $SCRATCH || exit $?
