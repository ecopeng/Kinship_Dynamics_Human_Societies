#!/bin/sh

wd="$(cd "$(dirname $0)" && pwd)"
dt=$(date '+%Y-%m-%d %H:%M:%S')
echo "[$dt]: GET DEMOGRAPHIC PARAMETERS" > "${wd}/log"

cd ./code/
Rscript "./DEMOGRAPHICS.R" ${wd}

dt=$(date '+%Y-%m-%d %H:%M:%S')
echo "[$dt]: RUN THE MODEL" >> "${wd}/log"

cd ../model_hybrid/
dt=$(date '+%Y-%m-%d %H:%M:%S')
echo "----> $dt: generating predictions by the hybrid model of kinship dynamics"
n_par_set=$(find ./par/ -maxdepth 1 -name 'NDM_*' -printf '.' | wc -c)
echo "----> number of parameter sets: $n_par_set"
n_cores=10
mpic++ -std=c++17 -O3 -o r r.cpp
app="./r"
CMD="mpirun -np $n_cores $app $n_par_set"
mpirun -np $n_cores ./r $n_par_set

cd ../code/
dt=$(date '+%Y-%m-%d %H:%M:%S')
echo "[$dt]: GET THE PREDICTIONS & OBSERVATIONS " >> "${wd}/log"
Rscript "./PREDICTIONS.R" ${wd}
Rscript "./OBSERVATIONS.R" ${wd}

dt=$(date '+%Y-%m-%d %H:%M:%S')
echo "[$dt]: GET THE TEST RESULTS" >> "${wd}/log"
Rscript "./TEST_RELATEDNESS.R" ${wd}
Rscript "./TEST_LOG_RATIO_UNRELATEDNESS.R" ${wd}

dt=$(date '+%Y-%m-%d %H:%M:%S')
echo "[$dt]: PLOT THE PREDICTIONS & OBSERVATIONS & TEST RESULTS & DEMOGRAPHICS" >> "${wd}/log"
Rscript "./PLOT_PREDICTIONS.R" ${wd}
Rscript "./PLOT_OBSERVATIONS.R" ${wd}
Rscript "./PLOT_TEST_RESULTS.R" ${wd}
Rscript "./PLOT_DEMOGRAPHICS.R" ${wd}

cd ../
rm -f Rplots.pdf
cd ./plots/
inkscape \
  --without-gui \
  --file=PRED.pdf \
  --export-plain-svg=PRED.svg
inkscape \
  --without-gui \
  --file=OBSE_KD_DATA.pdf \
  --export-plain-svg=OBSE_KD_DATA.svg
inkscape \
  --without-gui \
  --file=OBSE_KD_DIFF.pdf \
  --export-plain-svg=OBSE_KD_DIFF.svg
inkscape \
  --without-gui \
  --file=TEST.pdf \
  --export-plain-svg=TEST.svg
inkscape \
  --without-gui \
  --file=DEMO.pdf \
  --export-plain-svg=DEMO.svg
