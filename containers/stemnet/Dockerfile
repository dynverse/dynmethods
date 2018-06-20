FROM rocker/tidyverse

# install gsl
RUN apt-get install -y libgsl-dev

# install STEMNET
RUN R -e "setRepositories(ind=1:2);devtools::install_git('https://git.embl.de/velten/STEMNET/')"

ADD . /code

ENTRYPOINT Rscript code/run.R
