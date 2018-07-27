FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN apt-get install -y libgsl-dev

RUN R -e "devtools::install_git('https://git.embl.de/velten/STEMNET/')"

ADD . /code

ENTRYPOINT Rscript /code/run.R
