FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'devtools::install_github("kieranrcampbell/phenopath")'

ADD . /code

ENTRYPOINT Rscript /code/run.R
