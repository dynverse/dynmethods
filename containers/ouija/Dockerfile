FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'devtools::install_github("kieranrcampbell/ouija")'
RUN R -e 'install.packages("rstan")'
RUN R -e 'install.packages("coda")'

ADD . /code

ENTRYPOINT Rscript /code/run.R
