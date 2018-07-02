FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'devtools::install_github("dynverse/dynwrap")'



RUN R -e 'install.packages("rstan")'
RUN R -e 'install.packages("coda")'
RUN R -e 'install.packages("MCMCglmm")'
RUN R -e 'devtools::install_github("dynverse/dyndimred")'
RUN R -e 'devtools::install_github("kieranrcampbell/pseudogp")'

ADD . /code

ENTRYPOINT Rscript /code/run.R