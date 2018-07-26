FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'devtools::install_github("dynverse/dynwrap")'



RUN R -e 'devtools::install_github("dynverse/dynutils")'
RUN R -e 'install.packages("reshape2")'
RUN R -e 'install.packages("destiny")'

ADD . /code

ENTRYPOINT Rscript /code/run.R
