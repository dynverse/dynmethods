FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'devtools::install_github("dynverse/dynwrap")'



RUN R -e "devtools::install_github('dynverse/SLICER')"
RUN R -e "install.packages('lle')"
RUN R -e "install.packages('igraph')"

ADD . /code

ENTRYPOINT Rscript code/run.R