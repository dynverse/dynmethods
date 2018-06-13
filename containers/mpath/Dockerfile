FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'devtools::install_github("dynverse/dynwrap")'



RUN R -e "devtools::install_github('dynverse/Mpath')"
RUN R -e "devtools::install_github('dynverse/Mpath')"

ADD . /code

ENTRYPOINT Rscript code/run.R