FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'devtools::install_github("dynverse/dynwrap")'

RUN apt-get install -y libgsl-dev

RUN R -e "devtools::install_github('dynverse/cellTree@RELEASE_3_6')"

ADD . /code

ENTRYPOINT Rscript code/run.R