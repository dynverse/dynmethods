FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'devtools::install_github("dynverse/dynwrap")'



RUN R -e 'install.packages("scaterlegacy")'
RUN R -e 'devtools::install_github("dynverse/embeddr")'

ADD . /code

ENTRYPOINT Rscript /code/run.R