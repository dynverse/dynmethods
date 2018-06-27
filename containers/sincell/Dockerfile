FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'devtools::install_github("dynverse/dynwrap")'



RUN R -e "install.packages('sincell')"

ADD . /code

ENTRYPOINT Rscript /code/run.R