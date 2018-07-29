FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'install.packages("template")'

RUN R -e 'devtools::install_github("dynverse/template")'

ADD . /code

ENTRYPOINT Rscript /code/run.R
