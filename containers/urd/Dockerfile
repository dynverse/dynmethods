FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'install.packages("destiny")'

RUN R -e 'devtools::install_github("farrellja/URD")'

ADD . /code

ENTRYPOINT Rscript /code/run.R
