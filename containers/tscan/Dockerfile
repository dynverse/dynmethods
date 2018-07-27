FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'install.packages("TSCAN")'

ADD . /code

ENTRYPOINT Rscript /code/run.R
