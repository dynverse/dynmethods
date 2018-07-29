FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

ADD . /code

ENTRYPOINT Rscript /code/run.R
