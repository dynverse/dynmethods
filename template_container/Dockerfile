# R
FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile
RUN R -e 'install.packages("template")'
RUN R -e 'devtools::install_github("dynverse/template")'

# Python
FROM library/python
RUN pip install library1 library2

ADD . /code

ENTRYPOINT cat /code/run.R | R --no-save
ENTRYPOINT python /code/run.py
