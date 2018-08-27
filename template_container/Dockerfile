# R
FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile
RUN R -e 'devtools::install_cran("template")'
RUN R -e 'devtools::install_github("dynverse/template")'

# Python
FROM library/python
RUN pip install library1 library2

# Both
ADD . /code

ENTRYPOINT /code/run.sh
