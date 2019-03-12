############################################################
#                     Dockerfile for R                     #
############################################################
FROM dynverse/dynwrapr:v0.1.0

ARG GITHUB_PAT

RUN R -e 'devtools::install_cran("template")'
RUN R -e 'devtools::install_github("dynverse/template")'

COPY definition.yml run.R example.h5 /code/

ENTRYPOINT ["/code/run.R"]


############################################################
#                  Dockerfile for Python                   #
############################################################
FROM dynverse/dynwrappy:v0.1.0

RUN pip install library1 library2

COPY definition.yml run.py example.h5 /code/

ENTRYPOINT ["/code/run.py"]
