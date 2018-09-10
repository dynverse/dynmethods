############################################################
#                     Dockerfile for R                     #
############################################################
FROM dynverse/dynwrap:r
# or: FROM dynverse/dynwrap:bioc

RUN R -e 'devtools::install_cran("template")'
RUN R -e 'devtools::install_github("dynverse/template")'

ADD . /code

ENTRYPOINT Rscript /code/run.R


############################################################
#                  Dockerfile for Python                   #
############################################################
FROM dynverse/dynwrap:py3.6
# or FROM dynverse/dynwrap:py2.7

RUN pip install library1 library2

ADD . /code

ENTRYPOINT python /code/run.py
