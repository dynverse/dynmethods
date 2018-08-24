FROM dynverse/dynwrap:r

RUN R -e 'devtools::install_cran("SLICER")'

ADD . /code

ENTRYPOINT /code/run.sh
