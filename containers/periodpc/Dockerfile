FROM dynverse/dynwrap:r

RUN R -e 'devtools::install_cran("princurve")'

ADD . /code

ENTRYPOINT /code/run.sh
