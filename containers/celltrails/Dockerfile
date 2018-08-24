FROM dynverse/dynwrap:bioc

RUN R -e 'devtools::install_github("dcellwanger/CellTrails")'

ADD . /code
ENTRYPOINT /code/run.sh
