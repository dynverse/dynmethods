FROM dynverse/dynwrap:r

RUN R -e 'install.packages("princurve")'

ADD . /code

ENTRYPOINT /code/run.sh
