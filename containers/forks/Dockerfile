FROM dynverse/dynwrap:py3.6

LABEL version 0.1.0.1

RUN pip install seaborn hdbscan

RUN git clone https://github.com/macsharma/FORKS.git

ADD . /code
ENTRYPOINT /code/run.sh
