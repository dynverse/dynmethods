FROM dynverse/dynwrap:py3.6

RUN pip install sklearn seaborn Cython
RUN pip install hdbscan

RUN git clone https://github.com/macsharma/FORKS.git

ADD . /code
ENTRYPOINT /code/run.sh
