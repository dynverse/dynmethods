FROM dynverse/dynwrap:py3.6

LABEL version 0.1.0.1

RUN pip install tensorflow

RUN git clone https://github.com/GPflow/GPflow.git && cd GPflow && pip install GPflow
RUN git clone https://github.com/ManchesterBioinference/GrandPrix && cd GrandPrix && python setup.py install

ADD . /code

ENTRYPOINT /code/run.sh
