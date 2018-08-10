FROM dynverse/dynwrap:py3.6

RUN pip install tensorflow

RUN git clone https://github.com/GPflow/GPflow.git && cd GPflow && pip install GPflow
RUN git clone https://github.com/ManchesterBioinference/GrandPrix && cd GrandPrix && python setup.py install

ADD . /code

ENTRYPOINT python /code/run.py
