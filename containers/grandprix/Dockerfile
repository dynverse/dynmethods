FROM library/python
RUN pip install tensorflow

RUN git clone https://github.com/GPflow/GPflow.git
RUN cd GPflow; \
    pip install .; \
    cd ../;

RUN git clone https://github.com/ManchesterBioinference/GrandPrix
RUN cd GrandPrix; \
    python setup.py install; \
    cd ../;

ADD . /code
ENTRYPOINT python /code/run.py
