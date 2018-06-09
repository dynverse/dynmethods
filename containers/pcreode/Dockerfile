FROM library/python:2.7
RUN pip install scipy numpy pandas
RUN git clone git://github.com/KenLauLab/pCreode; \
    cd pCreode; \
    pip install .

ADD . /code
ENTRYPOINT python /code/run.py
