FROM library/python

RUN pip install scipy numpy pandas matplotlib

RUN pip install git+https://github.com/kieranrcampbell/ouijaflow.git --upgrade --upgrade-strategy only-if-needed
RUN pip install tensorflow==1.6 # temporary fix for edward https://github.com/blei-lab/edward/issues/882

ADD . /code
ENTRYPOINT python /code/run.py
