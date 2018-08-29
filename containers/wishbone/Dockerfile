FROM dynverse/dynwrap:py3.6

LABEL version 0.1.0.1

RUN pip install git+https://github.com/dynverse/pywishbone --upgrade --upgrade-strategy only-if-needed

ADD . /code
ENTRYPOINT /code/run.sh
