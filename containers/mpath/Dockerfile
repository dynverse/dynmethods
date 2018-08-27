FROM dynverse/dynwrap:r

RUN R -e 'devtools::install_github("dynverse/Mpath")'

# should be
# RUN R -e 'devtools::install_url("https://github.com/JinmiaoChenLab/Mpath/raw/master/Mpath_1.0.tar.gz")'
# but dynverse/Mpath is just so much easier to use than the tar gz above

ADD . /code

ENTRYPOINT /code/run.sh
