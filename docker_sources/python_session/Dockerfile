FROM tensorflow/tensorflow:2.7.0
LABEL authors="Carlos Ruiz-Arenas" \
      description="Docker image containing the python packages for Episignatures project"

## Install python dependencies
RUN pip install -U scikit-learn matplotlib scipy pandas
RUN apt-get update
RUN apt-get install -y git python3.8-venv

RUN git clone https://github.com/arnovanhilten/GenNet /usr/GenNet/

RUN cd /usr/GenNet/ && python3 -m venv env_GenNet
RUN sed -i 's/2.2/2.4/' /usr/GenNet/requirements_GenNet.txt
RUN . /usr/GenNet/env_GenNet/bin/activate && pip3 install -r /usr/GenNet/requirements_GenNet.txt
