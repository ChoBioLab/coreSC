FROM public.ecr.aws/c8r7r8g5/chobiolab/tools/coresc:0.1.7-beta

WORKDIR /app

ARG APP_PATH
ENV APP_PATH=$APP_PATH

ADD params.csv /app/config/params.csv
ADD samples.csv /app/config/samples.csv

RUN mkdir /app/output \
    && chmod 777 /app/output
