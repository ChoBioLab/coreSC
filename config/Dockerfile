FROM public.ecr.aws/c8r7r8g5/chobiolab/tools/coresc:0.1.7-beta

ADD params.csv /app/config/params.csv
ADD samples.csv /app/config/samples.csv

RUN mkdir /app/output
RUN usermod -u 1000 /app/output \
    && usermod -G 1000 /app/output
