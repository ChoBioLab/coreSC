FROM public.ecr.aws/c8r7r8g5/chobiolab/tools/coresc:0.1.7-beta

# edit PARAMS_HOST, SAMPLES_HOST, and OUT_HOST with the appropriate paths

ENV PARAMS_HOST=/ABSOLUTE/HOST/PATH/TO/params.csv
ENV SAMPLES_HOST=/ABSOLUTE/HOST/PATH/TO/params.csv
ENV OUT_HOST=/DESIRED/HOST/OUTPUT/PATH

# do not edit the values below this line
################################################################################

ENV PARAMS_CONTAINER=/app/config/params.csv
ENV SAMPLES_CONTAINER=/app/config/samples.csv
ENV OUT_CONTAINER=/app/output

VOLUME "${PARAMS_HOST}:${PARAMS_CONTAINER}"
VOLUME "${SAMPLES_HOST}:${SAMPLES_CONTAINER}"
VOLUME "${OUT_HOST}:${OUT_CONTAINER}"
