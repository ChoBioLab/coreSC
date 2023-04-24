FROM public.ecr.aws/chobiolab/seurat:v4-r2

WORKDIR /app
ADD . /app

MAINTAINER christopher.tastad@mssm.edu

ENTRYPOINT ["/app/main"]
