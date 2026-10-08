FROM sagemath/sagemath:9.7
USER root
RUN apt-get update && \
    apt-get install -y git && \
    apt-get clean

# Sage already provides numpy and scipy; lmfit is only needed by --fit.
RUN sage --pip install --no-cache-dir "lmfit>=1.3"

WORKDIR /app