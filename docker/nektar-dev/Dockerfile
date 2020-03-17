FROM nektarpp/nektar-env:latest AS build

LABEL maintainer="David Moxey <d.moxey@exeter.ac.uk>"

ARG INSTALL_PREFIX=/usr/local

USER root
COPY --from=nektarpp/nektar:latest ${INSTALL_PREFIX} ${INSTALL_PREFIX}

# Set up entrypoint for copying test files.
COPY docker/nektar/docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN chmod +x /usr/local/bin/docker-entrypoint.sh && \
    ln -s /usr/local/bin/docker-entrypoint.sh /
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

USER nektar
WORKDIR /home/nektar
CMD ["/bin/bash"]
