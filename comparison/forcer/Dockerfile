# FORM 5 with Forcer (B. Ruijl, T. Ueda, J.A.M. Vermaseren, https://arxiv.org/abs/1704.06650), for the comparison
# in comparison/forcer: the FORM image (comparison/form/Dockerfile) plus Forcer at a fixed commit, and FORM's own
# Forcer and Mincer checks (form-dev/form check/extra). Forcer runs as a separate program; nothing is linked or shipped.
#   docker build -t pqcf-form -f comparison/form/Dockerfile . && docker build -t pqcf-forcer comparison/forcer
FROM pqcf-form
ARG FORCER_REF=68a30bfb2aec3672a2f674ca8666bfc970105e24
RUN git clone -q https://github.com/benruijl/forcer.git /opt/forcer && cd /opt/forcer && git checkout -q "$FORCER_REF" \
 && mkdir -p /opt/formchecks && cd /opt/formchecks \
 && wget -q https://raw.githubusercontent.com/form-dev/form/master/check/extra/forcer.frm \
 && wget -q https://raw.githubusercontent.com/form-dev/form/master/check/extra/mincer.frm
ENV FORMPATH=/opt/forcer
