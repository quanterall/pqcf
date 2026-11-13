# PQCF with the variant-calling tools, so every PQCF worker can run a genome slice itself: minimap2 (alignment),
# samtools and bcftools (sorting, calling), tabix, and RTG Tools' vcfeval (evaluation against the GIAB truth set).
#   docker build -t pqcf-dev . && docker build -t pqcf-variants research/variants/docker
FROM pqcf-dev
RUN apt-get update && apt-get install -y --no-install-recommends minimap2=2.24+dfsg-3+b1 samtools=1.16.1-1 \
        bcftools=1.16-1 tabix=1.16+ds-3 openjdk-17-jre-headless unzip pigz && rm -rf /var/lib/apt/lists/*
RUN curl -sL -o /tmp/rtg.zip https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.12.1/rtg-tools-3.12.1-nojre.zip \
 && unzip -q /tmp/rtg.zip -d /opt && rm /tmp/rtg.zip \
 && printf '#!/bin/sh\n# RTG launcher for aarch64 (the shipped script accepts only x86_64): the jar itself is plain Java\nexec java -Xmx8g -Djava.awt.headless=true -jar /opt/rtg-tools-3.12.1/RTG.jar "$@"\n' > /usr/local/bin/rtg \
 && chmod +x /usr/local/bin/rtg
