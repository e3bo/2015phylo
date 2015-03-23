FROM r-base:3.1.2
MAINTAINER Eamon O'Dea <odea35@gmail.com>

RUN apt-get update && apt-get install -y -q --no-install-recommends \
  curl \
  git \
  mafft \
  openjdk-8-jre-headless \
  openssh-server \
  python-biopython \
  r-cran-rcurl \
  r-cran-xml

RUN curl -SL http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_Linux64_0.91b.tar.Z \
  | tar -xzC /opt \
  && ln -s /opt/Gblocks_0.91b/Gblocks /usr/bin/Gblocks
RUN curl -SL "http://tree.bio.ed.ac.uk/download.php?id=91&num=3" \
  | tar -xzC /opt \
  && for script in beast treeannotator; \
  do ln -s /opt/BEASTv1.8.1/bin/$script /usr/bin/$script; \
  done
RUN git clone --depth 1 --branch regression git://github.com/e3bo/phast-regression.git /tmp/phast \
  && cd /tmp/phast \
  && R CMD build rPackage \
  && R CMD INSTALL rPackage

RUN mkdir /var/run/sshd && echo 'root:screencast' | chpasswd \
  && sed -i 's/PermitRootLogin without-password/PermitRootLogin yes/' /etc/ssh/sshd_config \
  && sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd
EXPOSE 22

RUN install2.r --error \
  ape \
  knitr \
  ggplot2 \
  lubridate \
  phangorn \
  rentrez \
  && rm -rf /tmp/download_packages/ /tmp/*.rds

CMD ["/usr/sbin/sshd", "-D"]
