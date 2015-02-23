FROM r-base:3.1.2
MAINTAINER Eamon O'Dea <odea35@gmail.com>

RUN apt-get update && apt-get install -y -q --no-install-recommends \
  curl \
  mafft \
  openssh-server \
  python-biopython \
  r-cran-rcurl \
  r-cran-xml

RUN curl -SL http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_Linux64_0.91b.tar.Z \
  | tar -xzC /opt

RUN mkdir /var/run/sshd && echo 'root:screencast' | chpasswd \
  && sed -i 's/PermitRootLogin without-password/PermitRootLogin yes/' /etc/ssh/sshd_config \
  && sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd
EXPOSE 22

RUN install2.r --error \
  ape \
  lubridate \
  rentrez \
  && rm -rf /tmp/download_packages/ /tmp/*.rds

CMD ["/usr/sbin/sshd", "-D"]
