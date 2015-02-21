FROM r-base:3.1.2
MAINTAINER Eamon O'Dea <odea35@gmail.com>

RUN apt-get update && apt-get install -y -q --no-install-recommends \
  mafft \
  openssh-server \
  r-cran-rcurl \
  r-cran-xml

RUN mkdir /var/run/sshd && echo 'root:screencast' | chpasswd \
&& sed -i 's/PermitRootLogin without-password/PermitRootLogin yes/' /etc/ssh/sshd_config \
&& sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd
EXPOSE 22

RUN install2.r --error \
  rentrez \
&& rm -rf /tmp/download_packages/ /tmp/*.rds

CMD ["/usr/sbin/sshd", "-D"]
