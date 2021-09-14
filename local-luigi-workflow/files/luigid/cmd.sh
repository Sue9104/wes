sudo makedir /etc/luigi && cp luigi.cfg /etc/luigi/luigi.cfg
luigid --port 8082 --pidfile pid --logdir log --background
