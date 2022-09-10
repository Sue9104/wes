import logging
import logging.config
import yaml
import os

def set_logger():
    script_path = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(script_path,'..','logging.yaml'), 'r') as stream:
        config = yaml.load(stream, Loader=yaml.FullLoader)
    logging.config.dictConfig(config)
    rootlogger = logging.getLogger('root')
    cmdlogger = logging.getLogger('cmd')
    timelogger = logging.getLogger('time')
    return [rootlogger,cmdlogger,timelogger]
