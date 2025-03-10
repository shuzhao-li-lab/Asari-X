import logging, datetime

def setup_logger():
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    filename = f'log_asarix.log'
    logging.basicConfig(filename=filename, level=logging.INFO)
    return logging.getLogger(__name__)