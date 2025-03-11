import logging, datetime

def setup_logger():
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    filename = f'log_asarix.log'
    logging.basicConfig(filename=filename, level=logging.INFO, filemode='w')
    logger = logging.getLogger(__name__)
    logging.info(f"instantiating log: {timestamp}")
    return logger