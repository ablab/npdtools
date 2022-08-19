import sys
import logging
import config

log = None


def setup_log_handler(handler, level=logging.DEBUG):
    handler.setFormatter(logging.Formatter('%(asctime)s\t%(message)s', '%Y-%m-%d %H:%M:%S'))
    handler.setLevel(level)
    log.addHandler(handler)


def error(msg):
    if config.global_config is not None:
        from common import clean
        clean(config.global_config)
    for handler in list(log.handlers):  # do not print to stdout
        if type(handler) == logging.StreamHandler:
            log.removeHandler(handler)
    setup_log_handler(logging.StreamHandler(sys.stderr), level=logging.ERROR)
    log.error('ERROR: ' + msg)
    log.error("\n\nIn case you have troubles running our tool, you can write to npdtools.support@cab.spbu.ru\n")
    sys.exit(1)


def exception(e):
    if log is not None and log.handlers:
        log.error('')
        log.exception(e)
    else:
        sys.stderr.write(str(e) + '\n')


def info(msg, silent=False):
    if not silent:
        log.info('INFO: ' + msg)


def warning(msg):
    log.warning('WARNING: ' + msg)