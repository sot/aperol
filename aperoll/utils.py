from ska_helpers import logging

logger = logging.basic_logger("aperoll")


class AperollException(RuntimeError):
    pass
