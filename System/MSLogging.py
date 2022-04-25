import logging
import sys


def create_logger(name,log_file,level=logging.INFO):
    """
    创建日志logger对象
    :param name: logger对象的名字
    :param log_file: logger文件的路径名字
    :param level: 需要记录的log等级阈值
    :return:
    """
    handler=logging.FileHandler(log_file,mode='a+',encoding='utf-8')
    formatter=logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    logger=logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)
    return logger

def log_to_user(logger,information, print_on_screen=1):
    """
    打印相关信息到文件或窗口台
    :param logger: logger对象
    :param information: 需要打印的信息
    :param print_on_screen: 是否在窗口台上展示相关信息
    :return:
    """
    if print_on_screen:
        print(information)
    logger.info(information)


def log_get_error(information,log_path):
    """
    打印并记录错误
    :param information:
    :return: None
    """
    print(information)
    logging.basicConfig(filename=log_path,
                        filemode='a',
                        format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'
                        )
    logging.error(information)
    sys.exit(0)


def log_get_warning(information,log_path):
    """
    打印并记录警告
    :param information:
    :return: None
    """
    print(information)
    logging.basicConfig(filename=log_path,
                        filemode='a',
                        format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'
                        )
    logging.warning(information)
