import datetime

from System.MSData import CDataPack
from System.MSFlow import CFlow0, CFlow6, CFlow7
from Function.MSFunctionIO import CFunctionConfigParser
from System.MSLogging import log_to_user, log_get_error, log_get_warning,create_logger
from System.MSSystem import INFO_TO_USER_Staff, EXPIRATION_TIME


class CStaff:
    def __init__(self, input_args):
        self.data_pack = CDataPack()
        self.argv = input_args
        self.global_logger = create_logger('global_logger', '../SpotLink.log')

    def start(self):
        # 设置log对象

        # 版本信息
        log_to_user(self.global_logger,INFO_TO_USER_Staff[0])
        date_now = datetime.datetime.now()
        log_to_user(self.global_logger,date_now)

        # 检查过期
        self.__captain_check_time()

        # 运行流程
        self.__captain_run_flow()

        # 结束流程，打印信息
        log_to_user(self.global_logger,INFO_TO_USER_Staff[4])
        date_now = datetime.datetime.now()
        log_to_user(self.global_logger,date_now)

    def __captain_check_time(self):
        """
        检查软件许可和是否过期时间
        :return:
        """
        date_now = datetime.datetime.now()
        date_dead = datetime.datetime(EXPIRATION_TIME['Year'], EXPIRATION_TIME['Month'], EXPIRATION_TIME['Day'], 23, 59)
        dela_days = (date_dead - date_now).days
        if dela_days < 0:
            log_get_error(INFO_TO_USER_Staff[1], '../SpotLink.log')
        elif dela_days < 7:
            log_get_warning(INFO_TO_USER_Staff[2], '../SpotLink.log')

    def __captain_run_flow(self):
        """
        读取参数以运行不同的工作流
        :return: None
        """
        length_of_args = len(self.argv)
        if length_of_args == 1:
            log_to_user(self.global_logger,INFO_TO_USER_Staff[3])
            flow0 = CFlow0()
            flow0.run()
            log_to_user(self.global_logger,INFO_TO_USER_Staff[6])
        elif length_of_args == 2:
            function_config = CFunctionConfigParser()
            function_config.file_to_config(self.argv[1], self.data_pack.my_config)
            if self.data_pack.my_config.C_TYPE_SEARCH == 6:
                log_to_user(self.global_logger,INFO_TO_USER_Staff[7])
                date_now = datetime.datetime.now()
                log_to_user(self.global_logger,date_now)
                flow6 = CFlow6()
                flow6.run(self.data_pack)
            elif self.data_pack.my_config.C_TYPE_SEARCH == 7:
                log_to_user(self.global_logger,INFO_TO_USER_Staff[7])
                date_now = datetime.datetime.now()
                log_to_user(self.global_logger,date_now)
                flow7 = CFlow7()
                flow7.run(self.data_pack)
            else:
                log_get_error("Get wrong Flow" + str(self.data_pack.my_config.C_TYPE_SEARCH), '../SpotLink.log')
