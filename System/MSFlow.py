import datetime
import os

from System.MSData import CDataPack
from System.MSLogging import log_to_user
from System.MSSystem import CFLOW6_INFORMATION
from Task.MSTaskIO import CTaskPreparing, CTaskReadCxInput, CTaskDownloadUniprot, CTaskDownloadPDB, \
    CTaskPrepareLocalPDB, CTaskOutPut, CTaskReadCxBin
from Task.MSTaskQC import CTaskRank
from Task.MSTaskMatch import CTaskSequenceAlignment, CTaskStructureDistance


class CFlow0:
    def run(self):
        pass


class CFlow6:

    def run(self, data_package=CDataPack()):
        self.__prepareUserInputs(data_package)
        self.__downloadUniprot(data_package)
        self.__downloadPDB(data_package)
        self.__sequenceAlignment(data_package)
        self.__distanceCalculation(data_package)
        self.__outputAndClean(data_package)
        self.__qualityControl(data_package)

    @staticmethod
    def __prepareUserInputs(data_package):
        prepare_task = CTaskPreparing()
        prepare_task.work(data_package)
        log_to_user(data_package.my_config.E5_LOCAL_LOGGER, CFLOW6_INFORMATION[0])
        date_now = datetime.datetime.now()
        log_to_user(data_package.my_config.E5_LOCAL_LOGGER, date_now)
        input_task = CTaskReadCxInput()
        input_task.work(data_package)

    @staticmethod
    def __downloadUniprot(data_package):
        if data_package.my_config.C60_CALCULATION_TYPE:
            log_to_user(data_package.my_config.E5_LOCAL_LOGGER, CFLOW6_INFORMATION[1])
            date_now = datetime.datetime.now()
            log_to_user(data_package.my_config.E5_LOCAL_LOGGER, date_now)
            download_task = CTaskDownloadUniprot()
            download_task.work(data_package)
            log_to_user(data_package.my_config.E5_LOCAL_LOGGER, CFLOW6_INFORMATION[2])
            date_now = datetime.datetime.now()
            log_to_user(data_package.my_config.E5_LOCAL_LOGGER, date_now)

    @staticmethod
    def __downloadPDB(data_package):
        if data_package.my_config.C60_CALCULATION_TYPE:
            log_to_user(data_package.my_config.E5_LOCAL_LOGGER, CFLOW6_INFORMATION[3])
            date_now = datetime.datetime.now()
            log_to_user(data_package.my_config.E5_LOCAL_LOGGER, date_now)
            download_task = CTaskDownloadPDB()
            download_task.work(data_package)
            log_to_user(data_package.my_config.E5_LOCAL_LOGGER, CFLOW6_INFORMATION[4])
            date_now = datetime.datetime.now()
            log_to_user(data_package.my_config.E5_LOCAL_LOGGER, date_now)
        else:
            log_to_user(data_package.my_config.E5_LOCAL_LOGGER, CFLOW6_INFORMATION[10])
            date_now = datetime.datetime.now()
            log_to_user(data_package.my_config.E5_LOCAL_LOGGER, date_now)
            prepare_task = CTaskPrepareLocalPDB()
            prepare_task.work(data_package)
            log_to_user(data_package.my_config.E5_LOCAL_LOGGER, CFLOW6_INFORMATION[11])
            date_now = datetime.datetime.now()
            log_to_user(data_package.my_config.E5_LOCAL_LOGGER, date_now)

    @staticmethod
    def __sequenceAlignment(data_package):
        log_to_user(data_package.my_config.E5_LOCAL_LOGGER, CFLOW6_INFORMATION[5])
        date_now = datetime.datetime.now()
        log_to_user(data_package.my_config.E5_LOCAL_LOGGER, date_now)
        alignment_task = CTaskSequenceAlignment()
        alignment_task.work_cif(data_package)
        log_to_user(data_package.my_config.E5_LOCAL_LOGGER, CFLOW6_INFORMATION[6])
        date_now = datetime.datetime.now()
        log_to_user(data_package.my_config.E5_LOCAL_LOGGER, date_now)

    @staticmethod
    def __distanceCalculation(data_package):
        log_to_user(data_package.my_config.E5_LOCAL_LOGGER, CFLOW6_INFORMATION[7])
        date_now = datetime.datetime.now()
        log_to_user(data_package.my_config.E5_LOCAL_LOGGER, date_now)
        distance_task = CTaskStructureDistance()
        distance_task.work_cif(data_package)
        log_to_user(data_package.my_config.E5_LOCAL_LOGGER, CFLOW6_INFORMATION[8])
        date_now = datetime.datetime.now()
        log_to_user(data_package.my_config.E5_LOCAL_LOGGER, date_now)

    @staticmethod
    def __outputAndClean(data_package):
        log_to_user(data_package.my_config.E5_LOCAL_LOGGER, CFLOW6_INFORMATION[9])
        date_now = datetime.datetime.now()
        log_to_user(data_package.my_config.E5_LOCAL_LOGGER, date_now)
        output_task = CTaskOutPut()
        output_task.work(data_package)
        if data_package.my_config.C62_CLEAN_CALCULATE:
            for item in data_package.my_config.C63_PDBS_TO_DOWNLOAD:
                os.remove('./Data/{}.ent'.format(item.lower()))

    @staticmethod
    def __qualityControl(data_package):
        log_to_user(data_package.my_config.E5_LOCAL_LOGGER, CFLOW6_INFORMATION[13])
        date_now = datetime.datetime.now()
        log_to_user(data_package.my_config.E5_LOCAL_LOGGER, date_now)
        rank_task=CTaskRank()
        rank_task.work(data_package)


class CFlow7:

    def run(self, data_package=CDataPack()):
        self.__prepareUserInputs(data_package)
        self.__qualityControl(data_package)

    @staticmethod
    def __prepareUserInputs(data_package):
        prepare_task = CTaskPreparing()
        prepare_task.work(data_package)
        log_to_user(data_package.my_config.E5_LOCAL_LOGGER, CFLOW6_INFORMATION[0])
        date_now = datetime.datetime.now()
        log_to_user(data_package.my_config.E5_LOCAL_LOGGER, date_now)
        input_task = CTaskReadCxBin()
        input_task.work(data_package)

    @staticmethod
    def __qualityControl(data_package):
        log_to_user(data_package.my_config.E5_LOCAL_LOGGER, CFLOW6_INFORMATION[13])
        date_now = datetime.datetime.now()
        log_to_user(data_package.my_config.E5_LOCAL_LOGGER, date_now)
        rank_task=CTaskRank()
        rank_task.work(data_package)