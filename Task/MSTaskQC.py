from System.MSData import CDataPack
from Function.MSFunctionQC import CFunctionQCRank


class CTaskRank:

    def work(self, data_package=CDataPack()):
        rank_function=CFunctionQCRank()
        rank_function.rank(data_package)