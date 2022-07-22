import os
import sys
import datetime

from System.MSData import CDataPack
from System.MSLogging import create_logger,log_to_user
from System.MSSystem import MS2_OUTPUT_TYPE, PSM_OUTPUT_TYPE,CFLOW6_INFORMATION

from Function.MSFunctionIO import CFunctionFastaFileParser, CFunctionCxResultFileParser, CFunctionDownloadUniprotInfos, \
    CFunctionParseUniprot, CFunctionDownloadPDBInfos,CFunctionOutPut
from Function.MSFunctionStructure import CFunctionPPiPDB


class CTaskPreparing:
    def work(self, data_package=CDataPack()):
        file_name = self.__clearLogFile(data_package)
        self.__createLocalLogger(file_name, data_package)

    def __clearLogFile(self, data_package=CDataPack()):
        os.chdir(os.path.normpath(os.path.dirname(os.path.dirname(__file__))))
        file_name = os.path.join(data_package.my_config.E1_PATH_EXPORT, 'ComMap.log')
        data_package.my_config.E6_LOCAL_LOG_FILE_PATH = file_name
        with open(file_name, 'w') as f:
            f.write('\n')
        return file_name

    def __createLocalLogger(self, file_name, data_package=CDataPack()):
        data_package.my_config.E5_LOCAL_LOGGER = create_logger('local_logger', file_name)



class CTaskReadCxInput:
    def work(self, data_package=CDataPack()):
        self.__readFasta(data_package)
        self.__readCxResult(data_package)

    def __readFasta(self, data_package=CDataPack()):
        fasta_data = CFunctionFastaFileParser()
        fasta_data.read_fasta(data_package, data_package.my_ProteinIndex)

    def __readCxResult(self, data_package=CDataPack()):
        cx_data = CFunctionCxResultFileParser()
        cx_data.read_cx_result(data_package)



class CTaskDownloadUniprot:
    def work(self, data_package=CDataPack()):
        self.__getProteinSet(data_package)
        self.__checkAlreadyDownloaded(data_package)
        self.__downLoadUniprot(data_package)

    def __getProteinSet(self, data_package=CDataPack()):
        return_list = []
        for ppi in data_package.my_CxResult.CFLOW6_PPI:
            return_list.append(ppi.protein1)
            return_list.append(ppi.protein2)
        data_package.my_config.C62_UNIQUE_PROTEINS = list(set(return_list))

    def __checkAlreadyDownloaded(self, data_package):
        return_list = []
        for protein in data_package.my_config.C62_UNIQUE_PROTEINS:
            if protein + '.xml' not in os.listdir(data_package.my_config.C60_TYPE01_DATA_PATH):
                data_package.my_config.C62_PROTIENS_TO_DOWNLOAD.append(protein)
            else:
                data_package.my_config.C62_PROTEINS_TO_PDB.append(protein)

    def __downLoadUniprot(self, data_package):
        download_func = CFunctionDownloadUniprotInfos(proteins_list=data_package.my_config.C62_PROTIENS_TO_DOWNLOAD,
                                                      concurrent_threads=data_package.my_config.P1_NUMBER_THREAD)
        all_success_proteins, all_failure_proteins = download_func.robust_download_protein_info(
            retry_time=data_package.my_config.C62_MAX_DOWNLOAD_TRY_TIME,folder=data_package.my_config.C60_TYPE01_DATA_PATH)
        for item in all_success_proteins:
            data_package.my_config.C62_PROTEINS_TO_PDB.append(item)

class CTaskPrepareLocalUniprot:
    def work(self,data_package=CDataPack()):
        pass


class CTaskDownloadPDB:
    def work(self, data_package=CDataPack()):
        self.__getPDBInfomation(data_package)
        self.__checkAlreadyDownloaded(data_package)
        self.__downLoadPDB(data_package)


    def __getPDBInfomation(self, data_package=CDataPack()):
        parse_function = CFunctionParseUniprot()
        parse_function.buildUniprotPDBMap(data_package.my_config.C62_PROTEINS_TO_PDB, data_package)
        ppi_pdb_function = CFunctionPPiPDB()
        ppi_pdb_function.calculatePPiPDB(data_package)
        pdbs = [item2 for item1 in data_package.my_CxResult.CFLOW6_PPI_PDB for item2 in item1]
        data_package.my_config.C63_UNIQUE_PDBS = set(pdbs)

    def __checkAlreadyDownloaded(self, data_package):
        return_list = []
        for pdb_name in data_package.my_config.C63_UNIQUE_PDBS:
            if pdb_name.lower() + '.cif.gz' not in os.listdir(data_package.my_config.C60_TYPE01_DATA_PATH):
                data_package.my_config.C63_PDBS_TO_DOWNLOAD.append(pdb_name)
            else:
                data_package.my_config.C63_SUCCESS_PDBS.append(pdb_name)

    def __downLoadPDB(self, data_package):
        pdb_tasks = CFunctionDownloadPDBInfos(file_type='cif',
                                              proteins_list=data_package.my_config.C63_PDBS_TO_DOWNLOAD,
                                              concurrent_threads=data_package.my_config.P1_NUMBER_THREAD)
        all_success_proteins, all_failure_proteins = pdb_tasks.robust_download_protein_info(
            retry_time=data_package.my_config.C62_MAX_DOWNLOAD_TRY_TIME,folder=data_package.my_config.C60_TYPE01_DATA_PATH)
        data_package.my_config.C63_SUCCESS_PDBS.extend(all_success_proteins)
        if data_package.my_config.C63_SUCCESS_PDBS==[]:
            log_to_user(data_package.my_config.E5_LOCAL_LOGGER, CFLOW6_INFORMATION[12])
            date_now = datetime.datetime.now()
            log_to_user(data_package.my_config.E5_LOCAL_LOGGER, date_now)
            sys.exit(CFLOW6_INFORMATION[12])
        else:
            data_package.my_config.C63_SUCCESS_PDBS=[item.upper() for item in data_package.my_config.C63_SUCCESS_PDBS]
        num_max=len(data_package.my_CxResult.CFLOW6_PPI_PDB)
        for item in all_failure_proteins:
            for j in range(num_max):
                if item in data_package.my_CxResult.CFLOW6_PPI_PDB[j]:
                    data_package.my_CxResult.CFLOW6_PPI_PDB[j].remove(item)

class CTaskPrepareLocalPDB:

    def work(self,data_package=CDataPack()):
        map_result=self.__getNamePdbMap(data_package)
        self.__scanLocalPdb(data_package,map_result)
        for protein in data_package.my_config.C63_PROTEINS_PDB_MAP:
            for pdb_name in data_package.my_config.C63_PROTEINS_PDB_MAP[protein]:
                if pdb_name in data_package.my_config.C63_PDB_PROTEINS_MAP:
                    data_package.my_config.C63_PDB_PROTEINS_MAP[pdb_name].append(protein)
                else:
                    data_package.my_config.C63_PDB_PROTEINS_MAP[pdb_name]=[protein]
        for _ in range(len(data_package.my_CxResult.CFLOW6_PPI)):
            data_package.my_CxResult.CFLOW6_PPI_PDB.append(data_package.my_config.C63_SUCCESS_PDBS)


    def __getNamePdbMap(self,data_package):
        map_result=[[],{}]
        try:
            with open(data_package.my_config.C60_TYPE0_PROTEIN_STRUCTURE_MAP,'r') as f:
                for line in f:
                    protein_name,pdb_name=line.strip().split()
                    protein_name,pdb_name=protein_name.upper(),pdb_name.upper()
                    if protein_name and pdb_name:
                        map_result[0].append(pdb_name)
                        if protein_name in map_result[1]:
                            map_result[1][protein_name].append(pdb_name)
                        else:
                            map_result[1][protein_name]=[pdb_name]
        except:
            sys.exit(CFLOW6_INFORMATION[14])
        data_package.my_config.C63_PROTEINS_PDB_MAP=map_result[1]
        return map_result

    def __scanLocalPdb(self,data_package,map_result):
        for item in os.listdir(data_package.my_config.C60_TYPE0_PDB_STRUCTURE_PATH):
            if 'cif' in item and item.split('.')[0].upper() in map_result[0]:
                data_package.my_config.C63_SUCCESS_PDBS.append(item.split('.')[0].upper())
        if data_package.my_config.C63_SUCCESS_PDBS==[]:
            log_to_user(data_package.my_config.E5_LOCAL_LOGGER, CFLOW6_INFORMATION[12])
            date_now = datetime.datetime.now()
            log_to_user(data_package.my_config.E5_LOCAL_LOGGER, date_now)
            sys.exit(CFLOW6_INFORMATION[12])


class CTaskOutPut:

    def work(self,data_package=CDataPack()):
        output_task=CFunctionOutPut()
        output_task.writeTotalResult(data_package)
        output_task.writePDBResult(data_package)
        output_task.writeBinary(data_package)


class CTaskReadCxBin:
    def work(self,data_package=CDataPack()):
        output_task=CFunctionOutPut()
        output_task.readBinary(data_package)