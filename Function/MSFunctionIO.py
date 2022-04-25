import csv
import datetime
import os.path
import gzip
from lxml import etree
import concurrent.futures as cf
import urllib
import requests
import pandas as pd
import pickle

from System.MSData import CDataPack, CProteinIndex, Cppi
from System.MSLogging import log_get_error, log_to_user
from System.MSSystem import RESPONSE_RESULT, UNIPROT_BASE_URL, UNIPROT_FILE_TYPE, PDBJ_PDB_FORMAT_BASE_URL, \
    PDBJ_XML_FILE_TYPE, PDBJ_XML_FORMAT_BASE_URL, PDBJ_PDB_FILE_TYPE, PDBJ_CIF_FILE_TYPE,PDBJ_CIF_FORMAT_BASE_URL,\
    MAX_CONCURRENT_NUMBER,AMINO_ACIDS_DIC

from Tool.MSTool import toolGetWord, toolStr2List

'''
与参数文件相关的函数
'''


class CFunctionConfigParser:
    """
    从文件中读取config文件所存储的参数
    """

    @staticmethod
    def file_to_config(file_path, data_package_my_config):
        with open(file_path, 'r', encoding='utf8') as f:
            for line in f:
                line = line.strip()
                if '=' not in line:
                    continue
                if '#' in line:
                    line = line[:line.find('#')].strip()
                tmp = line.split('=')
                if len(tmp) < 2:
                    continue
                name = tmp[0].strip()
                value = tmp[1]
                if name == 'C_TYPE_SEARCH':
                    data_package_my_config.C_TYPE_SEARCH = int(value.split(';')[0])
                elif name == "I0_INI_PATH_ELEMENT":
                    data_package_my_config.I0_INI_PATH_ELEMENT = eval("r'%s'" % value.split(';')[0].replace(' ', ''))
                elif name == 'I4_INI_PATH_LINKER':
                    data_package_my_config.I4_INI_PATH_LINKER = eval("r'%s'" % value.split(';')[0].replace(' ', ''))
                elif name == 'I5_INI_PATH_ENZYME':
                    data_package_my_config.I5_INI_PATH_ENZYME = eval("r'%s'" % value.split(';')[0].replace(' ', ''))
                elif name == "I1_INI_PATH_AA":
                    data_package_my_config.I1_INI_PATH_AA = eval("r'%s'" % value.split(';')[0].replace(' ', ''))
                elif name == "I2_INI_PATH_MOD":
                    data_package_my_config.I2_INI_PATH_MOD = eval("r'%s'" % value.split(';')[0].replace(' ', ''))
                elif name == "A3_PATH_FASTA":
                    data_package_my_config.A3_PATH_FASTA = eval("r'%s'" % value.split(';')[0].replace(' ', ''))
                elif name == "E1_PATH_EXPORT":
                    data_package_my_config.E1_PATH_EXPORT = eval("r'%s'" % value.split(';')[0].replace(' ', ''))
                    if not os.path.exists(data_package_my_config.E1_PATH_EXPORT):
                        os.makedirs(data_package_my_config.E1_PATH_EXPORT)
                elif name == "P1_NUMBER_THREAD":
                    data_package_my_config.P1_NUMBER_THREAD = int(value.split(';')[0])
                elif name=="C60_CALCULATION_TYPE":
                    data_package_my_config.C60_CALCULATION_TYPE=int(value.split(';')[0])
                    if data_package_my_config.C60_CALCULATION_TYPE:
                        data_package_my_config.C60_TYPE01_DATA_PATH = os.path.normpath(
                            os.path.join(os.path.dirname(os.path.dirname(__file__)), r'Data'))
                elif name=="C60_PDB_STRUCTURE_PATH":
                    data_package_my_config.C60_TYPE0_PDB_STRUCTURE_PATH=eval("r'%s'" % value.replace(' ', ''))
                    data_package_my_config.C60_TYPE01_DATA_PATH =eval("r'%s'" % value.replace(' ', ''))
                elif name=="C60_TYPE0_PROTEIN_STRUCTURE_MAP":
                    data_package_my_config.C60_TYPE0_PROTEIN_STRUCTURE_MAP=eval("r'%s'" % value.replace(' ', ''))
                elif name == "C61_INPUT_LINK_RESULT_TYPE":
                    data_package_my_config.C61_INPUT_LINK_RESULT_TYPE = int(value.split(';')[0])
                elif name == "C61_INPUT_LINK_RESULT_FILE":
                    data_package_my_config.C61_INPUT_LINK_RESULT_FILE = [eval("r'%s'" % item.replace(' ', '')) for item
                                                                         in
                                                                         value.split(';') if item.strip() != '']
                elif name == "C62_MAX_DOWNLOAD_TRY_TIME":
                    data_package_my_config.C62_MAX_DOWNLOAD_TRY_TIME = int(value.split(';')[0])
                elif name == "C62_CLEAN_CALCULATE":
                    data_package_my_config.C62_CLEAN_CALCULATE = int(value.split(';')[0])
                elif name=="C65_LINKER_ARM_LENGTH":
                    data_package_my_config.C65_LINKER_ARM_LENGTH=float(value.split(';')[0])
                elif name=="C65_LINKER_MAX_ARM_LENGTH":
                    data_package_my_config.C65_LINKER_MAX_ARM_LENGTH=float(value.split(';')[0])
                elif name=='P1_NUMBER_THREAD':
                    data_package_my_config.P1_NUMBER_THREAD=int(value.split(';')[0])
                elif name=="FLOW_TYPE":
                    data_package_my_config.C_TYPE_SEARCH=int(value.split(';')[0])
                elif name=="BIN_FILE_PATH":
                    data_package_my_config.C66_BIN_FILE=eval("r'%s'" % value.replace(' ', ''))
                elif name=="C66_OUTPUT_DATE":
                    data_package_my_config.C66_OUTPUT_DATE=value.split(';')[0]


'''
与序列，鉴定结果文件相关的函数
'''


class CFunctionFastaFileParser:
    """
    读取fasta文件至内存以进行后续处理
    """

    def read_fasta(self, data_package=CDataPack(), indexes=CProteinIndex()):
        """
        读取fasta序列，并将其中信息保存到indexes中，同时生成反库信息
        :param data_package:
        :param indexes:
        :return:
        """
        protein_index = self.__proteinIndexWithDescription(data_package)
        if protein_index:
            for index, protein_name in enumerate(protein_index):
                indexes.CFLOW6_SEQUENCE_INDEX[protein_name] = protein_index[protein_name][0]
                indexes.CFLOW6_DESCRIPTION_INDEX[protein_name] = protein_index[protein_name][1]
        else:
            log_get_error('Cannot build protein index.', data_package.my_config.E5_LOG_FILE_PATH)

    @staticmethod
    def __proteinIndexWithDescription(data_package):
        """
        将fasta文件一次性的读入到内存当中
        :return: 一个字典，字典的键为蛋白的名称，字典的值为（蛋白所对应的序列，蛋白对应的描述）
        """
        with open(data_package.my_config.A3_PATH_FASTA, 'r') as file_to_read:
            protein_index_dictionary = {}
            temp_sequence = ''
            temp_protein_name = ''
            temp_protein_description = ''
            for item in file_to_read.readlines():
                if item.strip():
                    if '>' in item:
                        protein_index_dictionary[temp_protein_name] = (temp_sequence, temp_protein_description)
                        temp_protein_name = item.split('|')[1]
                        if item.count('|') >= 2:
                            temp_protein_description = item.split('|')[2]
                        else:
                            temp_protein_description = ''
                        temp_sequence = ''
                    else:
                        temp_sequence += item.strip()
        protein_index_dictionary[temp_protein_name] = (temp_sequence, temp_protein_description)
        del protein_index_dictionary['']
        return protein_index_dictionary


class CFunctionCxResultFileParser:

    def read_cx_result(self, data_package=CDataPack()):
        if data_package.my_config.C61_INPUT_LINK_RESULT_TYPE == 0:
            # 读取txt格式的结果
            for file_path in data_package.my_config.C61_INPUT_LINK_RESULT_FILE:
                self.__read_txt_result(file_path,data_package)
        elif data_package.my_config.C61_INPUT_LINK_RESULT_TYPE == 1:
            # 读取SpotLink的结果
            for file_path in data_package.my_config.C61_INPUT_LINK_RESULT_FILE:
                self.__read_spotlink_psm_result(file_path,data_package)
        elif data_package.my_config.C61_INPUT_LINK_RESULT_TYPE == 2:
            # 读取pLink2的工字形交联结果
            for file_path in data_package.my_config.C61_INPUT_LINK_RESULT_FILE:
                self.__read_plink2_psm_result(file_path, data_package)
        elif data_package.my_config.C61_INPUT_LINK_RESULT_TYPE==3:
            # 读取pLink2的loop交联结果
            for file_path in data_package.my_config.C61_INPUT_LINK_RESULT_FILE:
                self.__read_plink2_loop_result(file_path,data_package)
        elif data_package.my_config.C61_INPUT_LINK_RESULT_TYPE==4:
            # 读取xlinkx的结果
            for file_path in data_package.my_config.C61_INPUT_LINK_RESULT_FILE:
                self.__read_xlinkx_txt_result(file_path,data_package)

    def __read_txt_result(self,file_path,data_package=CDataPack()):
        counter=1
        with open(file_path,'r') as f:
            for line in f:
                if len(line.strip())==0:
                    continue
                else:
                    protein1=line.strip().split('\t')[0]
                    site1=int(line.strip().split('\t')[1])
                    protein2 = line.strip().split('\t')[2]
                    site2 = int(line.strip().split('\t')[3])
                    temp_ppi = Cppi(protein1, site1, protein2, site2)
                    if temp_ppi not in data_package.my_CxResult.CFLOW6_PPI:
                        data_package.my_CxResult.CFLOW6_PPI.append(temp_ppi)
                        data_package.my_CxResult.CFLOW6_SPECTRA.append([counter])
                        data_package.my_CxResult.CFLOW6_SCORE.append([1])
                    else:
                        ppi_index = data_package.my_CxResult.CFLOW6_PPI.index(temp_ppi)
                        data_package.my_CxResult.CFLOW6_SPECTRA[ppi_index].append(counter)
                        data_package.my_CxResult.CFLOW6_SCORE[ppi_index].append(1)
                    counter+=1

    def __read_spotlink_psm_result(self,file_path, data_package=CDataPack()):
        with open(file_path, 'r') as f:
            csv_f = csv.reader(f)
            for line in csv_f:
                if len(line) == 0 or line[0] == 'Order':
                    continue
                else:
                    spectra_title = line[1]
                    score = float(line[10])
                    ppis = line[5]
                    for ppi_item in ppis.split('/'):
                        if '-' in ppi_item:
                            pro1 = eval(line[8])[0]
                            pro2 = eval(line[9])[0]
                            if 'DECOY' in pro1 or 'DECOY' in pro2:
                                continue
                            seq1 = ppi_item.split('(')[0]
                            local_s1 = int(ppi_item.split('-')[0].split('(')[1][:-1])
                            seq2 = ppi_item.split('-')[1].split('(')[0]
                            local_s2 = int(ppi_item.split('-')[1].split('(')[1][:-1])
                            s1 = data_package.my_ProteinIndex.CFLOW6_SEQUENCE_INDEX[pro1].index(seq1) + local_s1
                            s2 = data_package.my_ProteinIndex.CFLOW6_SEQUENCE_INDEX[pro2].index(seq2) + local_s2
                            temp_ppi = Cppi(pro1, s1, pro2, s2)
                            if temp_ppi not in data_package.my_CxResult.CFLOW6_PPI:
                                data_package.my_CxResult.CFLOW6_PPI.append(temp_ppi)
                                data_package.my_CxResult.CFLOW6_SPECTRA.append([spectra_title])
                                data_package.my_CxResult.CFLOW6_SCORE.append([score])
                            else:
                                ppi_index = data_package.my_CxResult.CFLOW6_PPI.index(temp_ppi)
                                data_package.my_CxResult.CFLOW6_SPECTRA[ppi_index].append(spectra_title)
                                data_package.my_CxResult.CFLOW6_SCORE[ppi_index].append(score)

    def __read_plink2_psm_result(self, file_path, data_package=CDataPack()):
        with open(file_path, 'r') as f:
            csv_f = csv.reader(f)
            for line in csv_f:
                if len(line) == 0 or line[0] == 'Order':
                    continue
                else:
                    spectra_title = line[1]
                    score=float(line[10])
                    ppis = line[13]
                    for ppi_item in ppis.split('/'):
                        if '-' in ppi_item:
                            pro1 = ppi_item.split(')-')[0].split('|')[1]
                            pro2 = ppi_item.split(')-')[1].split('|')[1]
                            s1 = int(ppi_item.split(')-')[0].split('(')[1])
                            s2 = int(ppi_item.split(')-')[1].split('(')[1][:-1])
                            temp_ppi = Cppi(pro1, s1, pro2, s2)
                            if temp_ppi not in data_package.my_CxResult.CFLOW6_PPI:
                                data_package.my_CxResult.CFLOW6_PPI.append(temp_ppi)
                                data_package.my_CxResult.CFLOW6_SPECTRA.append([spectra_title])
                                data_package.my_CxResult.CFLOW6_SCORE.append([score])
                            else:
                                ppi_index = data_package.my_CxResult.CFLOW6_PPI.index(temp_ppi)
                                data_package.my_CxResult.CFLOW6_SPECTRA[ppi_index].append(spectra_title)
                                data_package.my_CxResult.CFLOW6_SCORE[ppi_index].append(score)

    def __read_plink2_loop_result(self,file_path,data_package=CDataPack()):
        with open(file_path, 'r') as f:
            csv_f = csv.reader(f)
            for line in csv_f:
                if len(line) == 0 or line[0] == 'Order':
                    continue
                else:
                    spectra_title = line[1]
                    score=float(line[10])
                    ppis = line[13]
                    for ppi_item in ppis.split('/'):
                        if '|' in ppi_item:
                            pro1 = ppi_item.split('(')[0].split('|')[1]
                            pro2 = ppi_item.split('(')[0].split('|')[1]
                            s1 = int(ppi_item.split('(')[1][:-1])
                            s2 = int(ppi_item.split('(')[2][:-1])
                            temp_ppi = Cppi(pro1, s1, pro2, s2)
                            if temp_ppi not in data_package.my_CxResult.CFLOW6_PPI:
                                data_package.my_CxResult.CFLOW6_PPI.append(temp_ppi)
                                data_package.my_CxResult.CFLOW6_SPECTRA.append([spectra_title])
                                data_package.my_CxResult.CFLOW6_SCORE.append([score])
                            else:
                                ppi_index = data_package.my_CxResult.CFLOW6_PPI.index(temp_ppi)
                                data_package.my_CxResult.CFLOW6_SPECTRA[ppi_index].append(spectra_title)
                                data_package.my_CxResult.CFLOW6_SCORE[ppi_index].append(score)

    def __read_xlinkx_txt_result(self,file_path,data_package=CDataPack()):
        score_scan_list=[]
        with open(file_path,'r') as f:
            for line in f:
                if line.strip() and 'Checked' in line:
                    continue
                else:
                    score_scan_list.append(float(line.split('\t')[10][1:-1]))
        max_score,min_score=max(score_scan_list),min(score_scan_list)
        with open(file_path,'r') as f:
            for line in f:
                if line.strip() and 'Checked' in line:
                    continue
                elif line.split('\t')[19][1:-1]=="True":
                    continue
                else:
                    spectra_title=line.split('\t')[40]+line.split('\t')[35]
                    score=1-(float(line.split('\t')[10][1:-1])/max_score)
                    pro1=line.split('\t')[21][1:-1].split('-')[0]
                    pro2=line.split('\t')[21][1:-1].split('-')[1]
                    seq1=line.split('\t')[22][1:-1]
                    local_s1=int(line.split('\t')[24][1:-1])
                    seq2=line.split('\t')[28][1:-1]
                    local_s2=int(line.split('\t')[30][1:-1])
                    try:
                        s1=data_package.my_ProteinIndex.CFLOW6_SEQUENCE_INDEX[pro1].index(seq1)+local_s1
                        s2=data_package.my_ProteinIndex.CFLOW6_SEQUENCE_INDEX[pro2].index(seq2)+local_s2
                        temp_ppi = Cppi(pro1, s1, pro2, s2)
                    except:
                        s1 = data_package.my_ProteinIndex.CFLOW6_SEQUENCE_INDEX[pro2].index(seq1) + local_s1
                        s2 = data_package.my_ProteinIndex.CFLOW6_SEQUENCE_INDEX[pro1].index(seq2) + local_s2
                        temp_ppi = Cppi(pro1, s2, pro2, s1)
                    if temp_ppi not in data_package.my_CxResult.CFLOW6_PPI:
                        data_package.my_CxResult.CFLOW6_PPI.append(temp_ppi)
                        data_package.my_CxResult.CFLOW6_SPECTRA.append([spectra_title])
                        data_package.my_CxResult.CFLOW6_SCORE.append([score])
                    else:
                        ppi_index = data_package.my_CxResult.CFLOW6_PPI.index(temp_ppi)
                        data_package.my_CxResult.CFLOW6_SPECTRA[ppi_index].append(spectra_title)
                        data_package.my_CxResult.CFLOW6_SCORE[ppi_index].append(score)


"""
下载相关的函数
"""

class CFunctionDownloadProteinInfo():
    def __init__(self, database_name='', base_url='', protein_name='', file_type=''):
        self.databaseName = database_name
        self.proteinName = protein_name
        self.fileType = file_type
        self.url = base_url + protein_name + file_type

    def get_protein_info(self, protein_url=''):
        '''
        Get protein infromation from input URL or self.url
        :param protein_url:
        :return:
        '''
        if protein_url == '':
            protein_url = self.url
        else:
            if protein_url[0:4] in ['http', 'HTTP']:
                # It the URL is a HTTP request
                try:
                    response = requests.get(protein_url)
                except requests.exceptions.HTTPError:
                    return RESPONSE_RESULT['Fail']
                else:
                    if response.status_code != 200:
                        return RESPONSE_RESULT['Fail']
                    return response.content
            elif protein_url[0:3] in ['ftp', 'FTP']:
                # If the URL is a FTP request
                try:
                    response = urllib.request.urlopen(protein_url)
                except urllib.error.URLError:
                    return RESPONSE_RESULT['Fail']
                else:
                    return response

    def save_protein_info(self, response_content, protein_name, method,folder):
        '''
        Save protein information from HTTP request or FTP request
        :param response_content:
        :param protein_name:
        :return:
        '''
        file_path = folder + '/' + protein_name + self.fileType
        if method == 'overall':
            with open(file_path, 'wb') as fp:
                fp.write(response_content)
        elif method == 'byline':
            with open(file_path, 'wb') as fp:
                for line in response_content:
                    fp.write(line)

    def download_protein_info(self, protein_url, protein_name,folder):
        protein_info = self.get_protein_info(protein_url)
        if protein_info == 0:
            return (RESPONSE_RESULT['Fail'], protein_name)
        else:
            if protein_url[0:4] in ['http', 'HTTP']:
                self.save_protein_info(protein_info, protein_name, method='overall',folder=folder)
            elif protein_url[0:3] in ['ftp', 'FTP']:
                self.save_protein_info(protein_info, protein_name, method='byline',folder=folder)
            return (RESPONSE_RESULT['Success'], protein_name)


class CFunctionDownloadProteinsInfo(CFunctionDownloadProteinInfo):
    def __init__(self, database_name='', base_url='', proteins_list=[], file_type='', concurrent_threads=5):
        self.databaseName = database_name
        self.proteinsList = proteins_list
        self.fileType = file_type
        self.baseUrl = base_url
        self.concurrentThreads = concurrent_threads

    def download_proteins_info(self, input_list=[],folder=''):
        success_proteins = []
        failure_proteins = []
        if input_list == []:
            input_list = self.proteinsList
        with cf.ThreadPoolExecutor(max_workers=min(self.concurrentThreads, MAX_CONCURRENT_NUMBER)) as executor:
            queue = {}
            for protein in input_list:
                if self.databaseName=='PDB':
                    protein_url = self.baseUrl+protein[1:-1]+'/'+ protein + self.fileType
                else:
                    protein_url = self.baseUrl + protein + self.fileType
                future = executor.submit(self.download_protein_info, protein_url, protein,folder)
                queue[future] = protein
            future_list = cf.as_completed(queue)
            for future_result in future_list:
                if future_result.result()[0] == 1:
                    success_proteins.append(future_result.result()[1])
                    print('[ComMap] Get {0} from {1} successfully!'.format(future_result.result()[1], self.databaseName))
                else:
                    failure_proteins.append(future_result.result()[1])
                    print('[ComMap] Get {0} from {1} failed!'.format(future_result.result()[1], self.databaseName))
        return success_proteins, failure_proteins

    def robust_download_protein_info(self, retry_time=3, input_list=[],folder=''):
        all_success_proteins = []
        if input_list == []:
            all_failure_proteins = self.proteinsList
        else:
            all_failure_proteins = input_list
        for i in range(retry_time):
            success_proteins, failure_proteins = self.download_proteins_info(input_list=all_failure_proteins,folder=folder)
            all_success_proteins += success_proteins
            all_failure_proteins = failure_proteins
            if all_failure_proteins == []:
                break
        return all_success_proteins, all_failure_proteins


class CFunctionDownloadUniprotInfo(CFunctionDownloadProteinInfo):
    def __init__(self, file_type='xml', protein_name='',folder=''):
        self.databaseName = 'Uniprot'
        self.proteinName = protein_name
        if file_type == 'xml':
            self.fileType = UNIPROT_FILE_TYPE
        self.url = UNIPROT_BASE_URL + protein_name + UNIPROT_FILE_TYPE
        self.folder=folder

    def download_protein_info(self):
        status, protein = super().download_protein_info(self.url, self.proteinName,self.folder)
        return status, protein


class CFunctionDownloadUniprotInfos(CFunctionDownloadProteinsInfo):
    def __init__(self, file_type='xml', proteins_list=[], concurrent_threads=5,folder=''):
        self.databaseName = 'Uniprot'
        self.proteinsList = proteins_list
        if file_type == 'xml':
            self.fileType = UNIPROT_FILE_TYPE
        self.baseUrl = UNIPROT_BASE_URL
        self.concurrentThreads = concurrent_threads
        self.folder=folder


class CFunctionDownloadPDBInfo(CFunctionDownloadProteinInfo):
    def __init__(self, file_type, protein_name='',folder=''):
        self.databaseName = 'PDB'
        self.proteinName = protein_name.lower()
        self.folder=folder
        if file_type == 'pdb':
            self.fileType = PDBJ_PDB_FILE_TYPE
            self.url = PDBJ_PDB_FORMAT_BASE_URL +r'/'+ protein_name.lower()[1:-1]+ r'/'+ 'pdb'+protein_name.lower()+self.fileType
        elif file_type == 'xml':
            self.fileType = PDBJ_XML_FILE_TYPE
            self.url = PDBJ_XML_FORMAT_BASE_URL + protein_name.lower() + self.fileType
        elif file_type=='cif':
            self.fileType=PDBJ_CIF_FILE_TYPE
            self.url=PDBJ_CIF_FORMAT_BASE_URL+r'/'+ protein_name.lower()[1:-1]+ r'/'+protein_name.lower()+self.fileType

    def download_protein_info(self):
        status, protein = super().download_protein_info(protein_url=self.url,
                                                        protein_name=self.proteinName,folder=self.folder)
        return status, protein


class CFunctionDownloadPDBInfos(CFunctionDownloadProteinsInfo):
    def __init__(self, file_type, proteins_list=[], concurrent_threads=5,folder=''):
        self.databaseName = 'PDB'
        self.proteinsList = [p.lower() for p in proteins_list]
        if file_type == 'pdb':
            self.fileType = PDBJ_PDB_FILE_TYPE
            self.baseUrl = PDBJ_PDB_FORMAT_BASE_URL
        elif file_type == 'xml':
            self.fileType = PDBJ_XML_FILE_TYPE
            self.baseUrl = PDBJ_XML_FORMAT_BASE_URL
        elif file_type == 'cif':
            self.fileType=PDBJ_CIF_FILE_TYPE
            self.baseUrl=PDBJ_CIF_FORMAT_BASE_URL
        self.concurrentThreads = concurrent_threads
        self.folder=folder

"""
Uniprot和PDB文件相关的函数
"""

class CFunctionParseUniprot:

    def buildUniprotPDBMap(self,uniprot_ids,data_package=CDataPack()):
        uniprot_files_path = [data_package.my_config.C60_TYPE01_DATA_PATH+r'/' + i + '.xml.gz' for i in uniprot_ids]
        htmlList = [etree.parse(filepath, etree.HTMLParser()) for filepath in uniprot_files_path]
        for protein,html in zip(uniprot_ids, htmlList):
            data_package.my_config.C63_PROTEINS_PDB_MAP[protein]=self.read_pdb_information(html)
        for protein in data_package.my_config.C63_PROTEINS_PDB_MAP:
            for pdb_name in data_package.my_config.C63_PROTEINS_PDB_MAP[protein].index:
                if pdb_name in data_package.my_config.C63_PDB_PROTEINS_MAP:
                    data_package.my_config.C63_PDB_PROTEINS_MAP[pdb_name].append(protein)
                else:
                    data_package.my_config.C63_PDB_PROTEINS_MAP[pdb_name]=[protein]

    def read_pdb_information(self, html):
        try:
            pdb_structure_list = html.xpath('//entry/dbreference[@type="PDB"]/@id')
            xpath = '//entry/dbreference[@type="PDB" and @id="{pdb_name}"]'
            method_dic = {}
            resolution_dic = {}
            chains_dic = {}
            for id in pdb_structure_list:
                new_xpath = xpath.format(pdb_name=id)
                id_method = html.xpath(new_xpath + '/property[@type="method"]/@value')
                id_resolution = html.xpath(new_xpath + '/property[@type="resolution"]/@value')
                id_chains = html.xpath(new_xpath + '/property[@type="chains"]/@value')
                if id_method != []:
                    method_dic[id] = id_method[0]
                if id_resolution != []:
                    resolution_dic[id] = id_resolution[0]
                if id_chains != []:
                    chains_dic[id] = id_chains[0].replace(',', ';')
            pdb_list = pd.DataFrame({'method': method_dic, 'resolution': resolution_dic, 'chains': chains_dic})
        except AssertionError:
            pdb_list = pd.DataFrame({'method': {}, 'resolution': {}, 'chains': {}})
        return pdb_list


class CFunctionParseCIF:

    def parseCifAtomInformation(self,pdb_name,folder):
        item_id_list = []
        result_list = []
        cif_file_path = folder + '\{}.cif.gz'.format(pdb_name.lower())
        atom_chunk=False
        with gzip.open(cif_file_path, 'rt') as f:
            for line in f:
                if line.startswith("_atom_site."):
                    item_id_list.append(line.strip().split('.')[1])
                    atom_chunk=True
                elif line.startswith("ATOM") and atom_chunk:
                    result_split=line.strip().split()
                    if result_split[0]=="ATOM":
                        result_list.append(result_split)
        return item_id_list, result_list

    def parseCIFSequenceFromAtom(self, pdb_name,folder):
        item_id_list, result_list = self.parseCifAtomInformation(pdb_name,folder)
        if 'pdbx_PDB_model_num' not in item_id_list or 'label_asym_id' not in item_id_list or 'label_atom_id' not in\
            item_id_list or 'group_PDB' not in item_id_list or 'label_comp_id' not in item_id_list or 'auth_seq_id' not\
                in item_id_list:
            print(f'\n[ComMap] Failed parsing structure {pdb_name}.')
            return [], []
        idx_model_num=item_id_list.index('pdbx_PDB_model_num')
        idx_asym_id=item_id_list.index('label_asym_id')
        idx_atom_id=item_id_list.index('label_atom_id')
        idx_group_pdb=item_id_list.index('group_PDB')
        idx_comp_id=item_id_list.index('label_comp_id')
        idx_seq_id=item_id_list.index('auth_seq_id')
        last_model,last_chain = result_list[0][idx_model_num],result_list[0][idx_asym_id]
        model_seq_list = []
        model_ord_list=[]
        chainseq_list = []
        chainord_list=[]
        sequence = ''
        order_list=[]
        for index in range(len(result_list)):
            if result_list[index][idx_atom_id] != 'CA' or result_list[index][idx_group_pdb] != 'ATOM':
                continue
            if result_list[index][idx_model_num] != last_model:
                chainseq_list.append((sequence,last_chain))
                model_seq_list.append(chainseq_list)
                chainseq_list = []
                chainord_list.append(order_list)
                model_ord_list.append(chainord_list)
                chainord_list=[]
                last_model = result_list[index][idx_model_num]
                sequence = ''
                order_list=[]
                last_chain = result_list[index][idx_asym_id]
                sequence+=AMINO_ACIDS_DIC[result_list[index][idx_comp_id]]
                order_list.append(int(result_list[index][idx_seq_id]))
            else:
                if result_list[index][idx_asym_id] != last_chain:
                    chainseq_list.append((sequence,last_chain))
                    chainord_list.append(order_list)
                    sequence = ''
                    order_list=[]
                    last_chain = result_list[index][idx_asym_id]
                sequence+=AMINO_ACIDS_DIC[result_list[index][idx_comp_id]]
                order_list.append(int(result_list[index][idx_seq_id]))
        chainseq_list.append((sequence,last_chain))
        chainord_list.append(order_list)
        model_seq_list.append(chainseq_list)
        model_ord_list.append(chainord_list)
        return model_seq_list, model_ord_list

    def parseCIFSequence(self, pdb_name,folder):
        cif_file_path = folder+'\{}.cif.gz'.format(pdb_name.lower())
        with gzip.open(cif_file_path, 'rt') as f:
            result_parser=[]
            pRd = PdbxReader(f)
            try:
                pRd.read(result_parser)
            except:
                print(f'\n[ComMap] Failed parsing structure {pdb_name}.')
                return [],[]
        last_model = result_parser[0].getObj('atom_site').getValue('pdbx_PDB_model_num', 0)
        last_chain = result_parser[0].getObj('atom_site').getValue('label_asym_id', 0)
        model_seq_list = []
        model_ord_list=[]
        chainseq_list = []
        chainord_list=[]
        sequence = ''
        order_list=[]
        for index in range(result_parser[0].getObj('atom_site').getRowCount()):
            if result_parser[0].getObj('atom_site').getValue('label_atom_id', index) != 'CA' or\
                    result_parser[0].getObj('atom_site').getValue('group_PDB', index) != 'ATOM':
                continue
            if result_parser[0].getObj('atom_site').getValue('pdbx_PDB_model_num', index) != last_model:
                chainseq_list.append((sequence,last_chain))
                model_seq_list.append(chainseq_list)
                chainseq_list = []
                chainord_list.append(order_list)
                model_ord_list.append(chainord_list)
                chainord_list=[]
                last_model = result_parser[0].getObj('atom_site').getValue('pdbx_PDB_model_num', index)
                sequence = ''
                order_list=[]
                last_chain = result_parser[0].getObj('atom_site').getValue('label_asym_id', index)
                sequence+=AMINO_ACIDS_DIC[result_parser[0].getObj('atom_site').getValue('label_comp_id', index)]
                order_list.append(int(result_parser[0].getObj('atom_site').getValue('auth_seq_id',index)))
            else:
                if result_parser[0].getObj('atom_site').getValue('label_asym_id', index) != last_chain:
                    chainseq_list.append((sequence,last_chain))
                    chainord_list.append(order_list)
                    sequence = ''
                    order_list=[]
                    last_chain = result_parser[0].getObj('atom_site').getValue('label_asym_id', index)
                sequence+=AMINO_ACIDS_DIC[result_parser[0].getObj('atom_site').getValue('label_comp_id', index)]
                order_list.append(int(result_parser[0].getObj('atom_site').getValue('auth_seq_id', index)))
        chainseq_list.append((sequence,last_chain))
        chainord_list.append(order_list)
        model_seq_list.append(chainseq_list)
        model_ord_list.append(chainord_list)
        return model_seq_list, model_ord_list

    def parseCIFCoordinateFromAtom(self, pdb_name,folder):
        item_id_list, result_list = self.parseCifAtomInformation(pdb_name, folder)
        if 'pdbx_PDB_model_num' not in item_id_list or 'label_asym_id' not in item_id_list or 'label_atom_id' not in \
                item_id_list or 'group_PDB' not in item_id_list or 'label_comp_id' not in item_id_list or 'auth_seq_id' not \
                in item_id_list or 'Cartn_x' not in item_id_list or 'Cartn_y' not in item_id_list or 'Cartn_z' not\
                in item_id_list or 'auth_asym_id' not in item_id_list:
            print(f'\n[ComMap] Failed parsing structure {pdb_name}.')
            return [], []
        idx_model_num = item_id_list.index('pdbx_PDB_model_num')
        idx_asym_id = item_id_list.index('label_asym_id')
        idx_atom_id = item_id_list.index('label_atom_id')
        idx_group_pdb = item_id_list.index('group_PDB')
        idx_comp_id = item_id_list.index('label_comp_id')
        idx_seq_id = item_id_list.index('auth_seq_id')
        idx_x=item_id_list.index('Cartn_x')
        idx_y=item_id_list.index('Cartn_y')
        idx_z=item_id_list.index('Cartn_z')
        idx_auth_asym_id=item_id_list.index('auth_asym_id')
        last_model, last_chain = result_list[0][idx_model_num], result_list[0][idx_asym_id]
        model_list = []
        coordinate_list = []
        chain_list = {}
        for index in range(len(result_list)):
            if result_list[index][idx_atom_id] != 'CA':
                continue
            if result_list[index][idx_model_num] != last_model:
                coordinate_list.append(chain_list)
                model_list.append(coordinate_list)
                coordinate_list = []
                last_model = result_list[index][idx_model_num]
                chain_list = {}
                last_chain = result_list[index][idx_asym_id]
                chain_list[int(result_list[index][idx_seq_id])]=(result_list[index][idx_x],result_list[index][idx_y],
                                   result_list[index][idx_z],result_list[index][idx_auth_asym_id])
            else:
                if result_list[index][idx_asym_id] != last_chain:
                    coordinate_list.append(chain_list)
                    chain_list = {}
                    last_chain = result_list[index][idx_asym_id]
                chain_list[int(result_list[index][idx_seq_id])] = (result_list[index][idx_x], result_list[index][idx_y],
                                                                   result_list[index][idx_z],result_list[index][idx_auth_asym_id])
        coordinate_list.append(chain_list)
        model_list.append(coordinate_list)
        return model_list

    def parseCIFCoordinate(self, pdb_name,folder):
        cif_file_path = folder + '\{}.cif.gz'.format(pdb_name.lower())
        with gzip.open(cif_file_path, 'rt') as f:
            result_parser = []
            pRd = PdbxReader(f)
            try:
                pRd.read(result_parser)
            except:
                print(f'\n[ComMap] Failed parsing structure {pdb_name}.')
                return []
        last_model = result_parser[0].getObj('atom_site').getValue('pdbx_PDB_model_num', 0)
        last_chain = result_parser[0].getObj('atom_site').getValue('label_asym_id', 0)
        model_list = []
        coordinate_list = []
        chain_list = {}
        for index in range(result_parser[0].getObj('atom_site').getRowCount()):
            if result_parser[0].getObj('atom_site').getValue('label_atom_id', index) != 'CA':
                continue
            if result_parser[0].getObj('atom_site').getValue('pdbx_PDB_model_num', index) != last_model:
                coordinate_list.append(chain_list)
                model_list.append(coordinate_list)
                coordinate_list = []
                last_model = result_parser[0].getObj('atom_site').getValue('pdbx_PDB_model_num', index)
                chain_list = {}
                last_chain = result_parser[0].getObj('atom_site').getValue('label_asym_id', index)
                chain_list[int(result_parser[0].getObj('atom_site').getValue('auth_seq_id',index))]=(result_parser[0].getObj('atom_site').getValue('Cartn_x', index),
                                   result_parser[0].getObj('atom_site').getValue('Cartn_y', index),
                                   result_parser[0].getObj('atom_site').getValue('Cartn_z', index))
            else:
                if result_parser[0].getObj('atom_site').getValue('label_asym_id', index) != last_chain:
                    coordinate_list.append(chain_list)
                    chain_list = {}
                    last_chain = result_parser[0].getObj('atom_site').getValue('label_asym_id', index)
                chain_list[int(result_parser[0].getObj('atom_site').getValue('auth_seq_id',index))]=(result_parser[0].getObj('atom_site').getValue('Cartn_x', index),
                                   result_parser[0].getObj('atom_site').getValue('Cartn_y', index),
                                   result_parser[0].getObj('atom_site').getValue('Cartn_z', index),
                                   result_parser[0].getObj('atom_site').getValue('auth_asym_id',index))
        coordinate_list.append(chain_list)
        model_list.append(coordinate_list)
        return model_list


class CFunctionParsePDB:

    def parsePDBSequence(self, pdb_name,folder):
        pdb_file_path = folder+'\{}.ent.gz'.format(pdb_name.lower())
        with gzip.open(pdb_file_path, 'rt') as f:
            file_lines = f.readlines()
        number_of_model = 0
        for item in file_lines:
            if item.startswith('MODEL'):
                number_of_model += 1
        if number_of_model == 0:
            pdb_sequence, pdb_sequence_order = [[]], [[]]
            last_domain, last_sequence, last_order = '', '', []
            for line in file_lines:
                if line.startswith('ATOM'):
                    amino_acid = line[17:20].strip()
                    atom = line[13:16].strip()
                    domain = line[21]
                    order = int(line[22:26])
                    if atom != 'CA' or amino_acid == 'UNK':
                        continue
                    one_letter_amino_acid = AMINO_ACIDS_DIC[amino_acid]
                    if last_domain != domain:
                        if last_sequence != '' and last_domain != '' and last_order != []:
                            pdb_sequence[0].append((last_sequence, last_domain))
                            pdb_sequence_order[0].append(last_order)
                        last_sequence = one_letter_amino_acid
                        last_domain = domain
                        last_order = [order]
                    else:
                        last_sequence += one_letter_amino_acid
                        last_order.append(order)
            if last_sequence != '' and last_domain != '' and last_order != []:
                pdb_sequence[0].append((last_sequence, last_domain))
                pdb_sequence_order[0].append(last_order)
        else:
            pdb_sequence = [[] for _ in range(number_of_model)]
            pdb_sequence_order = [[] for _ in range(number_of_model)]
            last_domain, last_sequence, last_order = '', '', []
            model_number = -1
            for line in file_lines:
                if line.startswith('MODEL'):
                    if last_sequence != '' and model_number != -1 and last_domain != '' and last_order != []:
                        pdb_sequence[model_number].append((last_sequence, last_domain))
                        pdb_sequence_order[model_number].append(last_order)
                        last_domain, last_sequence, last_order = '', '', []
                    model_number += 1
                if line.startswith('ATOM'):
                    amino_acid = line[17:20].strip()
                    atom = line[13:16].strip()
                    domain = line[21]
                    order = int(line[22:26])
                    if atom != 'CA' or amino_acid == 'UNK':
                        continue
                    one_letter_amino_acid = AMINO_ACIDS_DIC[amino_acid]
                    if last_domain != domain:
                        if last_sequence != '' and last_domain != '' and last_order != []:
                            pdb_sequence[model_number].append((last_sequence, last_domain))
                            pdb_sequence_order[model_number].append(last_order)
                        last_sequence = one_letter_amino_acid
                        last_domain = domain
                        last_order = [order]
                    else:
                        last_sequence += one_letter_amino_acid
                        last_order.append(order)
            if last_sequence != '' and last_domain != '' and last_order != []:
                pdb_sequence[model_number].append((last_sequence, last_domain))
                pdb_sequence_order[model_number].append(last_order)
        return pdb_sequence, pdb_sequence_order

    def parsePDBCoordinate(self, pdb_name,folder):
        pdb_file_path = folder+r'/{}.ent.gz'.format(pdb_name.lower())
        with gzip.open(pdb_file_path, 'rt') as f:
            file_lines = f.readlines()
        number_of_model = 0
        for item in file_lines:
            if item.startswith('MODEL'):
                number_of_model += 1
        if number_of_model == 0:
            pdb_coordinates = [[]]
            last_domain, last_coordinates = '', []
            for line in file_lines:
                if line.startswith('ATOM'):
                    amino_acid = line[17:20].strip()
                    atom = line[13:16].strip()
                    domain = line[21]
                    if atom != 'CA' or amino_acid == 'UNK':
                        continue
                    raw_coordinates = line[30:54]
                    x_coordinate = float(raw_coordinates[:8])
                    y_coordinate = float(raw_coordinates[8:16])
                    z_coordinate = float(raw_coordinates[16:24])
                    coordinates = (float(x_coordinate), float(y_coordinate), float(z_coordinate))
                    if last_domain != domain:
                        if last_domain != '' and last_coordinates != []:
                            pdb_coordinates[0].append(last_coordinates)
                        last_domain = domain
                        last_coordinates = [coordinates]
                    else:
                        last_coordinates.append(coordinates)
            if last_domain != '' and last_coordinates != []:
                pdb_coordinates[0].append(last_coordinates)
        else:
            pdb_coordinates = [[] for _ in range(number_of_model)]
            last_domain, last_coordinates = '', []
            model_number = -1
            for line in file_lines:
                if line.startswith('MODEL'):
                    if last_domain != '' and last_coordinates != [] and model_number != -1:
                        pdb_coordinates[model_number].append(last_coordinates)
                        last_domain, last_coordinates = '', []
                    model_number += 1
                if line.startswith('ATOM'):
                    amino_acid = line[17:20].strip()
                    atom = line[13:16].strip()
                    domain = line[21]
                    if atom != 'CA' or amino_acid == 'UNK':
                        continue
                    raw_coordinates = line[30:54]
                    x_coordinate = float(raw_coordinates[:8])
                    y_coordinate = float(raw_coordinates[8:16])
                    z_coordinate = float(raw_coordinates[16:24])
                    coordinates = (float(x_coordinate), float(y_coordinate), float(z_coordinate))
                    if last_domain != domain:
                        if last_domain != '' and last_coordinates != []:
                            pdb_coordinates[model_number].append(last_coordinates)
                        last_domain = domain
                        last_coordinates = [coordinates]
                    else:
                        last_coordinates.append(coordinates)
            if last_domain != '' and last_coordinates != []:
                pdb_coordinates[model_number].append(last_coordinates)
        return pdb_coordinates


"""
与输出相关的函数
"""

class CFunctionOutPut:

    def writeTotalResult(self,data_package=CDataPack()):
        today=datetime.date.today()
        output_path = os.path.join(data_package.my_config.E1_PATH_EXPORT, f'structure_alignment_result_total.{today}.csv')
        data_package.my_config.C66_OUTPUT_DATE=today
        with open(output_path, 'w') as f:
            f.write('PPI,PDB code,Model,Protein 1 Chain,Protein 2 Chain,CA Distance\n')
            for index in range(len(data_package.my_CxResult.CFLOW6_PPI)):
                for pdb_name in data_package.my_CxResult.CFLOW6_DISTANCE[index]:
                    for item in data_package.my_CxResult.CFLOW6_DISTANCE[index][pdb_name]:
                        f.write('{},{},{},{},{},{:.3f}\n'.format(data_package.my_CxResult.CFLOW6_PPI[index],
                                                                    pdb_name, item[1], item[2], item[3], item[4],))

    def writePDBResult(self,data_package=CDataPack()):
        pdb_ppi_result={}
        for index in range(len(data_package.my_CxResult.CFLOW6_PPI)):
            for pdb_name in data_package.my_CxResult.CFLOW6_DISTANCE[index]:
                for item in data_package.my_CxResult.CFLOW6_DISTANCE[index][pdb_name]:
                    if pdb_name in pdb_ppi_result:
                        pdb_ppi_result[pdb_name].append((data_package.my_CxResult.CFLOW6_PPI[index],
                                                                pdb_name, item[1], item[2], item[3], item[4],
                                                                str('dist {},/{}/{}/{}/{}/CA,/{}/{}/{}/{}/CA'.format(
                                                                    data_package.my_CxResult.CFLOW6_PPI[index],
                                                                    pdb_name, item[2],item[7], item[5], pdb_name, item[3], item[8],
                                                                    item[6]))))
                    else:
                        pdb_ppi_result[pdb_name]=[(data_package.my_CxResult.CFLOW6_PPI[index],
                                                                pdb_name, item[1], item[2], item[3], item[4],
                                                                str('dist {},/{}/{}/{}/{}/CA,/{}/{}/{}/{}/CA'.format(
                                                                    data_package.my_CxResult.CFLOW6_PPI[index],
                                                                    pdb_name, item[2], item[7],item[5], pdb_name, item[3], item[8],
                                                                    item[6])))]
        output_path = os.path.join(data_package.my_config.E1_PATH_EXPORT, f'structure_alignment_result_PDB.{data_package.my_config.C66_OUTPUT_DATE}.csv')
        with open(output_path, 'w') as f:
            f.write('PDB code,PPI,Model,Protein 1 Chain,Protein 2 Chain,CA Distance,Dist Command\n')
            for pdb_item in pdb_ppi_result:
                f.write('{},,,,,,\n'.format(pdb_item))
                for item in pdb_ppi_result[pdb_item]:
                    f.write(',{},{},{},{},{:.3f},"{}"\n'.format(item[0],item[2], item[3], item[4], item[5],item[6]))

    def writeBinary(self,data_package=CDataPack()):
        output_path = os.path.join(data_package.my_config.E1_PATH_EXPORT,
                                   f'cache.{data_package.my_config.C66_OUTPUT_DATE}.bin')
        with open(output_path, 'wb') as f:
            pickle.dump([data_package.my_CxResult.CFLOW6_PPI,data_package.my_CxResult.CFLOW6_SPECTRA,
                         data_package.my_CxResult.CFLOW6_PPI_PDB,data_package.my_CxResult.CFLOW6_DISTANCE,
                         data_package.my_CxResult.CFLOW6_SCORE],f)

    def readBinary(self,data_package=CDataPack()):
        with open(data_package.my_config.C66_BIN_FILE, 'rb') as f:
            load_result=pickle.load(f)
            data_package.my_CxResult.CFLOW6_PPI=load_result[0]
            data_package.my_CxResult.CFLOW6_SPECTRA=load_result[1]
            data_package.my_CxResult.CFLOW6_PPI_PDB=load_result[2]
            data_package.my_CxResult.CFLOW6_DISTANCE=load_result[3]
            data_package.my_CxResult.CFLOW6_SCORE=load_result[4]


if __name__=='__main__':
    download_func=CFunctionDownloadPDBInfo('cif','4mcn')
    status,protein=download_func.download_protein_info()
    a=1