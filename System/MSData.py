class CINI:
    MASS_ELECTRON = 0.0005485799
    MASS_PROTON_MONO = 1.00727645224
    MASS_PROTON_AVERAGE = 1.0025
    MASS_H2O_MONO = 18.010565
    DICT0_ELEMENT_MASS = {}
    DICT0_ELEMENT_ABDC = {}
    DICT1_AA_COM = {}
    DICT1_AA_MASS = {'I': 113.08406, 'E': 129.04259, 'A': 71.03711, 'K': 128.09496, 'G': 57.02146, 'R': 156.10111,
                     'D': 115.02694, 'T': 101.04768, 'Q': 128.05858, 'L': 113.08406, 'V': 99.06841, 'M': 131.004049,
                     'N': 114.04293, 'F': 147.06841}
    DICT2_MOD_COM = {}
    DICT3_GLYCO_COM = {}
    DICT4_LINKER_COM = {}
    ENZYME_CUT_LIST = []
    ENZYME_CUT_IGNORE_LIST = []
    ENZYME_TERMINAL = 'N'


class CConfig:
    A3_PATH_FASTA = r''
    C_TYPE_SEARCH = 6
    C60_CALCULATION_TYPE=1
    C60_TYPE0_PDB_STRUCTURE_PATH=''
    C60_TYPE0_PROTEIN_STRUCTURE_MAP=''
    C60_TYPE01_DATA_PATH=''
    C61_INPUT_LINK_RESULT_TYPE=0
    C61_INPUT_LINK_RESULT_FILE=[]
    C62_MAX_DOWNLOAD_TRY_TIME=3
    C62_CLEAN_CALCULATE=0
    C62_UNIQUE_PROTEINS=[]
    C62_PROTIENS_TO_DOWNLOAD=[]
    C62_PROTEINS_TO_PDB=[]
    C63_PROTEINS_PDB_MAP={}
    C63_PDB_PROTEINS_MAP={}
    C63_UNIQUE_PDBS=[]
    C63_PDBS_TO_DOWNLOAD = []
    C63_SUCCESS_PDBS=[]
    C64_SEQUENCE_ALIGNMENT_E_VALUE=1e-5
    C64_SEQUENCE_ALIGNMENT_RESULT={}
    C65_LINKER_ARM_LENGTH=11.4
    C65_LINKER_MAX_ARM_LENGTH=100
    C66_OUTPUT_DATE=''
    C66_BIN_FILE=''
    P1_NUMBER_THREAD = 10
    P2_TYPE_THREAD = 0
    I0_INI_PATH_ELEMENT = ''
    I1_INI_PATH_AA = ''
    I2_INI_PATH_MOD = ''
    I3_INI_PATH_GLYCO = r''
    I4_INI_PATH_LINKER = r''
    I5_INI_PATH_ENZYME = r''
    E1_PATH_EXPORT = r''
    E2_TYPE_EXPORT = 0
    E3_FLAG_CREATE_NEM_FOLDER = 0
    E4_FLAG_EXPORT_EVIDENCE = 0
    E5_LOCAL_LOGGER = r''
    E6_LOCAL_LOG_FILE_PATH = r''

class CProteinIndex:
    CFLOW6_DESCRIPTION_INDEX = {}
    CFLOW6_SEQUENCE_INDEX = {}

class Cppi:
    def __init__(self,pro1,site1,pro2,site2):
        self.protein1=pro1
        self.protein2=pro2
        self.site1=site1
        self.site2=site2

    def __contains__(self,other_item):
        if self.protein1==other_item.protein1 and self.protein2==other_item.protein2 and self.site1==other_item.site1 and self.site2==other_item.site2:
            return True
        elif self.protein1==other_item.protein2 and self.protein2==other_item.protein1 and self.site1==other_item.site2 and self.site2==other_item.site1:
            return True
        else:
            return False

    def __eq__(self,other_item):
        if self.protein1==other_item.protein1 and self.protein2==other_item.protein2 and self.site1==other_item.site1 and self.site2==other_item.site2:
            return True
        elif self.protein1==other_item.protein2 and self.protein2==other_item.protein1 and self.site1==other_item.site2 and self.site2==other_item.site1:
            return True
        else:
            return False

    def __repr__(self):
        return '{}({})-{}({})'.format(self.protein1,self.site1,self.protein2,self.site2)

class CCxResult:
    CFLOW6_PPI=[]
    CFLOW6_SPECTRA=[]
    CFLOW6_PPI_PDB=[]
    CFLOW6_DISTANCE=[]
    CFLOW6_SCORE=[]


class CDataPack:
    my_config = CConfig()
    my_ini = CINI()
    my_ProteinIndex=CProteinIndex()
    my_CxResult=CCxResult()
